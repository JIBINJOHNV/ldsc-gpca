"""LDSC result validation, compilation and reuse checks."""
import os
import math
import re
import pandas as pd
from .ldsc_export import BASE_RESULT_COLUMNS, HERITABILITY_COLUMNS, write_results_csv


def _trait_name(path):
    return os.path.basename(path).replace('.sumstats.gz', '')


def _detailed_correlations(lines, log):
    """Read rg/SE/Z/P at the precision printed in each named LDSC block.

    The final pandas display table can round a positive self-rg SE to zero.
    Match blocks by trait pair, never by their position in the summary table.
    These are logged estimates, not a reconstruction from rounded rg/Z.
    """
    details = {}
    reference = target = None
    active = False
    record = {}

    def finish_block():
        if not record:
            return  # Failed computations are rejected by table validation below.
        if reference is None or target is None or set(record) != {'rg', 'se', 'z', 'p'}:
            raise RuntimeError(f'Incomplete detailed LDSC correlation block in {log}')
        pair = (_trait_name(reference), _trait_name(target))
        if pair in details:
            raise RuntimeError(f'Duplicate detailed LDSC correlation block for {pair} in {log}')
        details[pair] = record.copy()

    for raw_line in lines:
        line = raw_line.strip()
        if line.startswith('Summary of Genetic Correlation Results'):
            finish_block()
            break
        if line.startswith('Computing rg for phenotype '):
            finish_block()
            target, record, active = None, {}, True
        elif line.startswith('Reading summary statistics from '):
            path = line.removeprefix('Reading summary statistics from ').removesuffix(' ...')
            if active:
                target = path
            elif reference is None:
                reference = path
        elif active and line.startswith('Genetic Correlation:'):
            match = re.fullmatch(r'Genetic Correlation:\s*(\S+)\s*\(([^)]+)\)', line)
            if match is None:
                raise RuntimeError(f'Malformed detailed LDSC correlation in {log}: {line}')
            record.update(rg=match[1], se=match[2].strip())
        elif active and line.startswith('Z-score:'):
            record['z'] = line.removeprefix('Z-score:').strip()
        elif active and line.startswith('P:'):
            record['p'] = line.removeprefix('P:').strip()
    else:
        finish_block()

    for pair, estimates in details.items():
        try:
            for value in estimates.values():
                float(value)
        except ValueError as error:
            raise RuntimeError(f'Non-numeric detailed LDSC result for {pair} in {log}') from error
    return details


def _restore_correlation_precision(rows, lines, log):
    """Prefer complete detailed results, with a safe legacy table-only fallback."""
    details = _detailed_correlations(lines, log)
    table_pairs = {(_trait_name(row['p1']), _trait_name(row['p2'])) for row in rows}
    unexpected = set(details) - table_pairs
    if unexpected:
        raise RuntimeError(f'Detailed LDSC pairs absent from summary table in {log}: {sorted(unexpected)}')
    for row in rows:
        pair = (_trait_name(row['p1']), _trait_name(row['p2']))
        if pair in details:
            row.update(details[pair])
        elif details and pd.notna(pd.to_numeric(row['rg'], errors='coerce')):
            raise RuntimeError(f'Missing detailed LDSC correlation for {pair} in {log}')
        se = pd.to_numeric(row['se'], errors='coerce')
        if not math.isfinite(se) or se <= 0:
            raise RuntimeError(
                f'Non-positive or non-finite LDSC rg SE for {pair} in {log}. '
                'A summary-table zero may be rounding; retain complete detailed logs '
                'to recover a reported positive SE. No SE was imputed.'
            )
    return rows


def _read_legacy_log(log):
    """Explicit legacy recovery only; managed regressions never use this path."""
    try:
        with open(log, 'r') as handle:
            lines = handle.readlines()
    except OSError as error:
        raise RuntimeError(f'Cannot read LDSC log {log}: {error}') from error
    rows = []
    for index, line in enumerate(lines):
        if "gcov_int_se" not in line:
            continue
        headers = line.split()
        for table_line in lines[index + 1:]:
            parts = table_line.split()
            if not parts or "Summary" in table_line:
                break
            if parts == headers:
                continue
            if len(parts) != len(headers):
                raise RuntimeError(f'Malformed LDSC result row in {log}: {table_line.strip()}')
            rows.append(dict(zip(headers, parts)))
        break
    if not rows:
        raise RuntimeError(f'No correlation table rows found in LDSC log: {log}')
    return pd.DataFrame(_restore_correlation_precision(rows, lines, log))


def _read_numerical_csv(path):
    """Read each machine-readable export once, without opening a readable log."""
    try:
        frame = pd.read_csv(path, float_precision='round_trip',
                            dtype={'p1': str, 'p2': str})
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise RuntimeError(
            f'Cannot read numerical LDSC CSV {path}: {error}. '
            'No readable-log fallback was attempted.'
        ) from error
    missing = set(BASE_RESULT_COLUMNS) - set(frame.columns)
    h2_columns = [pair for pair in HERITABILITY_COLUMNS if set(pair).issubset(frame.columns)]
    if missing or not h2_columns:
        raise RuntimeError(
            f'Missing LDSC result columns in {path}: {sorted(missing)}; '
            'require a complete h2_obs/h2_obs_se or h2_liab/h2_liab_se pair.'
        )
    if frame.empty:
        raise RuntimeError(f'No correlation results found in numerical LDSC CSV: {path}')
    for column in ['p1', 'p2']:
        if frame[column].isna().any() or frame[column].str.strip().eq('').any():
            raise RuntimeError(f'Missing LDSC trait identifier {column} in {path}')
    numeric_columns = [column for column in BASE_RESULT_COLUMNS if column not in ('p1', 'p2')]
    numeric_columns.extend(column for pair in h2_columns for column in pair)
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors='coerce')
        bad = ~frame[column].map(math.isfinite)
        if bad.any():
            pairs = list(zip(frame.loc[bad, 'p1'], frame.loc[bad, 'p2']))
            raise RuntimeError(f'Non-finite LDSC {column} in {path} for pairs: {pairs}')
    se_columns = ['se', 'h2_int_se', 'gcov_int_se'] + [pair[1] for pair in h2_columns]
    for column in se_columns:
        if frame[column].le(0).any():
            raise RuntimeError(
                f'Non-positive LDSC {column} in numerical CSV {path}; '
                'the original numerical result was exported. No SE was imputed.'
            )
    if (~frame['p'].between(0, 1)).any():
        raise RuntimeError(f'LDSC p-values outside [0,1] in {path}')
    return frame


def compile_results(output_folder, log_files, trait_metadata):
    """Compile fresh numerical CSVs; explicitly supplied .log paths retain legacy support.

    The historical log_files parameter name is retained for API compatibility.
    Missing/invalid CSVs fail, never silently falling back to rounded log text.
    """
    if not log_files:
        raise RuntimeError('No correlation results found in current LDSC outputs')
    print(f"  -> Reading {len(log_files)} LDSC result files...")
    frames = [(_read_legacy_log(path) if os.fspath(path).endswith('.log')
               else _read_numerical_csv(path)) for path in log_files]
    df = pd.concat(frames, ignore_index=True)
    for column in ['p1', 'p2']:
        df[column] = df[column].apply(_trait_name)
    df = df[df['p1'] != 'p1'].drop_duplicates()
    expected = {(ref, target) for ref in trait_metadata.loc[trait_metadata['ref'] == 'yes', 'gwas_name']
                for target in trait_metadata['gwas_name']}
    actual = set(zip(df['p1'], df['p2']))
    if expected != actual:
        raise RuntimeError(f'LDSC comparison mismatch: missing={sorted(expected - actual)}; unexpected={sorted(actual - expected)}')
    values = pd.to_numeric(df['rg'], errors='coerce')
    bad = ~values.map(math.isfinite)
    if bad.any():
        raise RuntimeError(f'Non-finite LDSC rg for pairs: {list(zip(df.loc[bad, "p1"], df.loc[bad, "p2"]))}')
    prevalence = trait_metadata.set_index('gwas_name')['pop_prevalence']
    # Preserve the historical scale mapping; no estimates are recomputed.
    for name in ['h2_obs', 'h2_obs_se', 'h2_liab', 'h2_liab_se']:
        if name not in df:
            df[name] = float('nan')
    no_conversion = df['p2'].map(prevalence).isna()
    for observed, liability in [('h2_obs', 'h2_liab'), ('h2_obs_se', 'h2_liab_se')]:
        df.loc[no_conversion, observed] = df.loc[no_conversion, observed].fillna(df.loc[no_conversion, liability])
        df.loc[no_conversion, liability] = float('nan')
    df['h2_scale'] = df['p2'].map(lambda name: 'NEF_unconverted' if pd.isna(prevalence[name]) else 'liability')
    write_results_csv(df, os.path.join(output_folder, 'ldsc_results.csv'))
    print(f"  -> Successfully compiled {len(df)} correlations.")


def check_saved_filters(metadata, filters, trait):
    """Reject changed filters; legacy files cannot establish filter provenance."""
    saved = metadata.get('filters')
    if saved is None:
        print(f'WARNING: {trait}: previous filter settings are unknown; rerun without --ldsc_only for verified filter provenance.')
    elif saved != filters:
        raise ValueError(f'{trait}: filter settings differ from the saved run; rerun without --ldsc_only')
