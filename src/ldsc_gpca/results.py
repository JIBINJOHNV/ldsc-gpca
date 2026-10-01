"""LDSC result validation, compilation and reuse checks."""
import os
import math
import re
import csv
import uuid
import pandas as pd
from .ldsc_export import BASE_RESULT_COLUMNS, HERITABILITY_COLUMNS, write_results_csv
from .result_qc import result_status, trait_status


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


def _read_numerical_csv(path, *, allow_failed=False):
    """Read each machine-readable export once, without opening a readable log."""
    try:
        with open(path, newline='', encoding='utf-8-sig') as handle:
            header = next(csv.reader(handle), [])
        if len(header) != len(set(header)):
            raise ValueError('Duplicate column headers')
        frame = pd.read_csv(path, float_precision='round_trip',
                            dtype={'p1': str, 'p2': str}, keep_default_na=False)
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
    # Extra annotation columns do not enter validation, deduplication or analysis.
    frame = frame[list(BASE_RESULT_COLUMNS) + [c for pair in h2_columns for c in pair]].copy()
    status = result_status(frame, path)
    if not allow_failed and status.Status.eq('failed_estimate').any():
        raise RuntimeError(f'{status.loc[status.Status.eq("failed_estimate"), "Reason"].iloc[0]} in {path}')
    frame.attrs['pair_status'] = status
    return frame


class EstimationFailure(RuntimeError):
    """Reportable numerical failure; its diagnostic outputs have been saved."""


def prepare_compilation_outputs(output_folder, result_failure_action='error'):
    os.makedirs(output_folder, exist_ok=True)
    # A failed rerun must never leave a previous success at the advertised path.
    # Retain old outputs for recovery rather than deleting them.
    for filename in ('ldsc_results.csv', 'ldsc_results_diagnostic.csv',
                     'LDSC_Pair_Status.csv', 'LDSC_Trait_Status.csv'):
        path = os.path.join(output_folder, filename)
        if os.path.exists(path):
            os.replace(path, path + '.previous-' + uuid.uuid4().hex)
    status_path = os.path.join(output_folder, 'LDSC_Compilation_Status.csv')
    write_results_csv(pd.DataFrame([{'Status': 'incomplete', 'Action': result_failure_action,
                                    'Result_File': '', 'Failed_Rows': ''}]), status_path)


def compile_results(output_folder, log_files, trait_metadata, *, result_failure_action='error'):
    if result_failure_action not in ('error', 'report'):
        raise ValueError('result_failure_action must be error or report')
    prepare_compilation_outputs(output_folder, result_failure_action)
    try:
        return _compile_results(output_folder, log_files, trait_metadata,
                                result_failure_action=result_failure_action)
    except EstimationFailure:
        raise
    except (RuntimeError, ValueError, OSError) as error:
        write_results_csv(pd.DataFrame([{'Status': 'structural_failure',
            'Action': result_failure_action, 'Result_File': '', 'Failed_Rows': '',
            'Error': str(error)}]), os.path.join(output_folder, 'LDSC_Compilation_Status.csv'))
        raise


def _compile_results(output_folder, log_files, trait_metadata, *, result_failure_action):
    """Compile fresh numerical CSVs; explicitly supplied .log paths retain legacy support.

    The historical log_files parameter name is retained for API compatibility.
    Structural failures remain fatal. Estimation failures are always audited;
    opt-in report mode exposes their diagnostic table to downstream trait QC.
    """
    if result_failure_action not in ('error', 'report'):
        raise ValueError('result_failure_action must be error or report')
    if not log_files:
        raise RuntimeError('No correlation results found in current LDSC outputs')
    print(f"  -> Reading {len(log_files)} LDSC result files...")
    frames, statuses = [], []
    for path in log_files:
        if os.fspath(path).endswith('.log'):
            frame = _read_legacy_log(path)
            status = result_status(frame, path)
        else:
            frame = _read_numerical_csv(path, allow_failed=True)
            status = frame.attrs.pop('pair_status')
        frames.append(frame)
        statuses.append(status)
    df = pd.concat(frames, ignore_index=True)
    pairs = pd.concat(statuses, ignore_index=True)
    for column in ['p1', 'p2']:
        df[column] = df[column].apply(_trait_name)
        pairs[column] = pairs[column].apply(_trait_name)
    df = df[df['p1'] != 'p1'].drop_duplicates()
    expected = {(ref, target) for ref in trait_metadata.loc[trait_metadata['ref'] == 'yes', 'gwas_name']
                for target in trait_metadata['gwas_name']}
    actual = set(zip(df['p1'], df['p2']))
    if expected != actual:
        raise RuntimeError(f'LDSC comparison mismatch: missing={sorted(expected - actual)}; unexpected={sorted(actual - expected)}')
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
    order = {name: i for i, name in enumerate(trait_metadata.gwas_name)}
    for table in (df, pairs):
        table.sort_values(['p1', 'p2'], key=lambda x: x.map(order), kind='stable', inplace=True)
        table.reset_index(drop=True, inplace=True)
    failed = pairs.Status.eq('failed_estimate')
    diagnostic = os.path.join(output_folder, 'ldsc_results_diagnostic.csv')
    write_results_csv(df, diagnostic)
    write_results_csv(pairs, os.path.join(output_folder, 'LDSC_Pair_Status.csv'))
    write_results_csv(trait_status(df, pairs, trait_metadata.gwas_name),
                      os.path.join(output_folder, 'LDSC_Trait_Status.csv'))
    output = diagnostic if failed.any() else os.path.join(output_folder, 'ldsc_results.csv')
    if not failed.any():
        write_results_csv(df, output)
    status_path = os.path.join(output_folder, 'LDSC_Compilation_Status.csv')
    write_results_csv(pd.DataFrame([{
        'Status': 'estimation_failures' if failed.any() else 'compiled',
        'Action': result_failure_action, 'Result_File': output, 'Failed_Rows': int(failed.sum()),
    }]), status_path)
    if failed.any():
        detail = '; '.join(dict.fromkeys(reason for text in pairs.loc[failed, 'Reason']
                                       for reason in text.split('; ')))
        message = (f'{int(failed.sum())} LDSC result row(s) failed QC: {detail}. '
                   f'All estimates and source locations saved in {diagnostic} and LDSC_Pair_Status.csv. '
                   'No estimates were imputed. Use --result_failure_action report for result collection; '
                   'genomicPCA remains strict unless --failed_ldsc_action drop_traits is explicitly selected.')
        if result_failure_action == 'error':
            raise EstimationFailure(message)
        print('WARNING: ' + message)
    if pairs.Warning.ne('').any():
        print('WARNING: LDSC diagnostic warnings recorded in LDSC_Pair_Status.csv.')
    return output


def check_saved_filters(metadata, filters, trait):
    """Reject changed filters; legacy files cannot establish filter provenance."""
    saved = metadata.get('filters')
    if saved is None:
        print(f'WARNING: {trait}: previous filter settings are unknown; rerun without --ldsc_only for verified filter provenance.')
    elif saved != filters:
        raise ValueError(f'{trait}: filter settings differ from the saved run; rerun without --ldsc_only')
