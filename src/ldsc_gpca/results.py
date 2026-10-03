"""LDSC result validation, compilation and reuse checks."""
import os
import math
import re
import csv
import uuid
import pandas as pd
from .ldsc_export import BASE_RESULT_COLUMNS, HERITABILITY_COLUMNS, NATIVE_COLUMNS, write_results_csv
from .result_qc import result_status, trait_status, _missing_value, self_rg_zero_se
from .trait_selection import select_complete_traits


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
    estimates = pd.DataFrame(rows)
    for column in ('rg', 'se', 'z', 'p'):
        estimates[column] = pd.to_numeric(estimates[column], errors='coerce')
    for row, allowed_zero in zip(rows, self_rg_zero_se(estimates)):
        se = pd.to_numeric(row['se'], errors='coerce')
        if not math.isfinite(se) or (se <= 0 and not allowed_zero):
            pair = (_trait_name(row['p1']), _trait_name(row['p2']))
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
    native_columns = [c for c in NATIVE_COLUMNS if c in frame]
    if native_columns and len(native_columns) != len(NATIVE_COLUMNS):
        raise RuntimeError(f'Incomplete native LDSC columns in {path}: require {NATIVE_COLUMNS}')
    # Preserve native exports; unrelated annotation columns remain excluded.
    frame = frame[list(BASE_RESULT_COLUMNS) + [c for pair in h2_columns for c in pair] + native_columns].copy()
    for column in native_columns:
        try:
            frame[column] = frame[column].where(~_missing_value(frame[column]), float('nan')).map(float)
        except (ValueError, TypeError) as error:
            raise RuntimeError(f'Non-numeric LDSC {column} in {path}') from error
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
                     'LDSC_Pair_Status.csv', 'LDSC_Trait_Status.csv',
                     'LDSC_Retained_Traits.csv', 'LDSC_Dropped_Traits.csv'):
        path = os.path.join(output_folder, filename)
        if os.path.exists(path):
            os.replace(path, path + '.previous-' + uuid.uuid4().hex)
    status_path = os.path.join(output_folder, 'LDSC_Compilation_Status.csv')
    write_results_csv(pd.DataFrame([{'Status': 'incomplete', 'Action': result_failure_action,
                                    'Result_File': '', 'Failed_Rows': ''}]), status_path)


def compile_results(output_folder, log_files, trait_metadata, *, result_failure_action='error',
                    input_manifest=None, execution_failures=None):
    if result_failure_action not in ('error', 'report', 'drop_traits'):
        raise ValueError('result_failure_action must be error, report or drop_traits')
    prepare_compilation_outputs(output_folder, result_failure_action)
    try:
        return _compile_results(output_folder, log_files, trait_metadata,
                                result_failure_action=result_failure_action,
                                input_manifest=input_manifest, execution_failures=execution_failures)
    except EstimationFailure:
        raise
    except (RuntimeError, ValueError, OSError) as error:
        write_results_csv(pd.DataFrame([{'Status': 'structural_failure',
            'Action': result_failure_action, 'Result_File': '', 'Failed_Rows': '',
            'Error': str(error)}]), os.path.join(output_folder, 'LDSC_Compilation_Status.csv'))
        raise


def _publish_retained_results(output_folder, frame, pairs, metadata, input_manifest):
    retained, excluded = select_complete_traits(frame, pairs, metadata.gwas_name)
    if input_manifest is None:
        input_manifest = metadata.rename(columns={'gwas_name': 'traitname',
                                                  'pop_prevalence': 'population_prevalence'}).copy()
        if 'sample_prevalence' not in input_manifest:
            input_manifest['sample_prevalence'] = float('nan')
    selected_manifest = input_manifest[input_manifest.traitname.isin(retained)].copy()
    if selected_manifest.traitname.tolist() != retained:
        raise RuntimeError('Retained manifest does not match LDSC trait names/order')
    manifest_path = os.path.join(output_folder, 'LDSC_Retained_Traits.csv')
    write_results_csv(selected_manifest, manifest_path)
    write_results_csv(excluded, os.path.join(output_folder, 'LDSC_Dropped_Traits.csv'))
    traits = trait_status(frame, pairs, metadata.gwas_name)
    traits['Retained'] = traits.Trait.isin(retained)
    traits['Exclusion_Reason'] = traits.Trait.map(excluded.set_index('Trait').Reason).fillna('')
    write_results_csv(traits, os.path.join(output_folder, 'LDSC_Trait_Status.csv'))
    status = {'Status': 'compiled_with_trait_exclusions' if len(excluded) else 'compiled',
              'Action': 'drop_traits', 'Result_File': '',
              'Failed_Rows': int(pairs.Status.ne('valid').sum()),
              'Retained_Traits': len(retained), 'Excluded_Traits': len(excluded),
              'Retained_Manifest': manifest_path}
    status_path = os.path.join(output_folder, 'LDSC_Compilation_Status.csv')
    if len(retained) < 2:
        status['Status'] = 'insufficient_traits'
        status['Error'] = f'Trait removal leaves {len(retained)} trait(s); at least two are required.'
        write_results_csv(pd.DataFrame([status]), status_path)
        raise EstimationFailure(status['Error'] + ' Exclusion and diagnostic reports were saved.')
    selected = frame[frame.p1.isin(retained) & frame.p2.isin(retained)].copy()
    output = os.path.join(output_folder, 'ldsc_results.csv')
    write_results_csv(selected, output)
    status['Result_File'] = output
    write_results_csv(pd.DataFrame([status]), status_path)
    print(f'  -> Trait selection: retained {len(retained)}, excluded {len(excluded)}. '
          f'Use {manifest_path} with {output} for genomicPCA.')
    if pairs.Warning.ne('').any():
        print('WARNING: LDSC diagnostic warnings recorded in LDSC_Pair_Status.csv.')
    return output


def _compile_results(output_folder, log_files, trait_metadata, *, result_failure_action,
                     input_manifest=None, execution_failures=None):
    """Compile fresh numerical CSVs; explicitly supplied .log paths retain legacy support.

    The historical log_files parameter name is retained for API compatibility.
    Structural failures remain fatal. Estimation failures are always audited;
    report preserves them and drop_traits selects a complete subset immediately.
    """
    dropping = result_failure_action == 'drop_traits'
    if dropping and not trait_metadata.ref.eq('yes').all():
        raise RuntimeError('drop_traits requires ref=yes for every trait so self-pair QC is available')
    if not log_files and not dropping:
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
    df = (pd.concat(frames, ignore_index=True) if frames else
          pd.DataFrame(columns=list(BASE_RESULT_COLUMNS) + [c for pair in HERITABILITY_COLUMNS for c in pair]))
    pairs = (pd.concat(statuses, ignore_index=True) if statuses else
             pd.DataFrame(columns=['p1', 'p2', 'Source_File', 'Source_Row', 'Status', 'Reason', 'Warning']))
    for column in ['p1', 'p2']:
        df[column] = df[column].apply(_trait_name)
        pairs[column] = pairs[column].apply(_trait_name)
    df = df[df['p1'] != 'p1'].drop_duplicates()
    expected = {(ref, target) for ref in trait_metadata.loc[trait_metadata['ref'] == 'yes', 'gwas_name']
                for target in trait_metadata['gwas_name']}
    actual = set(zip(df['p1'], df['p2']))
    if actual - expected or (expected - actual and not dropping):
        raise RuntimeError(f'LDSC comparison mismatch: missing={sorted(expected - actual)}; unexpected={sorted(actual - expected)}')
    if dropping and expected - actual:
        missing = pd.DataFrame([{'p1': a, 'p2': b, 'Source_File': '', 'Source_Row': float('nan'),
            'Status': 'missing_result', 'Warning': '',
            'Reason': 'Missing requested LDSC comparison' +
                      (': ' + execution_failures[(a, b)] if execution_failures and (a, b) in execution_failures else '')}
            for a, b in sorted(expected - actual)])
        pairs = pd.concat([pairs, missing], ignore_index=True)
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
    from .normalization import add_trait_wide_columns, DERIVED_COLUMNS
    import polars as pl
    annotations = add_trait_wide_columns(pl.DataFrame(
        {name: df[name].to_numpy() for name in ('p1', 'p2', 'rg', 'h2_obs', 'h2_liab', *NATIVE_COLUMNS) if name in df},
        schema_overrides={'p1': pl.String, 'p2': pl.String}))
    for name in DERIVED_COLUMNS:
        df[name] = annotations[name].to_numpy()
    failed = pairs.Status.ne('valid')
    diagnostic = os.path.join(output_folder, 'ldsc_results_diagnostic.csv')
    write_results_csv(df, diagnostic)
    write_results_csv(pairs, os.path.join(output_folder, 'LDSC_Pair_Status.csv'))
    write_results_csv(trait_status(df, pairs, trait_metadata.gwas_name),
                      os.path.join(output_folder, 'LDSC_Trait_Status.csv'))
    if dropping:
        return _publish_retained_results(output_folder, df, pairs, trait_metadata, input_manifest)
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
                   'use --result_failure_action drop_traits to exclude failed traits during LDSC collection.')
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
