"""Combine only the GWAMA outputs recorded for a successful current run."""
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

import pandas as pd

RESULT_SUFFIX = '.N_weighted_GWAMA.results.txt.gz'
LOG_SUFFIX = '.N_weighted_GWAMA.log'
GPCA_COLUMNS = ['SNPID', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'N_eff',
                'BETA', 'SE', 'Z', 'PVAL', 'INFO']


def validate_overrides(n_eff, info_value):
    if n_eff is not None and (not math.isfinite(n_eff) or n_eff <= 0):
        raise ValueError('--gwama-output-n-eff must be finite and > 0')
    if info_value is not None and (not math.isfinite(info_value) or not 0 <= info_value <= 1):
        raise ValueError('--gwama-output-info must be finite and in [0,1]')


def filename_component(value):
    if (not value or value in ('.', '..') or '/' in value or '\\' in value
            or any(ord(c) < 32 for c in value)):
        raise ValueError(f'Invalid output filename component: {value!r}')
    return value


def snapshot_outputs(outdir):
    """Capture file identities before R starts, to reject unchanged old results."""
    folder = Path(outdir)
    return {p.name: (p.stat().st_ino, p.stat().st_size, p.stat().st_mtime_ns)
            for p in folder.iterdir() if p.is_file()
            and (p.name == 'GWAMA_Run_Status.csv' or p.name.endswith(RESULT_SUFFIX))} if folder.is_dir() else {}


def current_run_files(outdir, previous_files):
    status_file = outdir / 'GWAMA_Run_Status.csv'
    current = snapshot_outputs(outdir)
    if not status_file.is_file():
        raise ValueError(f'Missing run status: {status_file}')
    if previous_files is not None and current.get(status_file.name) == previous_files.get(status_file.name):
        raise ValueError('GWAMA run status was not updated by the current run')
    status = pd.read_csv(status_file, dtype=str, keep_default_na=False)
    if not {'Success', 'Output'}.issubset(status.columns) or status.empty:
        raise ValueError('GWAMA run status is empty or missing Success/Output columns')
    if not status['Success'].str.lower().eq('true').all():
        raise ValueError('GWAMA was skipped or at least one run failed; post-processing stopped')
    names = [filename_component(v) for v in status['Output']]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate output prefixes in GWAMA run status')
    files = [outdir / (name + RESULT_SUFFIX) for name in names]
    for path in files:
        if not path.is_file():
            raise ValueError(f'Missing current-run GWAMA result: {path}')
        if previous_files is not None and current.get(path.name) == previous_files.get(path.name):
            raise ValueError(f'GWAMA result was not updated by the current run: {path}')
    logs = [outdir / (name + LOG_SUFFIX) for name in names]
    return files, [p for p in logs if p.is_file()]


def combine_results(files, n_eff, info_value):
    required = set(GPCA_COLUMNS) | {'Direction'}
    if n_eff is not None:
        required.remove('N_eff')
    if info_value is not None:
        required.remove('INFO')
    frames = []
    for path in files:
        frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        missing = required - set(frame.columns)
        if frame.empty or missing:
            raise ValueError(f'{path.name}: empty result or missing columns {sorted(missing)}')
        if not frame['Direction'].str.fullmatch(r'[+?\-]+').all():
            raise ValueError(f'{path.name}: Direction must contain only +, - and ? and be non-empty')
        for label, pattern in [('question', r'\?'), ('plus', r'\+'), ('minus', r'\-')]:
            frame[f'count_{label}'] = frame['Direction'].str.count(pattern)
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    if combined['SNPID'].str.strip().eq('').any() or combined['SNPID'].duplicated().any():
        raise ValueError('SNPID values must be non-empty and unique across current-run results')
    chromosomes = combined['CHR'].str.replace(r'^chr', '', regex=True).str.upper()
    chromosomes = chromosomes.replace({'X': '23', 'Y': '24', 'XY': '25', 'M': '26', 'MT': '26'})
    keys = pd.DataFrame({'chr': pd.to_numeric(chromosomes, errors='coerce'),
                         'bp': pd.to_numeric(combined['BP'], errors='coerce')})
    for column in keys:
        values = keys[column]
        if not (values.map(math.isfinite) & values.gt(0) & values.mod(1).eq(0)).all():
            raise ValueError('CHR/BP must provide valid positive integer sort positions (X/Y/XY/MT also accepted)')
    combined = combined.loc[keys.sort_values(['chr', 'bp'], kind='stable').index].reset_index(drop=True)
    summary = combined.copy()
    if n_eff is not None:
        summary['N_eff'] = n_eff
    if info_value is not None:
        summary['INFO'] = info_value
    return combined, summary[GPCA_COLUMNS].copy()


def save_outputs(combined, summary, combined_path, summary_path, audit_path, audit):
    """Stage all outputs, then publish without overwriting existing files."""
    staged, published = [], []
    try:
        for destination, data in [(combined_path, combined), (summary_path, summary), (audit_path, audit)]:
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=destination.parent, prefix='.gpca-', delete=False) as handle:
                temporary = Path(handle.name)
            staged.append((temporary, destination))
            if isinstance(data, pd.DataFrame):
                data.to_csv(temporary, sep='\t', index=False, compression='gzip')
            else:
                temporary.write_text(json.dumps(data, indent=2) + '\n')
        for temporary, destination in staged:
            os.link(temporary, destination)  # Atomic creation; fails if destination exists.
            published.append(destination)
    except Exception:
        for destination in published:
            destination.unlink()
        raise
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def process_gwama_results(outdir, harmonised_output, name=None, n_eff=None,
                          info_value=None, archive=False, previous_files=None):
    validate_overrides(n_eff, info_value)
    outdir, harmonised_output = Path(outdir).resolve(), Path(harmonised_output).resolve()
    name = filename_component(name if name is not None else outdir.name)
    combined_path = outdir / f'{name}_GWAMA_combined_results.txt.gz'
    summary_path = harmonised_output / f'{name}_GPCA_inputs.txt.gz'
    audit_path = outdir / f'{name}_postprocess.json'
    for path in [combined_path, summary_path, audit_path]:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}; choose a new --postprocess-name or output folder')
    files, logs = current_run_files(outdir, previous_files)
    archive_dir = outdir / 'chromosome_wise'
    if archive:
        for path in files + logs:
            if (archive_dir / path.name).exists():
                raise FileExistsError(f'Archive destination already exists: {archive_dir / path.name}')
    combined, summary = combine_results(files, n_eff, info_value)
    audit = {'sources': [str(p) for p in files], 'rows': len(combined),
             'combined_output': str(combined_path), 'summary_output': str(summary_path),
             'overrides': {'N_eff': n_eff, 'INFO': info_value},
             'override_scope': 'selected-column summary only; combined output preserves reported values',
             'archive_requested': archive,
             'archive_destination': str(archive_dir) if archive else None}
    save_outputs(combined, summary, combined_path, summary_path, audit_path, audit)
    if archive:
        archive_dir.mkdir(exist_ok=True)
        for path in files + logs:
            if (archive_dir / path.name).exists():
                raise FileExistsError(f'Outputs saved, but archive destination now exists: {archive_dir / path.name}')
            shutil.move(str(path), str(archive_dir / path.name))
    print(f'Combined GWAMA results: {combined_path}\nSelected-column summary: {summary_path}\nPost-processing audit: {audit_path}')
    return combined_path, summary_path
