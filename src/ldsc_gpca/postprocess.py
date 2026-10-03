"""Combine only the GWAMA outputs recorded for a successful current run."""
import json
import gzip
import csv
import math
import os
from pathlib import Path
import shutil
import tempfile
from time import perf_counter

import numpy as np
import pandas as pd
from .threads import export_workers
from .workers import run_parallel_jobs
from .compression import compression_settings, compressed_writer

RESULT_SUFFIX = '.N_weighted_GWAMA.results.txt.gz'
LOG_SUFFIX = '.N_weighted_GWAMA.log'
GPCA_COLUMNS = ['SNPID', 'CHR', 'BP', 'EA', 'OA', 'EAF', 'N_eff',
                'BETA', 'SE', 'Z', 'PVAL', 'INFO']


def validate_overrides(n_eff, info_value):
    if n_eff is not None and (not math.isfinite(n_eff) or n_eff <= 0):
        raise ValueError('--gwama_output_n_eff must be finite and > 0')
    if info_value is not None and (not math.isfinite(info_value) or not 0 <= info_value <= 1):
        raise ValueError('--gwama_output_info must be finite and in [0,1]')


def progress(message):
    print(f'[GWAMA export] {message}', flush=True)


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


def table_engine(n_cores=0):
    # Direct Python callers also get a bounded default if Polars is not yet
    # initialized. CLI startup sets this before any preparation work.
    os.environ.setdefault('POLARS_MAX_THREADS', str(export_workers(n_cores)))
    import polars as pl
    return pl


def read_chromosome(path, required):
    import polars as pl
    started = perf_counter()
    # Reject malformed headers rather than accepting a parser's silent renaming.
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        columns = next(csv.reader(stream, delimiter='\t'), [])
    if not columns or any(not c.strip() for c in columns) or len(columns) != len(set(columns)):
        raise ValueError(f'{path.name}: empty or duplicate column headers')
    missing = required - set(columns)
    if missing:
        raise ValueError(f'{path.name}: missing columns {sorted(missing)}')
    try:
        # Polars 0.20 can route names containing glob characters through its
        # compressed-CSV-incompatible scan path even with glob=False. Bytes
        # keep such filenames literal; ordinary paths retain the direct reader.
        source = path.read_bytes() if any(c in str(path) for c in '*?[') else path
        frame = pl.read_csv(source, separator='\t', infer_schema_length=0,
                            missing_utf8_is_empty_string=True, n_threads=1, glob=False)
    except pl.exceptions.PolarsError as error:
        raise ValueError(f'{path.name}: could not read GWAMA table: {error}') from error
    if frame.is_empty():
        raise ValueError(f'{path.name}: empty result')
    if frame.columns != columns:
        raise ValueError(f'{path.name}: parsed headers differ from the source headers')
    if frame['Direction'].null_count() or not frame['Direction'].str.contains(r'\A[+?\-]+\z').all():
        raise ValueError(f'{path.name}: Direction must contain only +, - and ? and be non-empty')
    frame = frame.with_columns([
        pl.col('Direction').str.count_matches(character, literal=True).alias(f'count_{label}')
        for label, character in [('question', '?'), ('plus', '+'), ('minus', '-')]])
    return frame, perf_counter() - started


def position_keys(combined):
    """Use integer keys without touching source strings or losing large integers."""
    import polars as pl
    chromosomes = combined['CHR'].str.replace(r'^chr', '').str.to_uppercase().replace(
        {'X': '23', 'Y': '24', 'XY': '25', 'M': '26', 'MT': '26'})
    keys = []
    for name, source in [('chr', chromosomes), ('bp', combined['BP'])]:
        values = source.cast(pl.UInt64, strict=False)
        if values.null_count():
            # Preserve the existing pandas interpretation of decimal/scientific
            # positions (e.g. 2.0, 2e1). Normal integer GWAMA positions stay in
            # Polars throughout; never round integer keys through float64.
            numeric = pd.to_numeric(pd.Series(source.to_list()), errors='coerce')
            if not (np.isfinite(numeric) & numeric.gt(0) & numeric.mod(1).eq(0)).all():
                raise ValueError('CHR/BP must provide valid positive integer sort positions (X/Y/XY/MT also accepted)')
            values = pl.Series(name, numeric.to_numpy())
        elif (values == 0).any():
            raise ValueError('CHR/BP must provide valid positive integer sort positions (X/Y/XY/MT also accepted)')
        keys.append(values.rename(name))
    return pl.DataFrame(keys)


def combine_results(files, n_eff, info_value, *, n_cores=0, status_path=None):
    pl = table_engine(n_cores)
    files = list(files)
    if not files:
        raise ValueError('No current-run GWAMA results to combine')
    required = set(GPCA_COLUMNS) | {'Direction'}
    if n_eff is not None:
        required.remove('N_eff')
    if info_value is not None:
        required.remove('INFO')
    workers = min(export_workers(n_cores), len(files))
    progress(f'Polars: {workers} concurrent file reader(s), 1 CSV parser thread per file; table pool {pl.thread_pool_size()} thread(s)')
    results = run_parallel_jobs([(path.name, read_chromosome, (path, required)) for path in files],
                               workers, stage='GWAMA export reads', status_path=status_path)
    frames = []
    for number, (path, (frame, elapsed)) in enumerate(zip(files, results), 1):
        frames.append(frame)  # Manifest order, regardless of completion/retry order.
        progress(f'Read/count {number}/{len(files)}: {path.name}; {len(frame):,} rows ({elapsed:.1f}s)')
    progress('Checking SNPIDs and chromosome/position keys...')
    combined = pl.concat(frames, how='diagonal', rechunk=False)
    del frames, results, frame
    if (combined['SNPID'].null_count() or (combined['SNPID'].str.strip_chars() == '').any()
            or combined['SNPID'].is_duplicated().any()):
        raise ValueError('SNPID values must be non-empty and unique across current-run results')
    keys = position_keys(combined)
    chromosomes, positions = keys['chr'].to_numpy(), keys['bp'].to_numpy()
    ordered = ((chromosomes[1:] > chromosomes[:-1]) |
               ((chromosomes[1:] == chromosomes[:-1]) & (positions[1:] >= positions[:-1]))).all()
    if ordered:
        progress(f'{len(combined):,} variants already in chromosome/position order; skipping redundant sort.')
    else:
        progress(f'Sorting {len(combined):,} variants by chromosome and position...')
        order = keys.with_row_index('row').sort(['chr', 'bp'], maintain_order=True)['row']
        combined = combined[order]
    # Select shared column buffers; replace only the explicitly supplied metadata.
    overrides = {'N_eff': n_eff, 'INFO': info_value}
    summary = combined.select([
        pl.lit(str(overrides[c])).alias(c) if overrides.get(c) is not None else pl.col(c)
        for c in GPCA_COLUMNS])
    return combined, summary


def write_table(data, stream):
    """Use the native TSV writer while retaining the previous empty-cell format."""
    import polars as pl
    empty = data.select(pl.col(pl.String).eq('').any()).row(0, named=True)
    columns = [name for name, contains_empty in empty.items() if contains_empty]
    if columns:
        data = data.with_columns([
            pl.when(pl.col(c) == '').then(None).otherwise(pl.col(c)).alias(c) for c in columns])
    try:
        data.write_csv(stream, separator='\t', batch_size=16384, null_value='')
    except pl.exceptions.PolarsError as error:
        raise OSError(f'Could not write GWAMA table: {error}') from error


def write_compressed_table(data, destination, settings):
    """Stream directly to gzip; pigz is optional and never needs a plain TSV."""
    with compressed_writer(destination, settings) as stream:
        write_table(data, stream)


def save_outputs(combined, summary, combined_path, summary_path, audit_path, audit,
                 *, n_cores=0, gzip_level=1):
    """Stage all outputs, then publish without overwriting existing files."""
    settings = compression_settings(n_cores, gzip_level)
    audit['compression'] = settings
    timings = audit.setdefault('timings_seconds', {})
    progress(f"Compression: {settings['backend']}, level {settings['level']}, {settings['workers']} worker(s)")
    if settings['backend'] == 'python_gzip':
        progress('pigz not found; using single-worker gzip. Install pigz to enable parallel compression.')
    staged, published = [], []
    try:
        for destination, data in [(combined_path, combined), (summary_path, summary), (audit_path, audit)]:
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=destination.parent, prefix='.gpca-', delete=False) as handle:
                temporary = Path(handle.name)
            staged.append((temporary, destination))
            if destination != audit_path:
                started = perf_counter()
                progress(f'Writing {destination.name}...')
                write_compressed_table(data, temporary, settings)
                key = 'write_combined' if destination == combined_path else 'write_summary'
                timings[key] = round(perf_counter() - started, 3)
                progress(f"Wrote {destination.name} ({timings[key]:.1f}s)")
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
                          info_value=None, archive=False, previous_files=None,
                          n_cores=0, gzip_level=1):
    validate_overrides(n_eff, info_value)
    compression_settings(n_cores, gzip_level)
    outdir, harmonised_output = Path(outdir).resolve(), Path(harmonised_output).resolve()
    name = filename_component(name if name is not None else outdir.name)
    combined_path = outdir / f'{name}_GWAMA_combined_results.txt.gz'
    summary_path = harmonised_output / f'{name}_GPCA_inputs.txt.gz'
    audit_path = outdir / f'{name}_postprocess.json'
    for path in [combined_path, summary_path, audit_path]:
        if path.exists():
            raise FileExistsError(f'Refusing to overwrite {path}; choose a new --dataset_id or output folder')
    files, logs = current_run_files(outdir, previous_files)
    archive_dir = outdir / 'chromosome_wise'
    if archive:
        for path in files + logs:
            if (archive_dir / path.name).exists():
                raise FileExistsError(f'Archive destination already exists: {archive_dir / path.name}')
    started = perf_counter()
    combined, summary = combine_results(files, n_eff, info_value, n_cores=n_cores,
        status_path=outdir/'GWAMA_Export_Worker_Attempts.csv')
    combine_seconds = perf_counter() - started
    progress(f'Combined and prepared both tables ({combine_seconds:.1f}s)')
    audit = {'dataset_id': name, 'sources': [str(p) for p in files], 'rows': len(combined),
             'combined_output': str(combined_path), 'summary_output': str(summary_path),
             'overrides': {'N_eff': n_eff, 'INFO': info_value},
             'override_scope': 'selected-column summary only; combined output preserves reported values',
             'table_engine': 'polars', 'read_workers': min(export_workers(n_cores), len(files)),
             'csv_threads_per_file': 1, 'polars_threads': table_engine(n_cores).thread_pool_size(),
             'timings_seconds': {'combine_and_prepare': round(combine_seconds, 3)},
             'archive_requested': archive,
             'archive_destination': str(archive_dir) if archive else None}
    save_outputs(combined, summary, combined_path, summary_path, audit_path, audit,
                 n_cores=n_cores, gzip_level=gzip_level)
    if archive:
        archive_dir.mkdir(exist_ok=True)
        for path in files + logs:
            if (archive_dir / path.name).exists():
                raise FileExistsError(f'Outputs saved, but archive destination now exists: {archive_dir / path.name}')
            shutil.move(str(path), str(archive_dir / path.name))
    print(f'Combined GWAMA results: {combined_path}\nSelected-column summary: {summary_path}\nPost-processing audit: {audit_path}')
    return combined_path, summary_path
