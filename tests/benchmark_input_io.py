"""Measure preparation/filtering and compare every output against a baseline.

Run this optional harness with a Python containing psutil; --python selects the
project's Python. No package imports occur in the monitoring process. Inputs are
never modified. Supply a fresh --output and keep --metrics outside that folder.
The config JSON contains manifest, hm3, munged_dir, traits (list), ld_ref and
optionally ld_weights. Only keys relevant to --workload are required.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def fingerprints(folder):
    result = {}
    for path in sorted(folder.rglob('*')):
        if not path.is_file():
            continue
        if path.suffix == '.gz':
            with gzip.open(path, 'rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        elif path.suffix in ('.csv', '.json'):
            data = path.read_text().replace(str(folder), '<OUTPUT>')
            if path.name.endswith('Worker_Attempts.csv'):
                lines = data.splitlines()
                data = '\n'.join(lines[:1]+sorted(lines[1:]))+'\n'
            digest = hashlib.sha256(data.encode()).hexdigest()
        else:
            with path.open('rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        result[str(path.relative_to(folder))] = digest
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--metrics', type=Path, required=True)
    parser.add_argument('--compare', type=Path)
    parser.add_argument('--workload', choices=['prepare-gpca','prepare-legacy','prepare-ldsc','prepare-both','auto','fixed'], required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--polars_threads', type=int)
    parser.add_argument('--clear_thread_env', action='store_true', help='Benchmark controlled defaults; otherwise preserve the caller environment.')
    args = parser.parse_args()
    import psutil  # Optional benchmark dependency, not a package dependency.
    psutil.Process().children(recursive=True)  # Check monitor permissions before launching work.
    config = json.loads(args.config.read_text())
    out = args.output.resolve()
    if out.exists():
        parser.error('--output must be a fresh path')
    if args.metrics.exists() or args.metrics.resolve().is_relative_to(out):
        parser.error('--metrics must be a fresh path outside --output')
    out.mkdir(parents=True)
    env = dict(os.environ, PYTHONPATH=str(args.source.resolve()/'src'), PYTHONDONTWRITEBYTECODE='1')
    if args.clear_thread_env:
        for key in ('POLARS_MAX_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
            env.pop(key, None)
    if args.polars_threads is not None:
        env['POLARS_MAX_THREADS'] = str(args.polars_threads)
    if args.workload.startswith('prepare-'):
        mode = args.workload.removeprefix('prepare-')
        command = [args.python, '-m', 'ldsc_gpca', 'prepare', '--input', config['manifest'],
                   '--outdir', str(out), '--splitby_chr', 'nosplit', '--n_cores', str(args.workers)]
        if mode == 'legacy':
            command += ['--write_munge_inputs', '--hm3', config['hm3']]
        elif mode != 'gpca':
            command += ['--mode', mode, '--raw_only', '--hm3', config['hm3']]
    else:
        cutoff = 'auto' if args.workload == 'auto' else 4
        code = ('from ldsc_gpca.threads import configure_numerical_threads; configure_numerical_threads(); '
                'import pandas as pd; from ldsc_gpca.munging import filter_munged_sumstats; '
                f'filter_munged_sumstats({args.workers!r}, pd.DataFrame({{"gwas_name":{config["traits"]!r}}}), '
                f'{config["munged_dir"]!r}, {str(out)!r}, {cutoff!r}, '
                f'ld_ref_dir={config["ld_ref"]!r}, ld_weights_dir={config.get("ld_weights", config["ld_ref"])!r})')
        command = [args.python, '-c', code]
    peak = samples = 0
    start = time.perf_counter()
    with args.metrics.with_suffix('.log').open('w') as log:
        child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        process = psutil.Process(child.pid)
        while child.poll() is None:
            rss = 0
            try:
                for member in [process, *process.children(recursive=True)]:
                    try:
                        rss += member.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak = max(peak, rss)
                samples += 1
            except psutil.NoSuchProcess:
                pass
            time.sleep(.1)
    record = dict(command=command, source=str(args.source.resolve()), workload=args.workload, workers=args.workers,
                  elapsed_seconds=time.perf_counter()-start, sampled_peak_tree_rss_bytes=peak, samples=samples,
                  returncode=child.returncode, fingerprints=fingerprints(out))
    args.metrics.write_text(json.dumps(record, indent=2)+'\n')
    if args.compare:
        baseline = json.loads(args.compare.read_text())
        assert record['returncode'] == baseline['returncode'], 'Exit status differs'
        assert record['fingerprints'] == baseline['fingerprints'], 'Output content differs'
        print('Exact output/status parity with', args.compare)
    print(json.dumps({k:v for k,v in record.items() if k not in ('fingerprints','command')}))
    return child.returncode


if __name__ == '__main__':
    raise SystemExit(main())
