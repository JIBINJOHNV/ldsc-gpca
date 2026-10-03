#!/usr/bin/env python3
"""Refresh R help shipped with the Python CLI, or check that it is current."""
import argparse
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable name or path.')
    parser.add_argument('--check', action='store_true', help='Fail if the checked-in help needs updating.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1] / 'src/ldsc_gpca/r'
    stale = []
    for backend, script in [('python_ldsc', 'gpsca_gwama_python_ldsc.r'), ('genomicsem', 'gpsca_gwama_v2.r')]:
        result = subprocess.run([args.rscript, str(root/script), '--help'],
            env={**os.environ, 'NO_COLOR': '1'}, capture_output=True, text=True, check=True)
        expected = result.stdout.rstrip() + '\n'
        path = root/backend/'help.txt'
        if args.check:
            if path.read_text() != expected: stale.append(str(path))
        else:
            path.write_text(expected)
    if stale: parser.exit(1, 'Stale CLI help; run scripts/update_cli_help.py:\n' + '\n'.join(stale) + '\n')
    return 0


if __name__ == '__main__': raise SystemExit(main())
