#!/usr/bin/env python3
"""Backward-compatible CLI; install this checkout with python -m pip install -e ."""
from ldsc_gpca.ldsc_cli import main

if __name__ == "__main__":
    raise SystemExit(main())
