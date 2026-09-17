# Managed numerical LDSC export validation

Validated 2026-09-16 for package version 0.5.3.

## Output contract

The existing isolated LDSC process calls native regression code once per batch.
The exporter observes the original DataFrame in `_get_rg_table`, saves it with
`%.17g` CSV formatting, and leaves readable rendering unchanged. The compiler
reads each numerical CSV once using pandas' round-trip float parser; no log is
read during managed compilation. Genuine invalid estimates are rejected, not
imputed. The raw `ldsc.py` passthrough and R statistical implementations are
unchanged. Existing `.log` API inputs retain explicit legacy recovery.

The upstream behavior was verified at the pinned revision:

- [Native result-table construction](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/sumstats.py)
- [Human display formatting](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldsc.py)

## Representative tests

- Exact round-trip of every correlation, heritability and intercept field,
  including `se = 1.9647123456789123e-06` and a scientific-notation p-value.
- Native table builder/rendering called once per batch; hooks restored after
  success or failure, including simulated disk-full errors.
- Original CSV survives failed replacement; temporary files are cleaned up.
- Unsupported schemas, missing/stale exports, invalid identifiers, missing
  pairs, non-finite values, non-positive SEs and out-of-range p-values fail.
- P-value boundaries 0/1 and the smallest positive binary64 SE are preserved.
- Three-trait shuffled complete pairs and observed/liability scale mapping.
- Exact duplicate rows collapse; conflicting estimates remain for the existing
  downstream duplicate-conflict validation, rather than being discarded.
- Legacy precision-recovery tests continue to pass.

The optional real-runtime tests can be reproduced with:

```bash
export LDSC_GPCA_TEST_LDSC_PREFIX=/absolute/path/to/ldsc/environment
export LDSC_GPCA_TEST_CONDA=/absolute/path/to/conda
PYTHONPATH=src python -B -m unittest discover -s tests -p 'test_ldsc_export_integration.py' -v
```

The real-runtime fixture uses 22,000 simulated variants across 22 chromosomes,
three quantitative traits, N=50,000, seed 294, and synthetic LD scores. These
are software test inputs, not biological findings. Expected: three regression
batches, nine pair rows, finite positive required SEs, near-unit self-pair rg,
and self-pair SEs below 0.00005 retained instead of rounded to zero.

Observed: the installed LDSC 3.0.2 runtime passed both integration tests. The
exact pinned LDSC 3.0.1 source at commit
`6c673952cee74bd5c57aef1555a03b1c015399a0` was independently run on the same
fixture using the child environment's Python 3.10.14, pandas 1.5.0, NumPy
1.23.3 and SciPy 1.9.2. All three batches/nine pair rows passed, with positive
tiny self-pair SEs and no legacy-log reads.

Parent-package tests used Python 3.10.21, pandas 2.0.3, NumPy 1.23.5 and
Polars 0.20.31. Wheel construction succeeded and includes the standalone exporter.

## Known unrelated validation limitations

The broader suite has three pre-existing failures, reproduced in the original
source: two R-help tests fail because the local R installation lacks `glue`,
and an installer-routing fixture fails its environment-creator preflight.
These do not exercise result export. Full R/PCA/GWAMA revalidation was not
performed in this turn; no R calculations or GWAMA implementation were changed.
The user's historical 43-trait analysis was not rerun by this implementation.
