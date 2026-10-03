# I/O acceleration validation

Validated on 2026-10-03 against baseline commit
`676ae359b12b0c6d6352a4bdf3b4534ce1dd219b`.

## Changes and compatibility

- Share the existing pigz/gzip output implementation between chi-square filtering
  and final GWAMA export. Filtering keeps gzip level 9; final export retains its
  existing configurable level. No dependency versions or defaults were upgraded.
- Detect pigz once per stage, fall back to Python gzip when unavailable, exclude
  GZIP/PIGZ environment overrides, and verify compressor exit status before
  publishing staged files. Broken pipes report the compressor error. Existing
  two-attempt trait-worker retries and temporary-file cleanup remain active.
- Divide compression workers between two output streams and active trait jobs,
  with a minimum of one worker per stream. Thus --n_cores is not a hard upper
  bound on all helper/I/O threads when many traits run concurrently. The combined
  exclusion writer runs after trait jobs and can use the full requested count.
- Read standard TSV LD-score SNP columns using Polars. Retain the pandas
  whitespace reader for padded, mixed-whitespace, quoted or empty-tab-field
  inputs; preserve pandas 2.x missing-value tokens and string identifiers.
  One chromosome is buffered at a time; the shared SNP set is still built once.
- Parse LDSC result CSVs as Polars strings, then apply the existing Python float
  conversion and numerical QC. Preserve leading-zero and NA-like trait names,
  precise decimals, subnormals and the zero-SE self-pair exception. Explicitly
  reject duplicate headers and malformed quoted input; Polars 0.20.31 alone can
  accept an unterminated final quoted field. No PyArrow conversion is needed.

Validated installed versions: Python 3.12.14, Polars **0.20.31**, pandas 2.3.3,
NumPy 1.26.4, pigz **2.8**, R 4.4.2, data.table 1.17.8, macOS ARM64. These Python
libraries satisfy pyproject.toml; the tested Polars and pigz versions also match
those pinned in environment.yml. The installed pandas version differs from that
file's 2.2.3 pin, but falls within the supported >=2,<3 range. This was not a new
installation or a test of every supported version/platform.

All 13 protected vendor/backend/adapter files matched baseline hashes, including
the supplied GWAMA implementation, GenomicSEM modules and LDSC adapters. No R
code, native regression, statistical formula, trait ordering or QC policy changed.

## Original-data equivalence

Saved complete LDSC tables:

| Dataset | Traits | Directed pair rows | Available real Z/N GWAMA subset |
| --- | ---: | ---: | ---: |
| T3_C0_C1_C2 | 43 | 1,849 | 521 SNPs |
| T3_C0_C1_NOISE | 49 | 2,401 | 259 SNPs |

Baseline and changed compilation agree exactly, including pair QC. All four
combinations of correlation/covariance PCA and pair/trait-wide normalization were
rerun per dataset. All 168 diagnostic CSVs matched exactly after round-trip
parsing, including matrices, eigenvalues and loadings. Every GWAMA output column
also matched exactly. Path-bearing run-status files were checked for successful
completion rather than compared as numerical tables.

The shared compression helper was additionally checked through final export in
all eight runs: all 16 decompressed exported files were byte-identical. The saved
GWAMA outputs lack INFO, so both export versions used an explicit INFO=1 test
override; that value is a fixture constant, not an estimated imputation quality.

Filtering used one full original saved munged file per dataset:
AGR2_O95994_OID20896_v1_Neurology (1,189,835 rows) and
ADAM8_P78325_OID21039_v1_Neurology (1,189,720 rows). Both pigz and the gzip fallback
matched baseline retained text, per-trait exclusions, combined exclusions and QC
exactly after decompression (apart from deliberately different output paths).
Checks covered cutoff 80, auto, and cutoff 4. The lower cutoff actually excluded
60,379 and 68,677 SNPs respectively.
Both auto cutoffs were 80. All 22 original LD-score files produced the identical
1,290,028-SNP reference set. Original input/reference hashes were unchanged.

These munged files were already filtered at 80; no removed SNPs were recovered.
The GPCA/GWAMA checks use saved complete LDSC tables and reconstructed real Z/N
subsets with shared reference EAF. They are not new whole-genome LDSC regressions,
full original GWAS GWAMA runs or an independent BETA/SE calibration assessment.

## Performance

Full fixed-80 filtering, including reading, filtering, compression and combined
exclusion output, with --n_cores 4 and one trait per run (two compressor workers
per output stream):

| Dataset | Baseline gzip | New pigz | Ratio |
| --- | ---: | ---: | ---: |
| T3_C0_C1_C2 | 4.136 s | 1.756 s | 2.36x |
| T3_C0_C1_NOISE | 4.039 s | 1.714 s | 2.36x |

These are medians of three alternating measured runs per version after warm-up.
They differ from the earlier four-workers-per-compressor prototype because the
implementation shares the budget between both output streams. OS caches were
not cleared. They measure this stage, not the whole pipeline or throughput over
many simultaneous traits.

Loading all 22 LD-score SNP files fell from 0.833 to
0.683 seconds, with exact SNP-set agreement (median of
three runs after warm-up). Result import plus QC fell from
0.033 to 0.023 seconds for C2 and
0.041 to 0.028 seconds for NOISE
(median of seven runs after warm-up). Small-result savings are milliseconds.

## Automated checks and reproduction

The full suite ran 194 tests: 189 passed and five optional tests were skipped.
Two skipped native munging tests were then run with the available pinned LDSC
source and both passed, giving **191 verified passes and three unavailable
optional LDSC integration tests**. The latter require an isolated child runtime.
The focused filter/result/export/restart suite passed all 67 tests.

Coverage includes normal rows, exact thresholds, missing Z and placeholders,
invalid numbers, malformed CSV and whitespace compatibility, literal trait/SNP
identifiers, exact numeric precision, source preservation, manifest order,
compressor failures and broken pipes, two attempts then stop, preservation of
existing files on failure, gzip fallback, ignored GZIP/PIGZ overrides, worker
allocation, restart content hashing, and final export behavior.

```bash
PYTHONDONTWRITEBYTECODE=1 POLARS_MAX_THREADS=4 \
  LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript PYTHONPATH=src \
  python -m unittest discover -s tests -v
LDSC_GPCA_TEST_MUNGE_SCRIPT=/path/to/pinned/munge_sumstats.py \
  PYTHONPATH=src:tests python -m unittest test_munge_precision.NativeMungingTests -v
```

Local original-data scripts, baseline snapshot, outputs, timing JSON, input and
protected-source hashes, and logs are in the sibling directory
`../io_speed_validation_20261003/`. Entry scripts are `validate_real_data.py`
(before, after_final, compare), `validate_filter.py`, `validate_exports.py` and
`quick_bench.py`. Private data and generated results are not committed.

API references: [Polars 0.20 read_csv](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.read_csv.html),
[pigz manual](https://www.zlib.net/pigz/pigz.pdf),
[pandas read_csv](https://pandas.pydata.org/pandas-docs/version/2.2/reference/api/pandas.read_csv.html).
