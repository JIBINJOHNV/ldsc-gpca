# Incremental preparation and filtering validation — 2026-10-04

This change implements performance candidates 1–3 only. Candidate 4 (directional
LDSC fit reuse) and candidate 5 (GWAMA matrix caching) are excluded. All bundled R
sources, pairwise fits, worker retries, compression implementation and statistical
estimators are unchanged; 36 protected implementation files were hash-compared.

## Changes

1. `prepare --mode both` writes its independent GPCA/LDSC original-record reports
   from one streaming VCF read. Only issue indices/actions/reasons are retained,
   not original VCF text. Available reports are flushed even on preparation failure.
2. The automatic-cutoff pass reduces Python missing-token overhead and reuses
   immediately repeated, validated N text. The preliminary complete LD-matched
   scan remains necessary. No temporary full-data spool or persistent cache was added.
3. Preparation releases the transformed GPCA frame before building the LDSC frame.
   Standalone CLI preparation sets its shared Polars pool to `--n_cores` (default 4)
   before importing Polars. Explicit `POLARS_MAX_THREADS` and BLAS settings remain
   respected. Existing worker/compression budgets and retries are reused.

## Original-data baseline and equivalence

The clean baseline was commit `20fe7bb`, archived and executed **before** production
code edits. Original inputs were read without modification. Their before/after
SHA-256 values and all numerical measurements are recorded in
[performance_incremental_metrics.json](performance_incremental_metrics.json).

| Original data | Input rows | Relevant unchanged output |
| --- | ---: | --- |
| `PGC3_SCZ_wave3_GRCh37_merged.vcf.gz` | 587,318 | 587,318 GPCA; 97,373 raw LDSC rows |
| `bip2024_eur_no23andMe_GRCh37_merged.vcf.gz` | 530,467 | 530,467 GPCA; 95,022 raw LDSC rows |
| `AGR2_O95994_OID20896_v1_Neurology.sumstats.gz` | 1,189,835 | At cutoff 4: 1,129,456 retained; 60,379 excluded |
| `AGRP_O00253_OID20658_v1_Inflammation.sumstats.gz` | 1,189,853 | At cutoff 4: 1,122,837 retained; 67,016 excluded |

The VCFs are the original local chromosome-1 GWAS datasets, tested with `nosplit`.
The preparation benchmark supplied explicit N overrides of 126,849 and 676,661,
sample prevalences 51,788/126,849 and 50,614/676,661, and population prevalences
0.01 and 0.02. The N/sample fractions come from the first original records and
exercise the preparation contract; they are **not validated study-wide scientific
sample definitions**. No estimates were derived from these benchmark manifests.
GPCA retained original NEF in both versions.

Both saved pQTL files were already filtered at 80. Automatic mode resolved to 80
and retained every original row. Maximum matched N was 32,733 and 34,049, with
1,176,251 and 1,176,265 complete matched rows, respectively. An independent stdlib
streaming oracle verified those counts/cutoffs, exact retained row text and both
per-trait and combined exclusions at cutoffs 4 and auto. The oracle uses exactly
`1.l2.ldscore.gz` through `22.l2.ldscore.gz`; an unused `6_old` backup discovered
during validation is excluded, matching the pipeline. Reference and weight paths
were identical for this original-data run; distinct paths are covered by tests.

Every repeated before/after output fingerprint matched. Comparisons include
decompressed gzip content, table row order, values, headers, QC, manifests and
settings. Only output-directory paths and nondeterministic worker-completion row
order are normalized. Compression container bytes need not match across backends.

Original-data GPCA-only, LDSC-only and `both` outputs matched, including one-worker
versus four-worker runs and an explicit seven-thread Polars override. Legacy
SCZ export failed in both versions on duplicate rsIDs and preserved identical
failure audits without publishing tables. A separate 22-autosome fixture matched
all **91** products across default, split, legacy, LDSC, both and N-override modes.

## Repeated measurements

Three alternating baseline/changed pairs per workload, sequentially executed on
macOS 26.6.2 arm64 (12 logical CPUs). Four workers were configured; one/two-trait
filtering uses one/two actual trait workers. Runs used fresh Python processes and
warm filesystem caches. No concurrent benchmark or test suite was run during
these repetitions. This is a small workload comparison, not a universal estimate.

| Workload | Median seconds before → after | Median sampled peak RSS MiB before → after |
| --- | ---: | ---: |
| Both preparation modes, two VCFs | 11.998 → 12.408 | 1,363.8 → 1,012.9 |
| Auto filter, one trait, Python gzip | 6.716 → 6.241 | 407.5 → 404.6 |
| Auto filter, two traits, Python gzip | 7.841 → 6.939 | 405.9 → 406.8 |
| Fixed cutoff 4, two traits, Python gzip | 5.024 → 5.075 | 94.4 → 94.3 |
| Auto filter, two traits, pigz | 7.383 → 6.564 | 414.9 → 412.4 |

Preparation reduced sampled peak memory by **25.7%**, with a **3.4%** increase in
median elapsed time. Automatic filtering reduced median elapsed time by **7.1%**
(one trait), **11.5%** (two traits/gzip), and **11.1%** (two traits/pigz). Fixed-cutoff
work was unchanged in implementation and showed no demonstrated speed benefit.
The original VCFs have little GPCA QC output, so they do not establish a speedup
from avoiding two long QC scans. The shared-scan behavior is separately tested.

Timings include process startup, reference loading, compression and audit writing,
but exclude output hashing. Memory is the maximum sampled sum of RSS for the
analysis process and descendants, including bcftools/pigz, sampled about every
100 ms plus monitor overhead. Brief peaks may be missed and shared pages may be
counted more than once. These are not cProfile cumulative timings or full-pipeline
costs. Individual timings and memory measurements are retained in the JSON.

## Regression and reproduction

Added tests cover shared/overlapping/empty original-record audits, CRLF and missing
final newline, unchanged field text and multiple reasons, incomplete sources,
header collisions, GPCA/LDSC failure audits and retries, repeated/changing N,
malformed N after valid rows, padded/missing extra fields, changed inputs/references,
and thread defaults/explicit overrides before Polars initialization. Existing tests
cover independent filters, duplicate handling, publication rollback, compression
failure, worker budgets, serial/parallel ordering and unchanged numerical outputs.

Full-suite result: **301 tests run, 296 passed, five optional Python LDSC
integrations skipped** (201.182 seconds). This included full installed pinned
GenomicSEM fixtures via R 4.3.2 and GenomicSEM commit
`6b65ca5db39fdade08b0d811477be1cdd57b5039`; it is not a full original-data regression
run. Command:

```bash
PYTHONPATH=src:tests \
PATH=/Users/JJOHN41/miniconda3/bin:/Users/JJOHN41/miniconda3/envs/rnaseq_env/bin:$PATH \
LDSC_GPCA_TEST_RSCRIPT=/private/tmp/ldsc-gpca-genomicsem-Rscript \
../.venv-ldsc-normalization/bin/python -m unittest discover -s tests -p 'test_*.py'
```

The installed-package wrapper selects the existing R 4.3.2 library; it does not
replace statistical functions. The complete log is retained locally at
`/private/tmp/ldsc-performance-full-suite.log`. The reusable benchmark harness was
also executed on the original pQTL and VCF datasets with exact baseline parity.

The reusable [benchmark_input_io.py](benchmark_input_io.py) harness compares all
output fingerprints and statuses, preserving failures. It requires `psutil` in
the monitoring Python; `--python` selects the project runtime. Example:

```bash
python tests/benchmark_input_io.py \
  --config /path/to/benchmark_config.json \
  --python /path/to/project/python --source /path/to/baseline \
  --output /path/to/fresh/before --metrics /path/to/before.json \
  --workload auto --workers 4 --clear_thread_env

python tests/benchmark_input_io.py \
  --config /path/to/benchmark_config.json \
  --python /path/to/project/python --source /path/to/changed \
  --output /path/to/fresh/after --metrics /path/to/after.json \
  --compare /path/to/before.json --workload auto --workers 4 --clear_thread_env
```

The config supplies `manifest`, `hm3`, `munged_dir`, `traits` (ordered list),
`ld_ref` and optional `ld_weights`, using absolute paths. Only keys required by
the workload are needed. Select `prepare-both`, `prepare-gpca`, `prepare-ldsc`,
`prepare-legacy`, `auto` or `fixed` (cutoff 4). Set PATH identically for compared
runs to select pigz or Python gzip. On macOS, process-tree monitoring requires
permission to enumerate processes. Individual local runs, the oracle, manifests,
baseline archive and scripts are retained in `../performance_validation_20261004/`.

This validates affected preparation/filtering products on original data and
regression fixtures. It does not rerun full original-data LDSC/GWAMA regression,
prove recovery of SNPs previously excluded at 80, or establish benefits on
many-trait whole-genome production workloads. Statistical sources and GenomicSEM
two-pass estimation remain unchanged.
