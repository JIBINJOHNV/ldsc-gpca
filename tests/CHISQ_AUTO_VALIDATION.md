# Automatic chi-square cutoff validation

Validated on 2026-10-02 against baseline commit
`19b035e1ca37a12248ffa7ef84317bee3f2efd02`.

The managed Python LDSC command accepts `--chisq_max auto` or a positive integer.
Omitting the option keeps filtering disabled. Automatic mode calculates
`max(80, 0.001 * maximum_matched_n)` separately for each trait before chi-square
filtering and before forming trait pairs. Missing summary-statistic rows and
SNPs absent from either the reference or regression-weight LD scores do not
contribute to the maximum. The result is not rounded.

The existing streaming filter writes unchanged row text into separate copies.
Reference SNPs are loaded once and shared among workers. Auto requires one
additional streaming pass to find each trait's maximum N; fixed cutoffs do not
read the references or make that pass. No upstream LDSC or GenomicSEM code is
changed. Decimal CLI cutoff spellings are deliberately rejected under the new
integer-or-auto contract; use `80` instead of `80.0`.

## Direct comparison with GenomicSEM

The complete filtering block (`all_y <- lapply(...)`) was extracted without
modification from [pinned GenomicSEM source](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R#L172-L195)
and executed using R/readr, the same summary statistics and references. The
source SHA-256 is
`9fddaee7d46d1faa06061ab023bbab4e3aa97d4365a902a40523221e09eac8e5`.
Every cutoff and every retained LD-matched SNP agreed exactly.

| Controlled trait | Maximum complete matched N | Expected / observed cutoff | Matched SNPs retained |
| --- | ---: | ---: | ---: |
| A | 50,000 | 80 / 80 | 1 |
| B | 80,000 | 80 / 80 | 1 |
| C | 144,000 | 144 / 144 | 2 |
| D | 500,250 | 500.25 / 500.25 | 2 |

Each fixture included unmatched SNPs, SNPs present in only one reference, and
rows missing Z or an allele, all with N=9,999,999. None inflated the cutoff.
The Z=-12 SNP was retained at the exact cutoff 144, excluded at 80, and retained
at 500.25. Fixed-80 copies and exclusions exactly matched the baseline filter.

## Full original trait files

One complete original trait file from each dataset was tested, along with all
22 chromosome LD-score files. These source sumstats had already been filtered
at 80; the checks validate retained-input compatibility, not recovery of SNPs
removed previously or a new full LDSC/GWAMA regression run.

| Dataset / trait | Source rows | Maximum matched N | Expected / observed cutoff | Matched SNPs retained by both implementations |
| --- | ---: | ---: | ---: | ---: |
| C2 / CCL25_O15444_OID20674_v1_Inflammation | 1,189,615 | 34,049 | 80 / 80 | 1,176,027 |
| NOISE / IDO1_P14902_OID30563_v1_Inflammation_II | 1,189,765 | 32,981 | 80 / 80 | 1,176,188 |

Source hashes were unchanged. Fixed-80 filtered content and exclusion content
matched the baseline exactly after decompression. Automatic filtering produced
the same LD-matched SNP sets as the original GenomicSEM block. Because both
traits have maximum N below 80,000, the controlled fixtures above supply the
required coverage of thresholds greater than 80.

With two trait workers, single-run filtering timings were:

| Dataset | Baseline fixed 80 | New fixed 80 | New auto |
| --- | ---: | ---: | ---: |
| C2 | 4.17 s | 4.17 s | 6.69 s |
| NOISE | 4.18 s | 4.17 s | 6.70 s |

Reference loading took 0.84 seconds once for both traits. These are observations
on this machine, not general performance guarantees. Auto's extra time reflects
its preliminary N scan. Filtering was benchmarked without rerunning regressions.

## Automated checks

The 31 focused filter/CLI/restart tests passed. They cover independent traits,
the 80 floor, inclusive boundaries, non-integer automatic cutoffs, maximum N in
a later row that is itself excluded, separate weight references, missing and
whitespace-padded missing fields, malformed and invalid inputs, no matches,
all-excluded inputs, source preservation, parallel order, and CLI validation.

The managed CLI integration runs real filtering and result compilation with a
mocked native regression command. Identical auto restarts launch no new commands;
switching to fixed 80 invalidates checkpoints and reruns the commands. The native
cross-product `--chisq-max` option is never forwarded. Existing fixed-mode tests
and PCA/GWAMA acceptance tests are included in the full suite.

The final full suite ran 160 tests in 159.70 seconds: 157 passed and three
optional Python LDSC CLI integrations were skipped because their child runtime
was unavailable. All PCA/GWAMA acceptance tests in the suite passed.

```bash
PYTHONPATH=src:tests python -m unittest test_chisq_auto test_munging_filter test_restart -v
PYTHONPATH=src LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript python -m unittest discover -s tests -v
```

Validation used Python 3.12.14, pandas 2.3.3 and R 4.4.2. Three optional Python
LDSC CLI integrations require an explicitly configured usable child runtime;
that runtime is unavailable in this workspace. The direct comparison above
executes the upstream GenomicSEM filtering block, not the whole estimator.

Private copies, pinned source/hash, `validate.py`, `oracle.R`, retained SNP lists
and `summary.json` are saved in the local sibling directory
`chisq_auto_validation_20261002/`; private inputs and results are not committed.
