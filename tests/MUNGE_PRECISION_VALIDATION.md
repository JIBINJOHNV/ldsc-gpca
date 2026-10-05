# Managed munging precision validation

Validated on 2026-10-03 against parent commit `2e83ff8`.

Managed munging now calls the pinned native `munge_sumstats(args, p=False)` and
writes its returned table once using `%.17g`. Native filtering, allele matching,
P-to-Z calculation, signs, sample sizes, row/column order and missing values are
retained. Only output formatting changes. The installed upstream code and raw
`ldsc-gpca munge_sumstats.py` passthrough remain unchanged.

The writer publishes a completed temporary gzip by atomic replacement. Failed
writes preserve an existing output and clean up the temporary file. Logs, runtime
metadata and prevalence sidecars identify the format. Restart fingerprints now
include the native munging source and wrapper; old checkpoints are invalidated.
`--ldsc_only` explicitly reports reused files with unspecified precision, because
it cannot recover digits already lost.

## Native and boundary checks

The unchanged [CBIIT source at 6c67395](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py)
was downloaded and installed with `--no-deps` in the local test virtual
environment. Native munging was executed directly and through the standalone
wrapper. The returned native values matched the new file exactly with a
round-trip reader. Exactly one table write occurred. The Python LDSC parser also
accepted the file, with ordinary floating-point reading tolerance of `1e-15`.

Fixtures cover positive/negative/zero and very small Z, fractional N, optional
frequency output, absent alleles, missing data and reference-only rows,
palindromic SNP and MAF exclusions, column order, invalid schema, missing native
runtime, native help/invalid arguments, write failure, cleanup, command routing,
format metadata and restart invalidation.

For computed Z approximately 8.9444, the native three-decimal output is 8.944:
its squared value is 79.995136 and passes cutoff 80. The new file retains the
computed value, whose square is approximately 80.00229136, and correctly excludes
it. The negative counterpart is also excluded; below-threshold and zero values
remain, and missing rows are preserved for LDSC to discard.

## Representative trait-sized data

One full saved trait file from each original dataset was used. These were already
filtered LDSC inputs, not original VCFs. Test-only munging fixtures reconstructed
P from the saved Z and used their SNP/allele lists as the merge reference;
`--n-min 0` prevented additional sample-size exclusion. This tests native
munging/output compatibility on realistic values and file sizes. It does not
recover original raw-data precision or constitute a full LDSC/PCA/GWAMA rerun.

| Dataset / trait | Rows | New output versus native computed values | Cutoff-80 result |
| --- | ---: | --- | --- |
| C2 / CCL25_O15444_OID20674_v1_Inflammation | 1,189,615 | Exact numerical and row-order match | All retained, as expected for these prefiltered inputs |
| NOISE / IDO1_P14902_OID30563_v1_Inflammation_II | 1,189,765 | Exact numerical and row-order match | All retained, as expected for these prefiltered inputs |

Original source SHA-256 hashes were unchanged. Native three-decimal rounding
changed Z by up to approximately 0.0005; the precision-preserving writer matched
the native in-memory values exactly. Python LDSC parsing of both output files
was checked separately: maximum absolute Z-reading difference was 8.88e-16
(the native pandas reader does not request round-trip conversion). Missing rows and threshold-crossing behavior are covered
by the controlled fixtures above.

Single runs took 12.11/12.54 seconds with native writing and 13.52/13.27 seconds
through the new wrapper for C2/NOISE. The latter timings include subprocess
startup; these are observations, not controlled performance claims. Gzip files
grew from about 7.06 MB to 15.95 MB. There is no second munged-file reading pass.

## Automated tests and packaging

The targeted run passed 48 tests. The full suite ran 175 tests: 172 passed and
three existing optional Conda-based regression integrations were skipped because
no usable child Conda environment was configured. The new native munging tests
ran successfully in the test virtual environment. The R/PCA acceptance checks
were included. After adding the native-parser assertion, all eight munging
precision tests passed again. A built wheel contains the exact new wrapper.

Validation used Python 3.12.14, NumPy 1.26.4, pandas 2.3.3, SciPy 1.13.1,
bitarray 2.9.3 and R 4.4.2. Reproduce with a matching native runtime:

```bash
PYTHONPATH=src:tests \
LDSC_GPCA_TEST_MUNGE_SCRIPT=/path/to/native/munge_sumstats.py \
python -m unittest test_munge_precision test_munging_filter test_restart -v

PYTHONPATH=src:tests \
LDSC_GPCA_TEST_MUNGE_SCRIPT=/path/to/native/munge_sumstats.py \
LDSC_GPCA_TEST_RSCRIPT=/path/to/Rscript \
python -m unittest discover -s tests -v
```

Private fixtures, `validate_real.py`, hashes, logs and JSON results are retained
in the sibling `munge_precision_validation_20261003/` directory, outside Git.
