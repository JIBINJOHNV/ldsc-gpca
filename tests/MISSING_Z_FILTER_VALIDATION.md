# Missing-Z filter validation

Validated on 2026-10-02 against parent commit
`df5c8e307bcb50ec47d7d910a4daa8b4d6c9c8d8` before committing.

The existing streaming chi-square filter now preserves blank, `NA`, `nan`,
`NaN`, and `.` Z placeholders without imputation. `variants_missing_z` counts
them separately. `variants_after` counts all output rows, so subtracting
`variants_missing_z` gives the number of finite-Z survivors. Invalid text,
infinity, Z-squared overflow, malformed rows, and no finite-Z survivors still
fail clearly. Threshold behavior and original upstream LDSC code are unchanged.

## Checks on original data

One complete original trait file from each dataset was compared against the
previous filter with threshold 80. These files were already filtered at 80.
Separate test copies received five additional missing-Z placeholder rows.

| Dataset / original trait | Original rows | Previous vs fixed output | Native parser after adding missing rows |
| --- | ---: | --- | --- |
| C2 / CCL25_O15444_OID20674_v1_Inflammation | 1,189,615 | Exact decompressed content match | Same 1,189,615 valid rows, exact values |
| NOISE / IDO1_P14902_OID30563_v1_Inflammation_II | 1,189,765 | Exact decompressed content match | Same 1,189,765 valid rows, exact values |

The old filter rejected the added blank row; the fixed filter preserved all five
missing rows and counted them correctly. The unchanged [pinned native LDSC parser](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/parse.py)
then discarded the placeholders, producing exactly the original finite input
tables. Source hashes were unchanged. This validates filtering/parser behavior
on representative full trait files, not a new full-dataset LDSC/GWAMA rerun.

Single-run filter timings were 4.20 versus 4.22 seconds for C2 and 4.22 versus
4.25 seconds for NOISE (fixed versus previous); these are observations, not a
general speedup claim. Streaming memory behavior is unchanged.

## Automated coverage

`test_munging_filter.py` covers all supported missing tokens, reference-only
placeholder rows, exact text and source preservation, zero and signed Z,
inclusive threshold boundaries, exclusions, summary counts, invalid values,
empty/all-missing/all-filtered files, unchanged previous outputs on failure,
temporary-file cleanup, CRLF and absent final newlines, whitespace input, and
parallel manifest order. `test_restart.py` exercises the managed CLI with actual
filtering and compilation plus a mocked LDSC command, checking missing counts
and unchanged results on restart.

The targeted 22 filter/restart tests passed. The full suite ran 151 tests:
148 passed and three optional native-CLI environment integrations were skipped.
Validation used Python 3.12.14, pandas 2.3.3, NumPy 1.26.4, Polars 0.20.31 and
R 4.4.2. Reproduce the targeted checks with:

```bash
PYTHONPATH=src:tests python -m unittest test_munging_filter test_restart -v
```

Private validation copies, the unmodified native parser and its source hash,
`validate_real_inputs.py`, and `summary.json` are stored in the local sibling
folder `missing_z_validation_20261002/`. They are not committed.
