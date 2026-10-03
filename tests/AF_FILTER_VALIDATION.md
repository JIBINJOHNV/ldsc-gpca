# Missing INFO allele-frequency filtering (review item 8.4)

Validated on 2026-10-03 with Python 3.12, bcftools/htslib 1.23, Bash, awk,
Polars 0.20.31 and R 4.4.2 for the existing R regression checks.

The extraction predicate now explicitly excludes missing INFO/AF as well as
missing INFO/EUR. INFO frequencies are never replaced with FORMAT/AF. Missing
header declarations remain errors. Missing-value expressions follow the
[official bcftools manual](https://samtools.github.io/bcftools/bcftools.html#EXPRESSIONS).

A shared streaming counter writes `munge_input/{traitname}_AF_Filter_QC.csv`
after successful extraction. It counts records after optional MHC selection,
before any other extraction filters. Missing AF and EUR counts may overlap;
`excluded_missing_either` counts their union. Reruns clear a stale audit before
starting and publish a replacement only after the pipeline succeeds. No second
VCF read or original Python LDSC source modification is needed.

## Validation

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:tests \
LDSC_GPCA_TEST_RSCRIPT=/Users/JJOHN41/miniconda3/bin/Rscript \
../.venv-ldsc-normalization/bin/python -m unittest \
  test_extraction_af test_ldsc_runtime test_worker_retries \
  test_munging_filter test_interface_aliases
```

Result: **44 tests passed**, no skips (55.751 seconds).

The new tests execute real bcftools extraction for both quantitative and
population-prevalence traits, rather than mocking the filtering commands.

| Check | Expected and observed |
| --- | --- |
| Mixed 15-record VCF, MHC exclusion enabled | 14 records reach AF QC; 6 lack INFO/AF, 3 lack INFO/EUR, 8 lack either; 3 records survive all filters |
| Missing record tags, explicit dots, partial missing vectors, both missing | Excluded and counted; a different INFO tag named EAF does not substitute for AF |
| FORMAT/AF disagrees with INFO/AF | AF-difference filter uses INFO/AF; FORMAT/AF still supplies the existing MAF check |
| Difference exactly equal to cutoff, zero INFO frequencies | Retained when the other filters pass |
| Valid inputs under previous and new predicates | Extracted TSVs are byte-for-byte identical in both routes, with palindromic filtering enabled and disabled |
| Binary counts and prevalence | Retained rows have N_TOTAL=1000 and derived sample prevalence=0.1 |
| INFO/AF or INFO/EUR header missing | Both extraction routes raise a bcftools error; stale audit removed and temporary counts cleaned up |
| Empty input and reruns | Empty audit counts are zero; reruns replace counts without accumulation |
| Missing header through parallel extraction | Two failed attempts recorded; pipeline stops before prevalence metadata/downstream work |

The fixtures are synthetic VCFs designed to cover these cases, including paths
with spaces. This validation did not rerun a whole-genome LDSC/PCA/GWAMA analysis
or benchmark production throughput. Existing munged files supplied through
`--ldsc_only` are not retroactively filtered; extraction and munging must be
rerun to apply the new rule to them.
