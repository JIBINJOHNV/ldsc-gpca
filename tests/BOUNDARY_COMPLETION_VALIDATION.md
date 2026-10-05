# Boundary fixes and preparation parity

Validated on 2026-10-04 with Python 3.12.14, bcftools/htslib 1.23,
R 4.3.2 (full pinned GenomicSEM library), and R 4.4.2 for focused checks.
The installed GenomicSEM commit is
`6b65ca5db39fdade08b0d811477be1cdd57b5039` (version 0.0.5).

## Scope and fixed contracts

- Every VCF route requires exactly one GWAS sample. Quantitative and binary
  Python extraction reject zero/multiple samples before creating or opening
  extraction tables or their AF-QC sidecars. Shared preparation and GenomicSEM
  conversion use the same guard. No sample-selection option was added.
- `vcf_common.py` contains the shared schema, sample check, row transformation,
  identifier checks and original-record audit helpers. GPCA writing and
  combined-mode orchestration stay in `prepare.py`; LDSC writing, metadata and
  munging stay in `prepare_ldsc.py`. The reverse import was removed.
- The existing R numeric guard now runs during reading and before heritability
  normalization. Malformed text identifies its field and original file row;
  direct in-memory calls identify their input row. Both policies reject it.
  Missing numeric tokens remain eligible for explicit trait exclusion.
- GWAMA uses the fixed **error** policy for unusable output, outside the vendor
  function. Raw output is retained, with per-output `GWAMA_QC_Issues.csv` and
  `GWAMA_QC_Summary.csv` suffixes. Checks cover zero available weight, BETA/Z/P
  finiteness, positive finite SE/N and P bounds. P=0 is accepted. Each weight
  is considered separately; a zero signed sum is not a zero-weight test.
- The exporter independently validates raw associations before N overrides,
  publication or archival. Its per-output suffixes are
  `GWAMA_Export_QC_Issues.csv` and `GWAMA_Export_QC_Summary.csv`.
- Existing worker retries/status capture, manifest rules, backend defaults,
  runtime/provenance modules and distinct scientific filters are retained.
  No generic pipeline framework or broad R QC reorganization was introduced.

## Coverage of the eight findings

| Finding | Regression coverage / outcome |
| --- | --- |
| 1. GenomicSEM reuse allele and metadata checks | Existing `test_genomicsem_inputs.py` and full installed-package `test_genomicsem_reuse_integration.py` pass. Same/swapped raw alleles are munged consistently; complemented/incompatible raw pairs follow the pinned munger's documented exclusion rule. Reused files must already match the verified reference, with valid sidecars/checksums/N/prevalences. See `GENOMICSEM_REUSE_VALIDATION.md`. |
| 2. Conflicting repeated self estimates | Existing Python/R `test_self_duplicates` pass for each self-only field, all policies, tolerances, missing/nonfinite values, negative/positive mixtures and valid reverse-row differences. |
| 3. Implicit first VCF sample | New `test_vcf_boundaries.py` passes for zero/two samples, both Python extraction routes and all preparation routes. Existing files/sidecars remain byte-identical; fresh invalid jobs publish no tables; failed attempts are recorded. |
| 4. Supplied munged-directory overwrite | Existing dispatch/guard tests pass: invalid fresh-mode combination fails without modifying supplied inputs; explicit reuse remains supported. |
| 5. Undefined GWAMA associations | New actual-vendor tests identify the C-only zero-loading SNP and preserve its raw NA Z/P plus positive N_eff=1,000. Both backend worker paths fail after their existing two attempts. Export numerical tests prevent publication/archival and reject attempts to conceal bad N through an override. |
| 6. Reverse-pair compilation conflicts | Existing common-field tests reject conflicts under every policy, while valid missing-estimate reporting remains covered. |
| 7. Malformed h2 text | New `test_numeric_input` and updated literal-ID tests fail at the reader/normalization boundary under both policies, with field/row errors. Valid numeric strings, missing tokens, observed/liability/mixed scales and matrix values pass. |
| 8. Literal trait identifiers | Existing literal-ID tests preserve quoted/unquoted NA, NaN, numeric-looking and ordinary names/order. Empty identifiers still fail. |

## Numerical and byte comparisons

`test_gwama_output_qc.R` calls the actual bundled GWAMA implementation using
three traits, signed loadings `(0.8,-0.6,0.2)`, N of 1,000/2,000/500 and CTI
with diagonal 1 and off-diagonal 0.1. It independently checks
`sum(w*Z)/sqrt(t(w)%*%CTI%*%w)`, where `w=sqrt(N)*loading`, restricted to
available, aligned inputs. Ordinary and partial availability, a swapped
allele pair and P underflow pass. The tolerances are absolute 1e-10 for Z/BETA,
1e-12 for SE/P and 1e-9 for effective N. All expected Direction strings match.
A separate block-diagonal genetic correlation (A/B=.5, C independent) yields
C's exact zero loading and the expected audited failure for a C-only SNP.
No estimator, vendor GWAMA source, PCA or covariance algorithm was changed.

`prepare_contract_fixture.py` was run against pre-change commit `056d076` and
the updated source, using a two-trait, 22-autosome fixture with duplicates,
invalid SE, missing SI, P-floor records and a manifest N override. All **91**
products matched: default whole-genome, split, legacy raw export, LDSC-only,
combined and combined-with-N-override outputs. Output directory strings were
normalized; worker-attempt files containing timings were excluded. SNP table
bytes, QC, settings and prepared manifests otherwise matched exactly.

The final discovery run executed **294 tests: 289 passed, 5 skipped** in
213.491 seconds. The five skipped optional Python LDSC integrations require
explicit runtime/script configuration. Full installed GenomicSEM regression
and bundled GWAMA checks ran. Generated Python/R CLI help parity and
`git diff --check` also passed.

## Reproduce

With bcftools on PATH and `LDSC_GPCA_TEST_RSCRIPT` pointing to an Rscript using
the full pinned GenomicSEM installation:

```bash
PYTHONPATH=src:tests python -m unittest discover -s tests -p 'test_*.py'
python scripts/update_cli_help.py --check
```

For preparation parity, use a fresh shared fixture directory and point
`PYTHONPATH` to each source snapshot in turn:

```bash
PYTHONPATH=/path/to/baseline/src python tests/prepare_contract_fixture.py /tmp/parity before
PYTHONPATH=src python tests/prepare_contract_fixture.py /tmp/parity after
cmp /tmp/parity/before.json /tmp/parity/after.json
```

These are synthetic boundary/integration checks, including actual pinned
GenomicSEM regression and actual bundled GWAMA. They do not establish suitable
N definitions for every study design or constitute production-cohort validation.
