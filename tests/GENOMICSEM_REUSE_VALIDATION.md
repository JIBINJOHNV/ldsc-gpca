# Preparation publication and GenomicSEM reuse validation

Validated on 2026-10-04 using Python 3.12.14, real bcftools, R 4.3.2 and the
**complete installed GenomicSEM 0.0.5 package**, with installed `RemoteSha`
`6b65ca5db39fdade08b0d811477be1cdd57b5039` verified. This replaces the earlier
munging-only harness limitation for the synthetic scenarios below. Neither the
upstream estimator nor its munging algorithm was changed.

## Changes and policy

- Preparation records completion atomically only after directory publication,
  atomic manifest update and artifact hashing. Failures retain raw inputs/logs;
  any already-published munged files remain incomplete and are rejected on reuse.
- New shared-preparation sidecars record both prevalences, N source/override/
  convention and a schema version. Completion covers the prepared manifest,
  munged files and sidecars. Python checks these new bundles while retaining its
  existing legacy-file policy.
- GenomicSEM preparation additionally records installed package/R versions,
  pinned source commit and a copied allele reference with its checksum.
- GenomicSEM reuse requires completed provenance, matching checksums/metadata,
  the same reference checksum across bundles, unique SNPs, valid N/Z and exact
  reference A1/A2 order. Swaps and either strand-complement orientation fail;
  there is no implicit realignment or exclusion. Older unverified GenomicSEM
  preparation outputs must be re-prepared from raw inputs.
- The input audit records accepted trait SNP counts/order hashes, or the first
  invalid trait/row/SNP and reason. It is not an exhaustive list of all invalid
  records after a fatal error.

## Tests and observations

| Test | Expected | Observed |
| --- | --- | --- |
| Publication failures at directory rename, manifest replacement and completion replacement | Never publish successful status; retain raw bytes/logs; reject incomplete reuse | Passed, including preservation of the previous manifest and cleanup of temporary writes |
| Existing aligned inputs, zero/negative finite Z, variable SNP availability, literal trait names, relocated bundle | Preserve bytes, names/order and signed products | Passed |
| Swapped, complemented, swapped-complemented, incompatible pairs; duplicate/unknown SNP IDs; missing/invalid N/Z | Fatal input error identifying the affected SNP; no regression or input edits | Passed |
| Missing/malformed sidecars, wrong backend/runtime, file/sidecar/reference/manifest checksums, N and prevalence contradictions | Fatal provenance error | Passed |
| Full pinned munger, two traits and 2,400 reference SNPs | A retains 2,400 SNPs; B retains 1,200 same/swapped raw pairs and removes complemented/incompatible raw pairs | Exact expected SNP sets; all retained alleles match reference order |
| Independent signed-product formula using actual prepared P text | `sign(BETA) * sqrt(chi2.isf(P,1))`, with a sign reversal for swapped raw alleles | Products agree within `rtol=1e-12, atol=1e-12`; no regression chi-square exclusions at 80 |
| Real two-pass wrapper versus direct installed `GenomicSEM::ldsc` | Unchanged numerical estimates/intercepts, S/I/V/S_Stand/V_Stand and SEs from V/V_Stand | Exact numeric equality; wrapper's existing matrix-name additions are checked separately from numerical values |
| Altering actual full-package outputs or their N/prevalence/checksum metadata | Reuse rejects them before estimation | Passed |

The regression fixture uses seed 731, one synthetic LD-score chromosome, N=20,000,
M=100,000 and 20 jackknife blocks. Its observed rounded h2 estimates are 0.1387
and 0.2368, genetic covariance 0.0717 and rg 0.3953. The comparison asserts exact
stored numerical values, not these rounded log summaries. Raw VCF FORMAT values
are subject to bcftools query precision, so the independent equation uses the
actual prepared P values rather than the higher-precision simulation source.

The complete Python discovery run executed **285 tests: 280 passed, 5 skipped**.
The skipped tests require a separately configured optional Python LDSC runtime;
the full GenomicSEM integration was executed. After two additional boundary tests
and expanded full-package rejection checks were added, the final targeted run
executed **14 tests, all passed**. The existing real Python preparation test also
passed, including rejection of a failed completion marker before Python regression.

## Reproduction

Install the pinned full package using `scripts/install_genomicsem.R` in a suitable
environment. Point `LDSC_GPCA_TEST_RSCRIPT` to that Rscript, or a wrapper selecting
its isolated library; place bcftools on PATH. From the repository:

```bash
PYTHONPATH=src:tests python -m unittest \
  test_preparation_publication test_prepare_modes \
  test_genomicsem_inputs test_genomicsem_reuse_integration

PYTHONPATH=src:tests python -m unittest discover -s tests -p 'test_*.py'
```

Boundary fixtures in `genomicsem_fixtures.py` use explicitly synthetic provenance.
The full-package integration generates real VCFs, calls actual preparation and
both real LDSC passes, and compares against a direct package call. There is no
injected estimator in that test. This validates the specified synthetic cases;
it does not establish appropriate scientific N/INFO definitions for every GWAS
design or substitute for validation on a production cohort.

Primary algorithm references:
[pinned GenomicSEM munging](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/munge_main.R),
[pinned GenomicSEM LDSC](https://github.com/GenomicSEM/GenomicSEM/blob/6b65ca5db39fdade08b0d811477be1cdd57b5039/R/ldsc.R),
[GWAS-VCF field definitions](https://github.com/MRCIEU/gwas-vcf-specification).
