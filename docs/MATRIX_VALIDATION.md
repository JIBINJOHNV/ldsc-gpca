# Separate PC1 and CTI validation

[README](../README.md) · [Python GPCA](PYTHON_GPCA.md) · [GenomicSEM GPCA](GENOMICSEM_GPCA.md)

## Choose the output you need

The same options apply to `ldsc-gpca gpca`, `ldsc-gpca genomicsem gpca`, and the
GPCA stage of `ldsc-gpca pipeline` with either LDSC backend.

| Mode | PC1 checks | CTI checks | GWAMA and SNP export |
| --- | --- | --- | --- |
| Default full run | Must pass | Must pass | Run after checks |
| `--validate_only` | Must pass | Must pass | Skipped |
| `--pc1_only` | Must pass | Reported; CTI numerical failure is nonblocking | Skipped |

Use either `--validate_only` or `--pc1_only`, not both. Neither mode prepares or
checks per-SNP GWAMA tables in standalone GPCA. Pipeline still performs its
upstream VCF extraction, munging and LDSC; these flags do not mean a dry run.
Its separate GWAMA preparation is skipped unless `--write_munge_inputs` requests it.

`--pc1_only` does not bypass the LDSC input contract: complete selected pairs,
valid required numerical estimates, self-pairs, consistent duplicates and
matching trait order are still required. It permits a non-positive-definite or
ill-conditioned CTI, not malformed or missing input. A successful PC1-only run
must not be described as GWAMA-ready unless its separate CTI assessment passes.

## Run checks on completed LDSC results

```bash
# Assess both matrices; return an error if a required check fails.
ldsc-gpca gpca \
  --input selected_traits.csv \
  --ldsc_results ldsc_results.csv \
  --outdir both_matrix_checks \
  --validate_only

# Compute PC1 without requiring CTI to pass the downstream numerical gate.
ldsc-gpca gpca \
  --input selected_traits.csv \
  --ldsc_results ldsc_results.csv \
  --outdir pc1_only_results \
  --pc1_only
```

For GenomicSEM, use `ldsc-gpca genomicsem gpca` and a final
`genomicsem_LDSC.RData` containing `LDSCoutput`. The Python CSV supplies marginal
SEs, not full sampling covariance; it is never converted to a fabricated
GenomicSEM object. Original inputs are not modified.

## What each assessment means

**PC1:** construct the requested correlation or covariance matrix after applying
the selected normalization, then compute its symmetric eigendecomposition.
PC1 requires a finite positive leading eigenvalue and finite, nonzero loadings
in manifest order. The report also gives PC2 and the PC1–PC2 eigenvalue gap.
A positive PC1 eigenvalue establishes computability, not scientific reliability.
A small gap can warrant a stability investigation; this package does not claim a
universal scientifically validated gap cutoff or calculate sampling SEs for PC1
from marginal LDSC SEs.

`--negative_eigen_action warn` (default) reports negative eigenvalues anywhere
in the genetic matrix and permits an otherwise valid PC1. `error` blocks the
run for substantive negative eigenvalues even if PC1 is computable. The
combined report distinguishes `PC1_Status` from `PCA_Matrix_Policy_Pass`.
This option never deletes traits, repairs CTI or clamps correlations. The
existing relative `--matrix_eigen_tolerance` controls this genetic-matrix
warning/error threshold; it is not the CTI threshold.

**CTI:** construct the intercept/error covariance matrix from the self-pair
intercepts and cross-trait intercepts. Report its eigenvalues, Cholesky result
and condition number. The shared numerical gate requires successful Cholesky,
minimum eigenvalue above `max(1e-8, 100*k*machine_epsilon*max(abs(eigenvalues)))`,
and condition number at most `1e8`. These are explicit numerical policy choices,
not biological quality cutoffs. They can reject a nearly singular matrix even
when Cholesky alone succeeds. A failed CTI is not evidence that PC1 itself failed.

These are **trait-level** matrices: exclusion removes a protein/trait and its
entire row and column, not individual participants or selected SNP pairs.

## Optional exploratory CTI selection

The default `--cti_action error` does not remove traits for CTI failure.
`--failed_ldsc_action drop_traits` remains a separate policy for missing or
invalid estimates; it does not enable numerical CTI selection.

To explicitly request exploratory selection:

```bash
ldsc-gpca gpca \
  --input selected_traits.csv \
  --ldsc_results ldsc_results.csv \
  --outdir exploratory_cti_subset \
  --validate_only \
  --cti_action explore_drop \
  --max_cti_drop_fraction 0.01
```

Here `0.01` permits at most `floor(0.01 * k)` exclusions from the traits entering
CTI assessment, after any earlier missing/failed-estimate exclusions. It is an
example user-selected limit, not a recommended universal threshold. The default
is zero; `explore_drop` requires an explicit positive fraction below 1. A small
fraction can permit zero exclusions on a small dataset. At least two traits must
remain. A nonzero limit with `cti_action error`, or `explore_drop` with
`--pc1_only`, is rejected rather than silently ignored.

At each iteration:

1. Assess PC1 and CTI separately and save their results.
2. Stop on invalid PC1 or a selected PCA-matrix error policy. Neither is used
   as a reason to delete a trait automatically.
3. If CTI fails, score squared eigenvector participation across its problematic
   directions, weighted by their eigenvalue deficits. The target eigenvalue is
   the larger of the CTI tolerance and `maximum_eigenvalue / 1e8`.
4. Remove the highest-scoring trait, with exact ties resolved in manifest order.
5. Subset both matrices using that same retained trait order and recompute both
   assessments. Stop when CTI passes or the exclusion limit is exhausted.

This greedy heuristic does not guarantee the smallest exclusion set or the
largest improvement at each step. It does not prove that excluded traits have
bad data. No matrix repair, covariance shrinkage, rg clipping or missing-pair
imputation is performed. Full GenomicSEM sampling-covariance matrices are subset
consistently with the retained traits; they are not reconstructed from marginal SEs.

PC1 loadings before/after selection are compared on retained traits after
aligning their overall sign. The report gives per-trait changes, RMS change and
maximum absolute change. These are descriptive sensitivity metrics, not
significance tests. No automatic scientific-acceptability threshold is imposed:
review the sensitivity before adopting the subset for a main analysis. CTI
selection does not diagnose upstream SNP QC, sample size or weighting differences
from a matrix file alone, and cannot guarantee that future datasets will pass.

## Outputs and failure semantics

| File | Meaning |
| --- | --- |
| `GenomicPCA_Validation_Settings.rds` | Exact parsed GPCA settings, including normalization, execution mode and exclusion policy. |
| `GenomicPCA_Matrix_Validation.csv` | Latest assessment: trait count, PC1 status/error, PC1/PC2 eigenvalues and gap, genetic-matrix policy result, CTI metrics, and `GWAMA_Matrix_Eligible`. This eligibility checks matrices, not SNP files or GWAMA results. |
| `GenomicPCA_Matrix_Validation_History.csv` | Initial and every post-exclusion assessment. |
| `GenomicPCA_CTI_Excluded_Traits.csv` | Exclusion step, exact trait, score, eigenvalue before/after, and exploratory reason. Header-only if none. |
| `GenomicPCA_Assessed_Traits.csv`, `GenomicPCA_CTI_Assessed.csv` | Last assessed subset, including a failed attempt. Not a promise of usability. |
| `Python_LDSC_Retained_Traits.csv` | Earlier Python missing/invalid-estimate QC checkpoint; it can include traits later excluded for CTI. Use the final GenomicPCA manifest below for subsequent analysis. |
| `GenomicPCA_Retained_Traits.csv` | Manifest in exact final order, published only after the requested matrix mode passes. In PC1-only mode CTI can still fail. |
| `GenomicPCA_PC1_Selection_Sensitivity.csv`, `GenomicPCA_PC1_Selection_Sensitivity_Summary.csv` | Sign-aligned changes and aggregate sensitivity. Without exclusions, changes are zero. |
| `Python_LDSC_Retained_Results.csv` (or `.tsv`/`.txt`) | Python backend only: complete source rows restricted to final retained traits, preserving original columns and numerical text. Matches the final manifest. Gzipped input is written decompressed. |
| `GenomicSEM_LDSC_Used.RData` | GenomicSEM backend: final consistently subset LDSC object, including its real sampling covariance. |
| Existing `GenomicPCA_*_Used.csv` and PC1 reports | Final matrices, loadings and eigenvalue diagnostics. |
| `GWAMA_Run_Status.csv` | Explicitly identifies `not_run_pc1_only` or `not_run_validate_only` when appropriate. |

Use a **fresh output directory**. An existing matrix assessment is not overwritten. Reports written during a failed selection are
not usable exports. Failure returns a nonzero exit code; no final retained
manifest or retained LDSC export is published for that failed attempt. Fatal
structural input errors can occur before matrix construction and therefore
before either matrix assessment exists. Preserve their original error logs.

For a later full GWAMA run, supply the final retained manifest and retained LDSC
results with matching per-trait SNP tables. Repeat the analysis settings, such
as normalization and PCA matrix choice, and omit the PC1/validation-only flag.
Do not pass an invalid CTI from a PC1-only run to GWAMA.

Standalone full commands needing automatic VCF preparation first run combined
matrix validation into `outdir/matrix_preflight/`. A failure blocks preparation.
On success, the ordinary full run revalidates after preparation; this is not a
cached/resumed analysis. Preparation still uses the supplied VCF manifest, so
use an already-retained manifest to avoid preparing excluded traits.
