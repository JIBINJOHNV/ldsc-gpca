# Python LDSC output file reference

[README](../README.md) · [Python LDSC guide](PYTHON_LDSC.md) · [Final result columns and real example](PYTHON_LDSC.md#results-and-next-step)

Use this reference when you need to understand an intermediate file, investigate
a failed run, or preserve enough information to restart. The main guide shows
separate folder layouts for [VCF input](PYTHON_LDSC.md#start-from-vcf-files) and
[munged input](PYTHON_LDSC.md#start-from-munged-files).

Paths below are relative to the standalone command's `--outdir`. For a complete
pipeline, these files are inside its `ldsc/` stage directory. CSV files are
comma-delimited with headers; `.tsv.gz` files are gzip-compressed tabular text.
The extracted quantitative `.tsv` is an exception: it uses spaces, as described
below. Conditional reports appear only when their stage is reached.

## In this reference

- [Which files should I open first?](#which-files-should-i-open-first)
- [Final estimates and collection status](#final-estimates-and-collection-status)
- [Pair and trait QC reports](#pair-and-trait-qc-reports)
- [Trait removal reports](#trait-removal-reports)
- [Resolved sample size and prevalence](#resolved-sample-size-and-prevalence)
- [VCF extraction files](#vcf-extraction-files)
- [Munged files and provenance](#munged-files-and-provenance)
- [Chi-square filtering files](#chi-square-filtering-files)
- [Regression batches and restart files](#regression-batches-and-restart-files)
- [Worker attempts, runtime and errors](#worker-attempts-runtime-and-errors)
- [Real-data verification](#real-data-verification)

## Which files should I open first?

| Your question | File to inspect |
| --- | --- |
| Did collection finish with usable estimates? | `LDSC_Compilation_Status.csv`; check `Status` and `Result_File`. |
| What estimates should I give to GPCA? | `ldsc_results.csv`, after checking QC and pair coverage. |
| Why did a particular comparison fail or warn? | `LDSC_Pair_Status.csv`; follow its source file/row to the numerical batch export. |
| Which traits have failed self-pairs or comparisons? | `LDSC_Trait_Status.csv`. |
| Which traits remain after `drop_traits`? | `LDSC_Retained_Traits.csv` and `LDSC_Dropped_Traits.csv`. |
| Did a command fail, retry or get reused? | `ldsc_results/LDSC_Batch_Status.csv`, worker attempt reports and available error logs. |
| Which SNPs did the extra Z² filter remove? | `LDSC_ChiSquare_Filter_Summary.csv` and `LDSC_ChiSquare_Excluded_Variants.tsv.gz`. |
| What N convention and filter settings were used? | `LDSC_Trait_Prevalence_Metadata.csv`, the munged `.prevalence.json` files and restart checkpoints. |

## Final estimates and collection status

| File | When written | Meaning |
| --- | --- | --- |
| `ldsc_results.csv` | Successful collection, including a valid retained subset under `drop_traits`. | Final pair/self estimates. See the [complete column definitions](PYTHON_LDSC.md#the-main-output-ldsc_resultscsv). |
| `ldsc_results_diagnostic.csv` | Collection reaches numerical QC, on successful as well as failed-estimate runs. | All collected estimates before trait selection. Same result columns; can contain failed or unusable values. |
| `LDSC_Compilation_Status.csv` | Initialized before analysis, then updated by collection. | Whether a final result was produced and which result policy was used. |

The compilation status has these base columns:

| Column | Meaning |
| --- | --- |
| `Status` | Current collection outcome; meanings below. |
| `Action` | Requested `error`, `report` or `drop_traits` policy. |
| `Result_File` | Published result path, diagnostic path, or empty when no result was published. |
| `Failed_Rows` | Number of non-valid result-status rows; empty when collection has not reached that check. |

`drop_traits` also writes `Retained_Traits`, `Excluded_Traits` and
`Retained_Manifest`. `Error` is added for structural failures or insufficient
retained traits.

| Status | How to read it |
| --- | --- |
| `incomplete` | Analysis/collection has not finished. An earlier stage may have stopped; inspect its logs. |
| `compiled` | No numerical estimate failures in the published selection. Diagnostic warnings can still exist. |
| `compiled_with_trait_exclusions` | A usable retained subset was published; use its retained manifest. |
| `estimation_failures` | Invalid estimates were found; `Result_File` points to diagnostics. No usable final CSV is published. |
| `insufficient_traits` | Trait removal left fewer than two traits. Removal reports exist, but no usable final CSV is published. |
| `structural_failure` | Collection failed because of malformed, inconsistent or otherwise structurally invalid input. Read `Error`. |

In `report` mode, a zero process exit code can accompany `estimation_failures`.
Read the status rather than treating file existence or command completion as
scientific validation.

When collection is prepared again, existing main/diagnostic result CSVs,
pair/trait reports and retained/dropped manifests are renamed with a
`.previous-<uuid>` suffix. These are backups. The unsuffixed files and current
compilation status describe the current attempt.

## Pair and trait QC reports

### `LDSC_Pair_Status.csv`

One record per collected source result row. Estimates can be valid while still
carrying a warning, such as a correlation outside the usual range.

| Column | Meaning |
| --- | --- |
| `p1`, `p2` | Trait names for this directed comparison. |
| `Source_File` | Numerical batch file from which the row came. |
| `Source_Row` | Source line number with the header counted as line 1. |
| `Status` | `valid` or `failed_estimate` for an exported estimate. |
| `Reason` | Failed numerical checks, separated by semicolons; empty for a valid row. |
| `Warning` | Diagnostic notes, including out-of-range rg, weak self h2/SE or the permitted self-pair zero-SE exception. |

Required finite values, positive SEs, positive h2 and P values in `[0,1]` are
checked on the applicable scale. The narrow self-correlation exception requires
`p1 == p2`, finite `rg` within `0.01` of 1, `se=0`, `z=+Inf` and `p=0` together.
It preserves the native values and records a warning. Other required SEs remain
strictly positive.

### `LDSC_Trait_Status.csv`

One record per manifest trait, preserving order. Heritabilities and intercepts
here come from each trait's **self-pair**, rather than a between-trait row.

| Column | Meaning |
| --- | --- |
| `Manifest_Order` | One-based position in the input manifest. |
| `Trait` | Exact trait name. |
| `Self_Status` | `valid`, `failed_estimate` or `not_requested` when no self-comparison was requested. |
| `Self_Reason` | Reasons the self-pair failed. |
| `Failed_Pair_Rows` | Failed directed result rows involving this trait. Each self-row counts once; A–B and B–A are separate rows. |
| `h2_obs`, `h2_obs_se` | Self-pair unconverted h2 and SE, when available. |
| `h2_liab`, `h2_liab_se` | Self-pair liability h2 and SE, when available. |
| `h2_int`, `h2_int_se` | Self-pair LDSC intercept and SE. |
| `Retained`, `Exclusion_Reason` | Added only with `drop_traits`; selection flag and removal explanation. |

A valid numerical result does not guarantee strong heritability or a positive
definite downstream matrix. GPCA performs its own matrix and input checks.

## Trait removal reports

These files are written only with `--result_failure_action drop_traits`, after
collection reaches trait selection. They may also be written when selection
ultimately leaves too few traits to continue.

**`LDSC_Retained_Traits.csv`** is the original input manifest restricted to
retained traits, in its original order and with its original column names.
Its columns therefore depend on the supplied manifest; it is not a fixed
numerical-results schema. Pass it to GPCA with the final retained results.

**`LDSC_Dropped_Traits.csv`** explains each removal:

| Column | Meaning |
| --- | --- |
| `Exclusion_Step` | Order in which the selection procedure removed traits. |
| `Manifest_Order`, `Trait` | Original manifest position and exact name. |
| `Reason` | Why the trait was removed. |
| `Failed_Pairs_At_Exclusion` | Number of failed pair connections at that removal step. |
| `Self_H2`, `Self_H2_SE`, `H2_Z` | Available self h2, its SE and their ratio used in selection. |

The selection removes failed self-pairs first, then removes traits to eliminate
failed pair connections. It is a deterministic greedy selection, not a claim
that the largest possible subset was found. See the
[selection rules](REFERENCE.md#python-ldsc-outputs) before using this mode.

## Resolved sample size and prevalence

**`LDSC_Trait_Prevalence_Metadata.csv`** is written by both input routes. It
retains manifest fields with two internal-name changes:
`traitname` becomes `gwas_name`, and `population_prevalence` becomes
`pop_prevalence`. Additional manifest columns are retained, except conflicting
internal-name columns.

| Column | Meaning |
| --- | --- |
| `gwas_name` | Trait name from `traitname`. |
| `ref` | Whether this trait is an LDSC reference (`yes`). |
| `pop_prevalence` | Resolved population prevalence, or empty. |
| `sample_prevalence` | Supplied, inferred or restored case fraction; empty when unused. |
| `vcf_files` | Source VCF path when present in the manifest; reuse does not read this path. |
| `sample_size_column` | `NEF` without population prevalence; `N_TOTAL` when population prevalence is supplied. |
| `sample_prevalence_source` | `not_used`, `provided`, `median_case_fraction` or `saved_metadata`, according to how the case fraction was resolved. |

This file records the resolved N convention, not a per-SNP list of N values.
Read N itself from the munged file.

## VCF extraction files

These are created only when extraction runs. They are under `munge_input/`.

### `{trait}_mungeinput.tsv`

Association statistics extracted after the requested VCF filters, before native
munging QC. Quantitative/NEF extraction currently writes **space-delimited**
text despite the `.tsv` extension. Binary total-N extraction writes tabs.

| Column | Source / meaning |
| --- | --- |
| `CHROM`, `POS` | VCF chromosome and base-pair position. |
| `ID` | VCF variant ID; used as the SNP identifier for munging. |
| `REF`, `ALT` | Reference and alternate alleles. Munging receives ALT as A1 and REF as A2. |
| `EZ` | Signed association Z from FORMAT/EZ; supplies effect direction to munging. |
| `P` | `10^(-LP)` from FORMAT/LP. |
| `AF` | FORMAT/AF; the frequency used by native munging. |
| `NEF` | FORMAT/NEF; chosen for N when population prevalence is absent. |
| `N_CASES`, `N_CONTROLS`, `N_TOTAL` | Binary extraction only: FORMAT/NC, FORMAT/NCO and their sum. Total N is chosen for munging. |

Column order is `CHROM ID POS REF ALT EZ P AF NEF`, followed by the three
count columns for binary extraction. The `prepare --p_min` option does not
apply here. Extremely small probabilities can underflow during conversion;
this is distinct from the full-precision export of the later munged Z values.

### `{trait}_AF_Filter_QC.csv`

Counts missing frequency annotations after optional MHC selection and before
the subsequent frequency-difference, SI/MAF and palindrome filters.

| Column | Meaning |
| --- | --- |
| `records_after_mhc` | Records entering this stage after optional MHC exclusion. |
| `missing_info_af` | Records with missing INFO/AF. |
| `missing_info_eur` | Records with missing INFO/EUR. |
| `excluded_missing_either` | Records missing either or both fields. Counted once per record. |

The two missing-field counts can overlap. This file does not count every record
removed by all extraction filters; use the extracted table and munging log to
follow subsequent stages.

## Munged files and provenance

New munged outputs go to `<outdir>/ldsc_input/`, or to `--munged_dir` when
specified in VCF mode. With `--ldsc_only`, these are **existing inputs** read
from that directory, not new outputs.

| File | Meaning |
| --- | --- |
| `{trait}.sumstats.gz` | Gzip-compressed, tab-delimited native munged summary statistics. |
| `{trait}.log` | Native munging commands, filtering messages and counts. Readable audit, not the final estimates. |
| `{trait}.prevalence.json` | Provenance sidecar used to check reuse. Preserve it alongside the sumstats file. |

The managed sumstats header is `SNP N Z A1 A2`, with tab separators:

| Column | Meaning |
| --- | --- |
| `SNP` | Variant ID retained by munging. |
| `N` | Per-SNP N from the selected `NEF` or `N_TOTAL` column. |
| `Z` | Native munging Z reconstructed from P and the sign of EZ. Exported with `%.17g` precision, without three-decimal rounding. |
| `A1`, `A2` | Alleles retained by native munging; supplied as ALT and REF respectively by this extraction route. |

The sidecar has keys `sample_size_column`, `sha256`, `filters`, `float_format`
and `sample_prevalence`. `sha256` covers the compressed sumstats file;
`float_format` describes the export. `filters` records MHC, INFO, extraction
MAF, munging MAF, frequency-difference and palindrome settings. Reuse checks
these settings when present; passing different settings does not re-filter an
existing munged file. Missing legacy filter metadata yields a warning.
Total-N reuse requires a sidecar; legacy NEF reuse without one is permitted
with a warning.

## Chi-square filtering files

Created for either input route only when `--chisq_max` is supplied. The filtered
copy keeps the source header and retained records, and becomes the regression
input. Source munged files are preserved.

**`LDSC_ChiSquare_Filter_Summary.csv`** has one row per trait:

| Column | Meaning |
| --- | --- |
| `gwas_name` | Trait name. |
| `chisq_max` | Numeric Z² cutoff actually applied to this trait. |
| `threshold_mode` | `fixed` or `auto`. |
| `maximum_matched_n` | Maximum N after complete-row matching to LD references and weights; populated for `auto`. |
| `complete_ld_matched_rows` | Number of complete matched rows used to determine automatic N; populated for `auto`. |
| `variants_before` | Source munged record count. |
| `variants_removed` | Records with Z² greater than the cutoff. |
| `variants_after` | Records written to the filtered copy, including missing-Z rows. |
| `variants_missing_z` | Rows retained for native missing-Z handling. |
| `maximum_chisq` | Largest finite Z² encountered before applying the cutoff. |
| `source_file` | Original munged file path. |
| `filtered_file` | Path to `ldsc_input_chisq_filtered/{trait}.sumstats.gz`. |

Both `ldsc_input_chisq_filtered/{trait}.chisq_excluded.tsv.gz` and the combined
`LDSC_ChiSquare_Excluded_Variants.tsv.gz` have:

| Column | Meaning |
| --- | --- |
| `gwas_name` | Trait whose record was removed. |
| `SNP` | Removed variant ID. |
| `Z` | Its original Z value. |
| `CHISQ` | Z², which exceeded that trait's threshold. |

They can contain only a header when no SNPs were excluded. A missing Z is
counted separately, not recorded as a large-Z exclusion.

## Regression batches and restart files

Files are written under `ldsc_results/`. A batch prefix has the form
`{reference_trait}_{batch_size}_ldsc_{part}`. The part index starts at zero;
batch size depends on the number of traits, references and requested workers.
Do not assume that every run uses the same batch filename or size.

| File | Meaning |
| --- | --- |
| `{prefix}.results.csv` | Native numerical result table used for collection. Original `p1,p2` entries are sumstats paths; collection resolves them to trait names. |
| `{prefix}.log` | Human-readable regression report; some printed numbers are rounded. |
| `{prefix}.results.csv.checkpoint.json` | Completion and content/provenance record for `--restart`. |
| `LDSC_Batch_Status.csv` | Every requested batch's current outcome and attempt count. |

Native batch CSVs contain `p1,p2,rg,se,z,p,h2_int,h2_int_se,gcov_int,gcov_int_se`,
the available `h2_obs,h2_obs_se` or `h2_liab,h2_liab_se` pair, and
`gcov_native_obs,h2_p1_pair_obs,h2_p2_pair_obs`. Collection adds the scale label
and derived covariance/normalization columns. Use numerical CSVs for precise
values rather than copying rounded values from logs.

`LDSC_Batch_Status.csv` columns are:

| Column | Meaning |
| --- | --- |
| `Reference` | Trait used as the first trait in the batch. |
| `Batch` | Zero-based part index for that reference. |
| `Status` | `pending`, `completed`, `reused` or `execution_failed`. |
| `Result_File` | Numerical result path expected or produced. |
| `Error` | Execution/export failure explanation, when present. |
| `Restart_Reason` | Why a saved batch could or could not be reused. |
| `Attempts` | Number of attempts in this invocation; zero for a verified reused batch. |

`completed` means the command produced a structurally readable export, not
that every numerical estimate passed later QC. Exhausted command failures
stop before final collection regardless of the numerical failure policy.

Checkpoint JSON has `status`, `request` and `result_sha256`. The request records
reference and runtime fingerprints, run parameters, command and sumstats
content hashes. Preserve the result CSV and checkpoint together. Changed input
contents or settings can force a rerun even when filenames are unchanged.

## Worker attempts, runtime and errors

The following reports share the columns `Job,Attempt,Success,Error`:

| Report | Stage / when written |
| --- | --- |
| `LDSC_Extraction_Worker_Attempts.csv` | VCF extraction only. |
| `LDSC_Munging_Worker_Attempts.csv` | New munging only. |
| `LDSC_ChiSquare_Worker_Attempts.csv` | Either input route when the managed chi-square filter runs. |

`Job` identifies the trait job; `Attempt` counts from 1; `Success` records that
attempt's outcome; `Error` explains failure. A failed first attempt can be
followed by a successful second attempt. The first attempt uses the worker
pool, then failed jobs are retried sequentially. These reports use two total
attempts, independently of the pairwise `--ldsc_retries` setting.

**`LDSC_Runtime.json`** records the selected Conda executable/environment or
prefix, bcftools setting, backend, reference paths, restart/retry policy,
munging/result export settings and managed chi-square configuration. Its
`munging_export.source` distinguishes new native munging from existing files.
It is written early and is not proof that the run finished.

**`execution_errors.log`**, when present, contains failed subprocess commands
and captured stderr. A malformed result export can instead be explained in the
batch report without this log being created. Use the worker/batch status and
current compilation status together to locate the failed stage.

## Real-data verification

The final-table example in the [main guide](PYTHON_LDSC.md#a-real-data-result-row)
was checked on **2026-10-04**, using saved numerical exports from two real
munged-data analyses:

| Dataset | Traits | Directed result rows | Current collection compared with saved compilation |
| --- | --- | --- | --- |
| `T3_C0_C1_C2` | 43 | 1,849 | Identical values, columns, row order and file SHA-256. |
| `T3_C0_C1_NOISE` | 49 | 2,401 | Identical values, columns, row order and file SHA-256. |

All source batch files were unchanged by this check. Both compiled files had
the 21 columns described in the main guide. Their SHA-256 values were:

```text
T3_C0_C1_C2     85a2d8247cdc64fc362fc1e671191a1f9da38096e352b3e833971bf2943df5a0
T3_C0_C1_NOISE  6c2b3a9d2438d8f090d72eefa56ec83092b80270fd1884e69c3009cb920dc60a
```

The saved exports came from earlier direct refits of pinned
[CBIIT Python LDSC regression code](https://github.com/CBIIT/ldsc/tree/6c673952cee74bd5c57aef1555a03b1c015399a0)
on real munged inputs, using the unchanged regression classes and native result
table builder. They were not produced by a new full CLI invocation during this
documentation check. The current check reran **collection**, not extraction,
munging or regression. Numerical QC warnings remained visible in the reports;
matching results do not mean every scientific diagnostic passed without warning.

The original VCFs were not available for this check. VCF folder names, schemas
and input-mode behavior were verified against the current package source.
The original input data and saved validation artifacts are not bundled with the
repository. Older saved runs can have fewer audit columns or files than this
current output reference describes.
