# Command and manifest naming

<a id="single-interface-specification--v060"></a>

[Choose a workflow](../README.md#choose-a-workflow) · [Input formats](INPUTS.md) · [Full reference](REFERENCE.md)

The canonical names below apply to v0.8.0. The naming changes began in v0.6.0;
see [pipeline migration](PIPELINE.md#migrate-existing-pipeline-commands) for the
additional VCF-only input restriction introduced in v0.8.0.

Managed `pipeline`, `prepare`, `ldsc`, `gpca`, `genomicsem ldsc` and `genomicsem gpca`
commands use one accepted name per option or manifest field. Multiword options
use underscores. Previous argument spellings and abbreviated options are rejected.
Extra manifest columns are allowed, including old alias names as annotations;
they do not substitute for required canonical headers.
Use these names when updating older scripts and manifests.

Raw `ldsc.py` and `munge_sumstats.py` passthrough commands retain the upstream
interface. Installer, Conda, bcftools and Docker command options are external to
this analysis interface. Variant-table headers and the original GWAMA code are
unchanged, while the GWAMA reader now accepts extra/reordered input columns and selects the
required canonical fields before calling the unchanged GWAMA function.

## Arguments

| Argument | Meaning / availability |
| --- | --- |
| `--input` | Manifest CSV; all managed analyses |
| `--outdir` | Output directory; all managed analyses |
| `--n_cores` | Workers; all managed analyses, retaining each mode's defaults and limits |
| `--ldsc_backend` | `python` (default) or `genomicsem`; `ldsc` and `pipeline`, selecting both munging and regression |
| `--munge_backend` | `python` (default) or `genomicsem`; standalone `prepare --mode ldsc`/`both` when munging runs |
| `--ld_ref` | LD-score reference directory; both LDSC backends |
| `--ld_weights` | Regression weights directory; both LDSC backends; default: `--ld_ref` |
| `--hm3` | HapMap reference for preparation or munging |
| `--munged_dir` | Existing munged directory; both LDSC backends |
| `--ldsc_results` | Pairwise text results for Python GPCA; genuine LDSCoutput RData for GenomicSEM GPCA |
| `--gpca_input_folder` | Prepared per-trait GWAMA inputs; standalone GPCA commands only |
| `--splitby_chr` | `split` or `nosplit` |
| `--validate_only` | GPCA: QC/PCA from saved estimates; pipeline: VCF preparation, LDSC and PCA. Both skip GWAMA/export. |
| `--dataset_id` | GWAMA export filename prefix |
| `--gwama_output_n_eff` | Exported effective sample-size override |
| `--gwama_output_info` | Exported INFO override |
| `--gzip_level` | Final-export gzip compression, 1–9; default 1 from v0.6.1 |
| `--prepare_workers` | Separate preparation-stage workers inside GPCA or pipeline |

Standalone `prepare` uses `--n_cores`; GPCA's `--n_cores` controls GWAMA,
concurrent chromosome reads, the default shared Polars thread pool (from v0.6.2),
and final-export pigz compression workers when pigz is installed;
`--prepare_workers` controls its earlier optional VCF preparation stage. These are
separate settings. Other managed options use the same underscore spelling; see
`--help`, `--prepare_help`, or `--postprocess_help` for the complete list.

Python `ldsc --munged_dir DIR` requires `--ldsc_only`; supplying the directory
alone fails before output creation. Fresh runs use `<outdir>/ldsc_input`.
For GenomicSEM LDSC, `--munged_dir DIR` selects existing munged files.
File types and statistical conventions are unchanged by their argument names.
`ldsc-gpca genomicsem ldsc` remains available as the GenomicSEM entry point;
it implies that backend and does not take `--ldsc_backend`.
Only the explicit `--exclude_mhc` switch enables MHC removal; omission keeps MHC.

Python `--ld_weights` goes to Python LDSC's `--w-ld-chr`, and `--ld_ref` goes to
`--ref-ld-chr`. Their resolved paths are recorded in `LDSC_Runtime.json`.
The upstream definitions are in [pinned CBIIT LDSC](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldsc.py).

## Manifest headers

| Header | Meaning |
| --- | --- |
| `traitname` | Unique, non-empty trait identifier; row order is authoritative |
| `vcf_files` | VCF paths |
| `ref` | Python LDSC reference selection; `yes` for every trait for complete GPCA coverage |
| `sample_prevalence` | Case fraction in the study, when applicable |
| `population_prevalence` | Population disease prevalence, when applicable |
| `sumstats_file` | Unmunged summary-statistic file paths |
| `munged_file` | Already munged file paths; GenomicSEM file-path mode |
| `N` | Optional constant sample size for GenomicSEM munging |

| Analysis/mode | Required headers |
| --- | --- |
| Preparation or automatic GPCA preparation | `traitname,vcf_files` |
| Python LDSC from VCFs | `traitname,vcf_files,ref,population_prevalence,sample_prevalence` |
| Python LDSC with `--ldsc_only` | `traitname,ref,population_prevalence,sample_prevalence` |
| GenomicSEM LDSC from unmunged files | `traitname,sumstats_file,sample_prevalence,population_prevalence` |
| GenomicSEM LDSC with `--munged_dir` | `traitname,sample_prevalence,population_prevalence` |
| GenomicSEM LDSC with `--munged_input` | `traitname,munged_file,sample_prevalence,population_prevalence` |
| Either GPCA with prepared inputs or `--validate_only` | `traitname` |

A master CSV can contain all the listed fields. Modes require only their relevant
fields; other metadata can be present. Old manifest names may be extra metadata
but do not substitute for required canonical fields. Headers are case-sensitive and unique.
Trait names stay strings, including leading zeros; row order is preserved.
Existing backend/audit report names are retained. This specification applies to
CLI options and input manifests, not variant tables or output report schemas.

The biological rules are unchanged: Python uses NEF when population prevalence
is absent, and NC+NCO when supplied; GenomicSEM requires both prevalences
or neither. Python LDSC reads `{traitname}.sumstats.gz` from the munged directory;
it does not use the `munged_file` column. Use absolute VCF paths for Python LDSC;
preparation/GenomicSEM resolve relative paths beside the manifest.

A tab-separated HapMap table with `SNP,A1,A2` works with all `--hm3` options;
preparation uses only SNP identifiers. See [REFERENCE.md](REFERENCE.md) for
unchanged variant-table formats and statistical requirements.

## Your GPCA command

```bash
ldsc-gpca gpca \
  --input "${input_csv}" \
  --ldsc_results "${python_ldsc_csv}" \
  --gpca_input_folder "${gpca_inputs}/gpca_inputs" \
  --outdir "${out_folder}" \
  --splitby_chr split \
  --n_cores 22 \
  --gwama_output_n_eff 33000 \
  --gwama_output_info 0.9 \
  --failed_ldsc_action drop_traits \
  --dataset_id "${cluster_id}"
```

This interface update does not itself accelerate analysis. The earlier single-
thread numerical-library defaults remain enabled unless the environment overrides
them. Install the updated package in the environment that runs your commands.

<a id="one-command-pipeline-v070"></a>

## One-command pipeline with VCF inputs

`ldsc-gpca pipeline --ldsc_backend python|genomicsem` reuses the managed
interface and writes canonical per-stage manifests. Python `--chisq_max` accepts
a positive integer or `auto`; GenomicSEM accepts a positive finite number.
Omitting it preserves the selected backend's default. Python
`--rg_normalization pair|trait_wide` selects the downstream correlation source.
`--n_cores` stays positive for pipeline (Python default 5, GenomicSEM default 1), and
`--prepare_workers` remains a separate setting. There are no compatibility aliases.
The pipeline accepts only VCF study inputs. It generates both LDSC results and
GWAMA tables internally; users do not supply `--ldsc_results` or `--gpca_input_folder`.
Use standalone commands for raw/munged inputs or completed estimates.
GenomicSEM VCF conversion is automatic; omit `--vcf_input` from pipeline.
It requires a fresh output directory; `--restart` remains a standalone LDSC option.
See the [pipeline user guide](PIPELINE.md).
