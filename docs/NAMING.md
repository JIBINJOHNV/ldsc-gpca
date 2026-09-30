# Single interface specification — v0.6.0

Managed `prepare`, `ldsc`, `gpca`, `genomicsem ldsc` and `genomicsem gpca`
commands use one accepted name per option or manifest field. Multiword options
use underscores. Previous spellings and abbreviated options are rejected.
Update existing scripts and manifests before using this breaking release.

Raw `ldsc.py` and `munge_sumstats.py` passthrough commands retain the upstream
interface. Installer, Conda, bcftools and Docker command options are external to
this analysis interface. Variant-table headers and the original GWAMA code are
unchanged, including the pre-existing GWAMA variant-column handling.

## Arguments

| Argument | Meaning / availability |
| --- | --- |
| `--input` | Manifest CSV; all managed analyses |
| `--outdir` | Output directory; all managed analyses |
| `--n_cores` | Workers; all managed analyses, retaining each mode's defaults and limits |
| `--ld_ref` | LD-score reference directory; both LDSC backends |
| `--ld_weights` | Regression weights directory; both LDSC backends; default: `--ld_ref` |
| `--hm3` | HapMap reference for preparation or munging |
| `--munged_dir` | Existing munged directory; both LDSC backends |
| `--ldsc_results` | Pairwise text results for Python GPCA; genuine LDSCoutput RData for GenomicSEM GPCA |
| `--gpca_input_folder` | Prepared per-trait GWAMA inputs |
| `--splitby_chr` | `split` or `nosplit` |
| `--validate_only` | GPCA QC/PCA without GWAMA |
| `--dataset_id` | GWAMA export filename prefix |
| `--gwama_output_n_eff` | Exported effective sample-size override |
| `--gwama_output_info` | Exported INFO override |
| `--gzip_level` | Final-export gzip compression, 1–9; default 1 from v0.6.1 |
| `--prepare_workers` | Separate preparation-stage workers inside GPCA only |

Standalone `prepare` uses `--n_cores`; GPCA's `--n_cores` controls GWAMA,
concurrent chromosome reads, the default shared Polars thread pool (from v0.6.2),
and final-export pigz compression workers when pigz is installed;
`--prepare_workers` controls its earlier optional VCF preparation stage. These are
separate settings. Other managed options use the same underscore spelling; see
`--help`, `--prepare_help`, or `--postprocess_help` for the complete list.

Python `ldsc --munged_dir DIR` needs `--ldsc_only` to skip extraction/munging.
For native GenomicSEM LDSC, `--munged_dir DIR` selects existing munged files.
File types and statistical conventions are unchanged by their argument names.
Only the explicit `--exclude_mhc` switch enables MHC removal; omission keeps MHC.

Python `--ld_weights` goes to native LDSC's `--w-ld-chr`, and `--ld_ref` goes to
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
| `munged_file` | Already munged file paths; native GenomicSEM file-path mode |
| `N` | Optional constant sample size for native GenomicSEM munging |

| Analysis/mode | Required headers |
| --- | --- |
| Preparation or automatic GPCA preparation | `traitname,vcf_files` |
| Python LDSC from VCFs | `traitname,vcf_files,ref,population_prevalence,sample_prevalence` |
| Python LDSC with `--ldsc_only` | `traitname,ref,population_prevalence,sample_prevalence` |
| Native GenomicSEM LDSC from unmunged files | `traitname,sumstats_file,sample_prevalence,population_prevalence` |
| Native GenomicSEM LDSC with `--munged_dir` | `traitname,sample_prevalence,population_prevalence` |
| Native GenomicSEM LDSC with `--munged_input` | `traitname,munged_file,sample_prevalence,population_prevalence` |
| Either GPCA with prepared inputs or `--validate_only` | `traitname` |

A master CSV can contain all the listed fields. Modes require only their relevant
fields; other metadata can be present. Removed manifest names are rejected even
if their new counterparts are also present. Headers are case-sensitive and unique.
Trait names stay strings, including leading zeros; row order is preserved.
Existing backend/audit report names are retained. This specification applies to
CLI options and input manifests, not variant tables or output report schemas.

The biological rules are unchanged: Python uses NEF when population prevalence
is absent, and NC+NCO when supplied; native GenomicSEM requires both prevalences
or neither. Python LDSC reads `{traitname}.sumstats.gz` from the munged directory;
it does not use the `munged_file` column. Use absolute VCF paths for Python LDSC;
preparation/native GenomicSEM resolve relative paths beside the manifest.

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
