# ldsc-gpca

Python LDSC and genomicPCA/GWAMA through a unified command-line interface.

* `ldsc-gpca ldsc`: VCF extraction, munging, and pairwise Python LDSC.
* `ldsc-gpca gpca`: genomicPCA/GWAMA from Python LDSC tables using the bundled R script.

This package does not replace LDSC, convert the R method to Python, or fabricate
GenomicSEM sampling covariance matrices. The R implementation is bundled unchanged.
This initial package version preserves the existing analysis defaults and P-value
conversion. Packaging tests are not a substitute for scientific validation on real data.

## Installation

Use Python 3.10–3.12 in an isolated environment (local validation used Python 3.10).
Clone the repository, then install the package:

```bash
git clone https://github.com/JIBINJOHNV/ldsc-gpca.git
cd ldsc-gpca
python -m venv .venv
source .venv/bin/activate
python -m pip install .
ldsc-gpca --help
ldsc-gpca ldsc --help
ldsc-gpca gpca --help
```

For development, use `python -m pip install -e .` instead.
You can also run `python -m ldsc_gpca ...`.
Use the installed commands above for both workflows. The R script is packaged
under `src/ldsc_gpca/r/` and located automatically by `ldsc-gpca gpca`.

### External dependencies

* LDSC: Docker on PATH, a running Docker daemon, and `jibinjv/ldsc:v3`.
  Docker must be able to mount all input/output/reference folders. Use absolute paths.
  The current extraction code assumes a POSIX host (Linux/macOS); native Windows
  execution is not supported. No image is pulled by package installation.
* GPCA: `Rscript` on PATH, R packages `argparse`, `data.table`, and `glue`
  (`parallel` is supplied with R). Install them in R with
  `install.packages(c("argparse", "data.table", "glue"))`.
* The modified Fürtjes GWAMA source must be supplied via `--source_path` along with
  any dependencies it requires. It is not included or altered by this package.
* Supply your own LD-score reference files, SNP list, and GWAS data.

Python dependency ranges are in `pyproject.toml`. For reproducible analyses, record
`pip freeze`, `sessionInfo()` from R, and the Docker image digest; the image tag is
not an immutable scientific environment.

## Python LDSC

```bash
ldsc-gpca ldsc \
  --input_file /absolute/path/traits.csv \
  --output_folder /absolute/path/ldsc_output \
  --ld_ref_snp_file /absolute/path/w_hm3.snplist \
  --ld_ref /absolute/path/eur_w_ld_chr \
  --n_cores 5 --ldsc-retries 1
```

Required manifest headers: `gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence`.
Names must be unique/non-empty. Blank prevalence values are allowed, but the headers
are required. `ref=yes` requests that trait against every listed trait, including
itself. For a complete matrix for GPCA, set `ref=yes` for **every selected trait**.

VCFs must contain the fields consumed by the original extraction filters:
INFO/AF, INFO/EUR, FORMAT/SI, FORMAT/AF, FORMAT/EZ, FORMAT/LP and FORMAT/NEF.
Traits with population prevalence additionally require FORMAT/NC and FORMAT/NCO.
Missing population prevalence keeps NEF and does not use sample prevalence.
Supplied population prevalence selects NC+NCO and fills missing sample prevalence
with the median per-variant case fraction. If all population prevalences are blank,
no prevalence flags are passed to LDSC.

MHC exclusion defaults to off. Extraction MAF defaults to 0.01 and munging MAF
to 0.005. Standard LDSC removes all palindromic SNPs in munging regardless of the
optional extraction-stage palindromic filter. See `ldsc-gpca ldsc --help`.

### Failure handling and reuse

`--ldsc-retries 1` means one initial attempt plus one retry, with identical inputs
and settings. Successful batches are not repeated during this run. Exhausted retries
stop the workflow before compilation; other submitted jobs finish before exit.
Validation failures are not retried. Retries do not repair files or statistical issues.
Execution failures are appended to `execution_errors.log` with attempt labels.

`--ldsc_only` reuses munged files and reruns **all** requested LDSC batches, not just
failed batches. Filters are not reapplied; recorded filter mismatches fail and unknown
legacy filter settings cause a warning. No selective cross-run batch resume is provided.
Extraction failures stop the run; handled munging failures remove affected traits.
Missing/malformed results or non-finite rg fail compilation. Existing result files
from earlier runs are not deleted on failure and must not be mistaken for new output.

Outputs include `ldsc_results.csv`, `LDSC_Trait_Prevalence_Metadata.csv`, per-trait
`.prevalence.json` sidecars, munged files, and LDSC logs. The compiled table reports
target-trait heritabilities and their scale. Quantitative and binary NEF-unconverted
estimates should not be assumed to be on a common liability scale.

## genomicPCA/GWAMA

```bash
ldsc-gpca gpca \
  --input /absolute/path/selected_traits.csv \
  --python_ldsc /absolute/path/ldsc_output/ldsc_results.csv \
  --gpca_input_folder /absolute/path/gpca_inputs \
  --source_path /absolute/path/my_GWAMA_26032020.R \
  --outdir /absolute/path/gpca_output \
  --splitby_chr split
```

The GPCA manifest is a separate CSV requiring `traitname`; its row order defines
matrix and GWAMA order. The LDSC file must supply every selected self-pair and pair,
with rg/diagnostic/intercept columns and supported observed or liability heritability
columns. Extra traits are allowed. The bundled R help describes QC, scale selection,
correlation/covariance modes, input file naming, and audit outputs in detail.

The externally supplied GWAMA function must already use the Fürtjes modification:
`W <- t(t(sqrt(N)) * h2)` and `sqrt_W <- W`, where its `h2` argument receives PC
loadings. The package does not modify this function. Python pairwise LDSC output
does not supply GenomicSEM's full `V`/`V_Stand`, so it cannot support paLDSC here.

### Automatic GWAMA post-processing

After successful GWAMA, results are combined automatically. No `--postprocess`
flag is needed. The selected-column summary defaults to `--outdir/harmonised/`.
Optionally append these settings to the GPCA command above:

```bash
--harmonised-output /absolute/path/harmonised \
--dataset-id cluster1 \
--gwama-output-n-eff 330000 \
--gwama-output-info 0.9
```

Post-processing always runs after successful GWAMA. It requires successful R execution
and a successful, updated `GWAMA_Run_Status.csv`. Validation-only runs skip it.
Only current-run output prefixes in that status file are read; other matching files
in the folder are ignored. The expected source filename suffix is
`.N_weighted_GWAMA.results.txt.gz`; custom GWAMA functions must use that convention.

| Output | Location |
| --- | --- |
| `{name}_GWAMA_combined_results.txt.gz` | `--outdir` |
| `{name}_GPCA_inputs.txt.gz` | `--harmonised-output` (default: `--outdir/harmonised/`) |
| `{name}_postprocess.json` (sources, row count, overrides) | `--outdir` |

`--dataset-id` sets `{name}` and defaults to the name of the `--outdir` folder. Both compressed
tables are tab-delimited. Results are sorted numerically by chromosome and position;
the combined table includes `count_question`, `count_plus`, and `count_minus` from
`Direction`. Missing required columns, empty results, invalid direction strings or
positions, duplicate SNPIDs, missing files and unchanged old results cause an error.

The selected-column summary contains:

```text
SNPID CHR BP EA OA EAF N_eff BETA SE Z PVAL INFO
```

`--gwama-output-n-eff` and `--gwama-output-info` are optional, finite overrides,
**not automatic defaults**. They affect only the exported GWAMA selected-column
summary, never LDSC/GPCA calculations or variant filtering.
N_eff must be positive and INFO must lie in [0,1]. Without an override, that column
must exist and its reported values are preserved. An override can also supply an
absent column. Overrides affect only the selected-column summary, not the combined
results. They are user-supplied constants, not values estimated by this operation.
Despite its filename, this summary is not the nine-column per-trait GWAMA input
format and does not perform allele harmonisation or genome-build conversion.

Original files are kept by default. Add `--archive-chromosomes` to move the current
run's source results and matching logs into `--outdir/chromosome_wise/` only after
both tables and the audit have been saved. Existing output/archive paths are never
overwritten. An archive failure leaves saved outputs in place and is reported as an
error; some originals may already have moved. Choose a new name/output directory
for repeat exports. No automatic post-processing retry is performed.

Show these options without requiring R:

```bash
ldsc-gpca gpca --postprocess-help
```

## Package layout

```text
ldsc-gpca/
├── README.md
├── pyproject.toml
├── MANIFEST.in
└── src/ldsc_gpca/
    ├── cli.py              # Command dispatch
    ├── ldsc_cli.py         # LDSC options and orchestration
    ├── extraction.py       # VCF filtering and prevalence preparation
    ├── munging.py          # Summary-statistic munging
    ├── pairwise.py         # LDSC batches and retries
    ├── results.py          # Result validation and compilation
    ├── utils.py            # Shared helpers
    ├── gpca.py             # R launcher
    ├── postprocess.py      # Automatic GWAMA combination and export
    └── r/gpsca_gwama_python_ldsc.r
```

## Attribution

Method: [Anna Fürtjes genomicPCA tutorial](https://annafurtjes.github.io/genomicPCA/25082021_geneticPCA_explanation.html).
Upstream: [Bulik-Sullivan LDSC](https://github.com/bulik/ldsc).
Packaging follows [setuptools entry points](https://setuptools.pypa.io/en/latest/userguide/entry_point.html)
and [package data](https://setuptools.pypa.io/en/latest/userguide/datafiles.html).

## License

No redistribution license has been selected for this repository. Public availability
does not grant an open-source license. Third-party software retains its own licensing.
