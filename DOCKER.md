# Docker installation

From the repository directory, with Docker Desktop or Docker Engine running:

```bash
docker build --platform linux/amd64 --progress=plain -t ldsc-gpca:0.5.0 .
docker run --rm --platform linux/amd64 ldsc-gpca:0.5.0
docker run --rm --platform linux/amd64 ldsc-gpca:0.5.0 genomicsem ldsc --help
docker run --rm --platform linux/amd64 ldsc-gpca:0.5.0 gpca --help
```

The default command displays help. Arguments after the image name are passed
to `ldsc-gpca`. The build reuses `scripts/setup_environments.sh`, installing
bcftools, the Python/R package environment, pinned GenomicSEM, and the isolated
CBIIT LDSC environment. Its existing installation checks must all pass.
Docker is not required inside the container.

Use Linux amd64 initially; native arm64 dependency resolution has not been
validated. On Apple Silicon this requests emulation and may be substantially
slower. Allow sufficient disk space and memory for two Conda environments and
compilation; installation requires internet access.

## Data and paths

The image includes the reviewed modified GWAMA v1.2.6 script as the default for
both GPCA backends. GWAS data and LD reference files must be mounted explicitly.
Mount a custom R script only if overriding `--source_path`. For example:

```bash
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v /absolute/path/to/project:/work -w /work \
  ldsc-gpca:0.5.0 gpca \
  --input /work/selected_traits.csv \
  --python_ldsc /work/all_pairwise_python_ldsc.csv \
  --gpca_input_folder /work/gpca_inputs \
  --outdir /work/results --splitby_chr split
```

Replace the host directory and filenames with your own. Paths inside manifests
must also resolve inside the container. The mounted output directory must be
writable by the selected user. Supply scientifically justified export overrides
only when needed; this example does not invent INFO or sample sizes. The bundled
GWAMA source does not output INFO, so automatic export requires a scientifically
justified `--gwama-output-info` override or will stop after GWAMA completes.

## Reproducibility and verification

The base image tag and environment specifications are not complete dependency
locks. GenomicSEM and CBIIT LDSC source commits are pinned by the existing
installer/specification. Actual Conda package URLs, pip versions and R session
information are recorded inside `/opt/environments/`. Archive the built image
digest and these records for reproducible analyses. A successful build checks
software availability, not scientific correctness of an analysis dataset.

## Verified local build

The Linux amd64 image `ldsc-gpca:0.5.0` was successfully built and run locally.
The image must be rebuilt after package changes, including bundled-source updates.
The bundled-source rebuild succeeded with image ID
`sha256:a2cfb478bffa304f31095aa6a9f5e13ce111d205190a06a45de983f32e7b0b0d`.
The installed wheel's GWAMA SHA-256 matches the maintainer-supplied file exactly.
Both installed GPCA backends passed default/override selection, missing-source
failure and signed-loading numerical tests. The updated source also passed
102 Python regression tests and 12 R-module checks.
All installer checks passed. Container tests passed for the CLI help commands,
non-root access to the isolated LDSC environment, GenomicSEM loading and UTC
timezone detection. All 27 preparation tests passed, including synthetic VCF
extraction, output columns, boundary conditions and invalid/missing inputs.
GenomicSEM emits upstream namespace-import warnings but loads successfully.

A full real-data LDSC-to-GWAMA analysis and native arm64 builds were not tested.
The image is local only; it has not been published to a registry.

See the [official Dockerfile documentation](https://docs.docker.com/build/concepts/dockerfile/).
