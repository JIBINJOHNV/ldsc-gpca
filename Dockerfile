# syntax=docker/dockerfile:1
# Build with --platform linux/amd64 for the older pinned LDSC dependencies.
FROM ubuntu:24.04
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl bzip2 build-essential gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    && apt-get clean
WORKDIR /opt/ldsc-gpca
COPY environment.yml environment.ldsc.yml pyproject.toml README.md ./
COPY scripts/setup_environments.sh scripts/install_genomicsem.R ./scripts/
COPY src/ ./src/
# Includes import, R-package, executable and CLI checks; fails on any error.
RUN bash scripts/setup_environments.sh --no-activate /opt/environments
ENV PATH="/opt/environments/main/bin:/opt/environments/conda/bin:${PATH}" \
    LDSC_GPCA_LDSC_PREFIX="/opt/environments/ldsc" \
    CONDA_EXE="/opt/environments/conda/bin/conda" \
    TZ="Etc/UTC"
RUN apt-get update && apt-get install -y --no-install-recommends procps && apt-get clean \
    && mkdir /work
# Nextflow supplies the task command; tools are available without activation.
USER root
WORKDIR /work
ENTRYPOINT []
CMD ["ldsc-gpca", "--help"]
