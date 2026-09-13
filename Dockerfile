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
COPY --chmod=755 scripts/docker-entrypoint.sh /usr/local/bin/ldsc-gpca-entrypoint
ENV PATH="/opt/environments/main/bin:/opt/environments/conda/bin:${PATH}" \
    LDSC_GPCA_LDSC_PREFIX="/opt/environments/ldsc" \
    CONDA_EXE="/opt/environments/conda/bin/conda" \
    TZ="Etc/UTC"
RUN apt-get update && apt-get install -y --no-install-recommends procps && apt-get clean \
    && useradd --create-home --uid 10001 analysis && mkdir /work && chown analysis /work
USER analysis
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/ldsc-gpca-entrypoint"]
CMD ["ldsc-gpca", "--help"]
