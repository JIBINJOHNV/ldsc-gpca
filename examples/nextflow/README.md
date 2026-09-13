# Nextflow: local Docker and Google Cloud Batch

These examples use the same image for local and cloud tasks:

- `smoke.nf`: checks the installed tools; no scientific input required.
- `gpca_qc.nf`: runs Python-LDSC QC and PCA only, with strict package defaults.
  It does not run LDSC estimation, GWAMA or automatic export.

Nextflow and Java run on your launcher (workstation or VM), not inside the analysis
image. Install Nextflow using its [official guide](https://docs.seqera.io/nextflow/install).
For these examples, use Nextflow 26.04.6 and Java 17–26. Pin your launcher version
and image digest for reproducibility. Java 24 was available for local testing.

## 1. Local container test

From the repository root, build the image as described in the main README, then:

```bash
nextflow run examples/nextflow/smoke.nf -profile local \
  --image ldsc-gpca:0.5.0 --outdir nf-smoke
```

Success produces `nf-smoke/tools.txt`. The local profile requires a running Docker
daemon and explicitly requests linux/amd64. The examples default to container UID
0 so Batch staging directories are writable. Locally, use your UID/GID if desired:

```bash
nextflow run examples/nextflow/smoke.nf -profile local \
  --container_options "--user $(id -u):$(id -g)" --outdir nf-smoke-user
```

This is a Linux filesystem identity, not Google IAM. Avoid running untrusted images.

## 2. Run QC/PCA with staged files

```bash
nextflow run examples/nextflow/gpca_qc.nf -profile local \
  --manifest /absolute/path/selected_traits.csv \
  --ldsc /absolute/path/all_pairwise_python_ldsc.csv.gz \
  --outdir nf-qc
```

The manifest is comma-separated with unique, non-empty `traitname` values. The
LDSC file must meet the package's complete pair/self-pair and numeric QC contracts;
see the main README. Nextflow stages both files using `path` inputs. Outputs are
published under `nf-qc/qc/`. Missing input arguments stop before task submission.

Resource defaults in `nextflow.config`: 2 CPUs, 8 GB memory, 50 GB disk and 2 hours.
These are starting values, not sizing guarantees. Override `--task_cpus`,
`--task_memory`, `--task_disk` and `--task_time` for your dataset.

## 3. Publish your image for GCP

The local tag is not accessible to cloud workers. Authenticate your Google Cloud
CLI and choose a project with billing enabled. The following commands are examples
for you to run; no cloud resources or image publication are performed automatically:

```bash
GCP_PROJECT=your-project-id
GCP_REGION=us-central1
GCP_REPOSITORY=analysis-images
gcloud artifacts repositories create "$GCP_REPOSITORY" \
  --repository-format=docker --location="$GCP_REGION" --project="$GCP_PROJECT"
gcloud auth configure-docker "$GCP_REGION-docker.pkg.dev"
IMAGE="$GCP_REGION-docker.pkg.dev/$GCP_PROJECT/$GCP_REPOSITORY/ldsc-gpca:0.5.0"
docker tag ldsc-gpca:0.5.0 "$IMAGE"
docker push "$IMAGE"
```

Create the repository only if it does not already exist. Use the registry digest
reported after pushing (`.../ldsc-gpca@sha256:...`) as `--image` for reproducible runs.
See [Artifact Registry push/pull instructions](https://cloud.google.com/artifact-registry/docs/docker/pushing-and-pulling).

## 4. Google Batch prerequisites and smoke test

Before submitting jobs, configure Batch/Compute/Storage/Artifact Registry APIs,
quotas, a GCS bucket and a Batch task service account. The launcher must have
permission to submit jobs and use that account; workers need image-pull and
input/output bucket permissions. Follow the
[Nextflow Google Batch setup](https://docs.seqera.io/nextflow/google) and your
organization's IAM policy. Do not bake service-account keys into the image.

On a workstation, configure Application Default Credentials:

```bash
gcloud auth application-default login
```

On a GCP launcher VM, prefer its attached service account. Then run:

```bash
nextflow run examples/nextflow/smoke.nf -profile gcp \
  --project YOUR_PROJECT_ID --region us-central1 \
  --service_account YOUR_BATCH_SERVICE_ACCOUNT_EMAIL \
  --image REGION-docker.pkg.dev/PROJECT/REPOSITORY/ldsc-gpca@sha256:YOUR_DIGEST \
  -work-dir gs://YOUR_BUCKET/ldsc-gpca/work \
  --outdir gs://YOUR_BUCKET/ldsc-gpca/smoke
```

Replace every placeholder. The GCP profile selects the x86-64 N2 machine family
(`--machine_type 'n2-*'`). Override only with another compatible x86-64 family;
do not select ARM for this linux/amd64 image. Ensure enough disk for image and data.
The GCS work URI must include a subdirectory, not just a bucket name. Nextflow's
Google plugin and cloud data-transfer helpers require the relevant network access.
The analysis image itself does not need a Docker daemon or embedded credentials.

After smoke succeeds, replace `smoke.nf` with `gpca_qc.nf` and add
`--manifest gs://YOUR_BUCKET/inputs/selected_traits.csv` and
`--ldsc gs://YOUR_BUCKET/inputs/all_pairwise_python_ldsc.csv.gz`.
Set a separate output prefix for each analysis. Use `-resume` to reuse successful
tasks; it does not correct invalid scientific inputs.

For your own GWAMA/preparation processes, declare every file or directory as a
Nextflow `path` input and pass the staged paths to `ldsc-gpca`. Merely embedding
`gs://` strings in a manifest does not stage the referenced VCFs. Limit package
parallelism (e.g. GWAMA `--cores`) to `task.cpus`. Do not call `conda activate`
or nest `docker run` inside a process script.

## Validation boundary

Locally verified with Nextflow 26.04.6, Java 24 and the linux/amd64 image:

- `smoke.nf` completed successfully, including resource tracing and LDSC/R checks.
- `gpca_qc.nf` completed on a synthetic three-trait input. Matrix ordering, CTI,
  eigenvalues and PC1 loadings matched independent NumPy calculations.
- Missing parameters stopped before submission; an absent trait produced a failed
  task (exit 1), not a successfully published analysis.
- The GCP profile parsed with the Nextflow-selected `nf-google` plugin 1.27.3.

Local image/Nextflow checks are distinct from GCP deployment validation. No live
GCP job has been submitted, no registry push has been performed, and project IAM,
networking and quotas have not been tested here. The cloud examples require a
successful smoke run in your own project before a production analysis.
