#!/usr/bin/env bash
# Sourcing activates Bash/zsh; installation runs separately to protect the shell.
if [ -n "${ZSH_VERSION:-}" ]; then
  eval '_ldsc_setup_file=${(%):-%x}'
  _ldsc_setup_sourced=1
elif [ "${BASH_SOURCE[0]}" != "$0" ]; then
  _ldsc_setup_file=${BASH_SOURCE[0]}
  _ldsc_setup_sourced=1
else
  _ldsc_setup_sourced=0
fi
if [ "$_ldsc_setup_sourced" = 1 ]; then
  bash "$_ldsc_setup_file" --no-activate "$@" || return $?
  case " ${*} " in *' --help '*|*' -h '*|*' --no-activate '*) return 0 ;; esac
  _ldsc_setup_root=${1:-"$(cd "$(dirname "$_ldsc_setup_file")/.." && pwd)/.environments/named"}
  . "$_ldsc_setup_root/activate.sh"
  return $?
fi
set -euo pipefail
activate=1
named=1
root=''
for argument in "$@"; do
  case "$argument" in
    --help|-h)
      echo 'Usage: bash scripts/setup_environments.sh --no-activate'
      echo 'Then: conda activate ldsc-gpca (initialize Conda in your shell first).'
      echo 'Or: bash scripts/setup_environments.sh [--no-activate] [ENV_ROOT]'
      echo 'Default: named environments ldsc-gpca and ldsc-gpca-ldsc; records in <repository>/.environments/named.'
      echo 'Advanced: an explicit ENV_ROOT keeps legacy prefix mode (ENV_ROOT/main and ENV_ROOT/ldsc).'
      echo 'Installs main + isolated LDSC, checks tools, and activates after success.'
      echo 'Linux x86_64/aarch64 and macOS Intel/Apple Silicon; subject to dependency availability.'
      echo 'Bootstraps checksum-verified Miniforge if Conda is absent. No sudo or shell-profile edits.'
      echo 'Sourcing activates this Bash/zsh terminal; bash execution opens an activated Bash shell on a terminal.'
      echo '--no-activate: install/check only (also used automatically without a terminal).'
      exit 0 ;;
    --no-activate) activate=0 ;;
    -*) echo "Unknown option: $argument" >&2; exit 2 ;;
    *) if [[ -n "$root" ]]; then echo 'Expected at most one ENV_ROOT argument.' >&2; exit 2; fi
       root=$argument; named=0 ;;
  esac
done
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=${root:-"$repo/.environments/named"}
os=$(uname -s)
cpu=$(uname -m)
case "$os:$cpu" in
  Darwin:arm64|Darwin:x86_64) installer_os=MacOSX ;;
  Linux:x86_64|Linux:aarch64) installer_os=Linux ;;
  *) echo "Unsupported platform: $os/$cpu. Use a supported 64-bit Linux or macOS host." >&2; exit 2 ;;
esac
for required in environment.yml environment.ldsc.yml scripts/install_genomicsem.R pyproject.toml; do
  [[ -f "$repo/$required" ]] || { echo "Missing repository file: $required" >&2; exit 2; }
done
mkdir -p "$root"
root=$(cd "$root" && pwd)
main_prefix="$root/main"
ldsc_prefix="$root/ldsc"
if [[ -e "$main_prefix" || -e "$ldsc_prefix" || -e "$root/activate.sh" ]]; then
  echo 'Environment target already exists. Choose a fresh ENV_ROOT; no environments were changed.' >&2
  exit 2
fi
trap 'status=$?; echo "Installation/check failed (exit $status) on $os/$cpu. Partial files remain in $root. No activation performed. Dependency pins were not relaxed." >&2; exit "$status"' ERR
conda_bin=${CONDA_EXE:-conda}
if ! command -v "$conda_bin" >/dev/null; then
  if [[ -n "${CONDA_EXE:-}" ]]; then echo "CONDA_EXE is invalid: $CONDA_EXE" >&2; exit 2; fi
  version=26.7.2-0
  asset="Miniforge3-${version}-${installer_os}-${cpu}.sh"
  url="https://github.com/conda-forge/miniforge/releases/download/$version/$asset"
  command -v curl >/dev/null || { echo 'Install curl first to bootstrap Miniforge.' >&2; exit 2; }
  downloads=$(mktemp -d "$root/download.XXXXXX")
  curl --fail --location --retry 3 "$url" --output "$downloads/$asset"
  curl --fail --location --retry 3 "$url.sha256" --output "$downloads/$asset.sha256"
  expected=$(awk 'NR==1 {print $1}' "$downloads/$asset.sha256")
  if command -v sha256sum >/dev/null; then
    actual=$(sha256sum "$downloads/$asset" | awk '{print $1}')
  else
    actual=$(shasum -a 256 "$downloads/$asset" | awk '{print $1}')
  fi
  [[ ${#expected} == 64 && "$actual" == "$expected" ]] || { echo 'Miniforge checksum mismatch; refusing to execute.' >&2; exit 2; }
  bash "$downloads/$asset" -b -p "$root/conda"
  conda_bin="$root/conda/bin/conda"
fi
conda_bin=$(command -v "$conda_bin")
base=$("$conda_bin" info --base)
[[ -f "$base/etc/profile.d/conda.sh" ]] || { echo 'Conda shell activation script is missing.' >&2; exit 2; }
manager=$conda_bin
if [[ -x "$base/bin/mamba" ]]; then manager="$base/bin/mamba"; fi
echo "Installing on $os/$cpu using $manager. Older pinned LDSC dependencies must resolve for this platform."
if [[ $named == 1 ]]; then
  # Refuse either existing name before creating anything; preserve user environments.
  existing=$("$conda_bin" env list --json)
  if ! printf '%s' "$existing" | "$base/bin/python" -c 'import json,os,sys; names={os.path.basename(p) for p in json.load(sys.stdin)["envs"]}; sys.exit(bool(names & {"ldsc-gpca","ldsc-gpca-ldsc"}))'; then
    echo 'Could not verify fresh environment names, or ldsc-gpca/ldsc-gpca-ldsc already exists. Inspect conda env list; nothing was overwritten.' >&2
    exit 2
  fi
  "$manager" env create --name ldsc-gpca --file "$repo/environment.yml" --yes
  "$manager" env create --name ldsc-gpca-ldsc --file "$repo/environment.ldsc.yml" --yes
  main_prefix=$("$conda_bin" run --name ldsc-gpca python -c 'import sys; print(sys.prefix)')
  ldsc_prefix=$("$conda_bin" run --name ldsc-gpca-ldsc python -c 'import sys; print(sys.prefix)')
else
  "$manager" env create --prefix "$main_prefix" --file "$repo/environment.yml" --yes
  "$manager" env create --prefix "$ldsc_prefix" --file "$repo/environment.ldsc.yml" --yes
fi
main_run() { "$conda_bin" run --no-capture-output --prefix "$main_prefix" "$@"; }
child_run() { "$conda_bin" run --no-capture-output --prefix "$ldsc_prefix" "$@"; }
main_run Rscript "$repo/scripts/install_genomicsem.R"
main_run python -m pip install "$repo"
echo 'Checking installed tools and package imports...'
main_run python -m pip check
child_run python -m pip check
main_run python -c 'import ldsc_gpca,numpy,pandas,polars; print("ldsc-gpca",ldsc_gpca.__version__)'
child_run python -c 'import numpy,pandas,scipy,bitarray,pysam,pybedtools,nose,ldscore; print("CBIIT LDSC imports OK")'
main_run Rscript -e 'for (p in c("argparse","data.table","glue","GenomicSEM")) {library(p,character.only=TRUE); cat(p,as.character(packageVersion(p)),"\n")}; sessionInfo()'
main_run bcftools --version
main_run bash --version
main_run awk 'BEGIN { print "awk OK"; exit 0 }'
main_run git --version
main_run ldsc-gpca --help
for script in ldsc.py munge_sumstats.py make_annot.py; do child_run "$script" --help; done
"$conda_bin" env config vars set --prefix "$main_prefix" "LDSC_GPCA_LDSC_PREFIX=$ldsc_prefix"
"$conda_bin" list --explicit --prefix "$main_prefix" > "$root/main-conda-explicit.txt"
"$conda_bin" list --explicit --prefix "$ldsc_prefix" > "$root/ldsc-conda-explicit.txt"
main_run python -m pip freeze > "$root/main-pip-freeze.txt"
child_run python -m pip freeze > "$root/ldsc-pip-freeze.txt"
main_run Rscript -e 'write.table(installed.packages()[,c("Package","Version")],stdout(),sep="\t",row.names=FALSE); sessionInfo()' > "$root/R-packages-session.txt"
printf '. %q\nconda activate %q\n' "$base/etc/profile.d/conda.sh" "$main_prefix" > "$root/activate.sh"
echo "All checks passed. Activation file: $root/activate.sh"
echo 'Modified GWAMA v1.2.6 is bundled; reference/data files must still be supplied by the user.'
if [[ $activate == 1 && -t 0 && -t 1 ]]; then
  echo 'Opening an activated Bash shell. Type exit to return to your previous shell.'
  exec bash --rcfile "$root/activate.sh" -i
fi
if [[ $named == 1 ]]; then
  echo 'Activate later with: conda activate ldsc-gpca'
  echo 'The LDSC environment is launched internally; do not activate it yourself.'
  printf 'If Conda is not initialized in this terminal, first run: source %q\n' "$base/etc/profile.d/conda.sh"
  printf 'Optional persistent setup: %q init bash   (use zsh instead for zsh), then reopen your terminal.\n' "$conda_bin"
else
  printf 'Activate later with: source %q\n' "$root/activate.sh"
fi
