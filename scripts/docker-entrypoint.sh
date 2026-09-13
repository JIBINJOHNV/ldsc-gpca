#!/usr/bin/env bash
set -e
# Accept Nextflow's /bin/bash task command as well as the short package CLI.
source /opt/environments/activate.sh
if [[ $# == 0 ]]; then set -- ldsc-gpca --help; fi
case "$1" in
  prepare|ldsc|gpca|genomicsem|-*) set -- ldsc-gpca "$@" ;;
esac
exec "$@"
