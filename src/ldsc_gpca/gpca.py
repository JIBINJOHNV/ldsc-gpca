"""Run the bundled, unchanged genomicPCA R script."""
from importlib.resources import as_file, files
import shutil
import subprocess
import sys


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    rscript = shutil.which("Rscript")
    if rscript is None:
        print("ERROR: Rscript was not found on PATH. Install R and the R packages argparse, data.table and glue.", file=sys.stderr)
        return 127
    script = files("ldsc_gpca").joinpath("r/gpsca_gwama_python_ldsc.r")
    with as_file(script) as path:
        return subprocess.run([rscript, str(path), *argv], check=False).returncode
