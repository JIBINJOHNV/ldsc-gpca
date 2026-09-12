import os
import math
import gzip
import subprocess
import pandas as pd

def optional_prevalence(value, label):
    if pd.isna(value) or str(value).strip().lower() in ('', '.', 'na', 'nan', 'none'):
        return None
    number = float(value)
    if not math.isfinite(number) or not 0 < number < 1:
        raise ValueError(f"{label} must be finite and strictly between 0 and 1")
    return number

DEFAULT_FILTERS = {
    'exclude_mhc': False, 'mhc_chr': '6', 'mhc_start': 25000000, 'mhc_end': 35000000,
    'info_min': 0.7, 'maf_min': 0.01, 'munge_maf_min': 0.005, 'max_af_difference': 0.2,
    'remove_palindrome': False, 'pal_lower': 0.45, 'pal_upper': 0.55,
}
def run_command(command, error_msg="Command failed", *, output_folder):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # Save error to a log file instead of flooding screen
        log_path = os.path.join(output_folder, "execution_errors.log")
        with open(log_path, "a") as f:
            f.write(f"--- ERROR: {error_msg} ---\n{result.stderr}\n")
        raise RuntimeError(error_msg)

def is_valid_gz(filepath):
    if not os.path.exists(filepath) or os.path.getsize(filepath) < 50:
        return False
    try:
        with gzip.open(filepath, 'rt') as f:
            f.read(1)
        return True
    except: return False
