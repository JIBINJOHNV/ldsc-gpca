"""Launch LDSC in an isolated environment without changing the caller's environment."""
import os
import shutil
import subprocess
import sys


RAW_LDSC_SCRIPTS = frozenset({'ldsc.py', 'munge_sumstats.py'})


def ldsc_command(script, *, conda='conda', environment='ldsc-cbiit', prefix=None):
    target = ['--prefix', prefix] if prefix else ['--name', environment]
    return [conda, 'run', '--no-capture-output', *target, script]


def run_ldsc_script(script, arguments, *, conda=None, environment=None, prefix=None):
    """Run an allow-listed CBIIT LDSC console script with unchanged arguments."""
    if script not in RAW_LDSC_SCRIPTS:
        raise ValueError(f'Unsupported raw LDSC script: {script}')
    conda = conda or os.environ.get('CONDA_EXE', 'conda')
    environment = environment or os.environ.get('LDSC_GPCA_LDSC_ENV', 'ldsc-cbiit')
    prefix = prefix if prefix is not None else os.environ.get('LDSC_GPCA_LDSC_PREFIX')
    command = ldsc_command(script, conda=conda, environment=environment, prefix=prefix)
    try:
        return subprocess.run(command + list(arguments)).returncode
    except OSError as error:
        print(f'Unable to launch {script} through Conda: {error}', file=sys.stderr)
        return 127


def check_runtime(*, conda='conda', environment='ldsc-cbiit', prefix=None, bcftools=None):
    if not shutil.which(conda):
        raise ValueError(f'Conda executable not found: {conda}; use --conda-executable.')
    if bcftools is not None:
        for tool in (bcftools, 'bash', 'awk'):
            if not shutil.which(tool):
                raise ValueError(f'Required extraction executable not found: {tool}')
    for script in ('ldsc.py', 'munge_sumstats.py'):
        result = subprocess.run(ldsc_command(script, conda=conda, environment=environment, prefix=prefix) + ['--help'],
                                capture_output=True, text=True)
        if result.returncode:
            raise ValueError(f'{script} failed in the LDSC environment. Run scripts/setup_environments.sh.\n{result.stderr}')
