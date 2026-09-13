"""Launch LDSC in an isolated environment without changing the caller's environment."""
import shutil
import subprocess

def ldsc_command(script, *, conda='conda', environment='ldsc-cbiit', prefix=None):
    target = ['--prefix', prefix] if prefix else ['--name', environment]
    return [conda, 'run', '--no-capture-output', *target, script]

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
