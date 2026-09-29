"""CLI startup, subprocess inheritance and single worker-option names."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from ldsc_gpca.threads import THREAD_ENV_VARS, configure_numerical_threads


ROOT = Path(__file__).resolve().parents[1]


def clean_environment():
    env = {name: value for name, value in os.environ.items() if name not in THREAD_ENV_VARS}
    env['PYTHONPATH'] = str(ROOT / 'src')
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return env


class ThreadDefaultsTests(unittest.TestCase):
    def test_defaults_preserve_unrelated_environment_and_are_idempotent(self):
        with patch.dict(os.environ, {'UNRELATED_SETTING': 'keep'}, clear=True):
            expected = dict.fromkeys(THREAD_ENV_VARS, '1')
            self.assertEqual(configure_numerical_threads(), expected)
            self.assertEqual(configure_numerical_threads(), expected)
            self.assertEqual(os.environ['UNRELATED_SETTING'], 'keep')

    def test_explicit_values_are_preserved_and_reported(self):
        # Existing values belong to the caller; even invalid/empty values are
        # left to the numerical library rather than silently replaced.
        for value in ('1', '4', '', 'invalid'):
            with self.subTest(value=value), patch.dict(os.environ, {'OMP_NUM_THREADS': value}, clear=True):
                output = io.StringIO()
                with contextlib.redirect_stderr(output):
                    settings = configure_numerical_threads(report=True)
                self.assertEqual(settings['OMP_NUM_THREADS'], value)
                self.assertTrue(all(settings[name] == '1' for name in THREAD_ENV_VARS if name != 'OMP_NUM_THREADS'))
                self.assertIn('OMP_NUM_THREADS=' + value, output.getvalue())

    def test_importing_cli_does_not_change_environment_or_load_numerical_libraries(self):
        code = '''import os, sys
before = dict(os.environ)
import ldsc_gpca.cli
assert dict(os.environ) == before
assert not {'numpy', 'pandas', 'polars'}.intersection(sys.modules)
'''
        result = subprocess.run([sys.executable, '-c', code], env=clean_environment(), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_defaults_are_set_before_numerical_imports_for_every_managed_workflow(self):
        code = '''import sys, os
from ldsc_gpca.threads import THREAD_ENV_VARS
seen = []
def check_import(event, args):
    if event == 'import' and args[0] in ('numpy', 'pandas', 'polars'):
        assert all(os.environ.get(name) == '1' for name in THREAD_ENV_VARS)
        seen.append(args[0])
sys.addaudithook(check_import)
from ldsc_gpca.cli import main
try:
    status = main(sys.argv[1:])
except SystemExit as error:
    status = error.code
assert status == 0
if sys.argv[1:3] != ['genomicsem', 'ldsc']:
    assert seen, 'No numerical import was checked'
'''
        for command in (['ldsc'], ['gpca'], ['prepare'], ['genomicsem', 'gpca'], ['genomicsem', 'ldsc']):
            with self.subTest(command=command):
                result = subprocess.run([sys.executable, '-c', code, *command, '--help'],
                                        env=clean_environment(), capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('--n_cores', result.stdout)
                self.assertNotIn('Numerical thread environment:', result.stderr)

    def test_raw_ldsc_child_inherits_defaults_and_overrides(self):
        with tempfile.TemporaryDirectory() as folder:
            fake = Path(folder) / 'conda'
            fake.write_text('#!' + sys.executable + '\nimport json, os, sys\n'
                            'print(json.dumps({k: os.environ[k] for k in ' + repr(THREAD_ENV_VARS) + '}))\n'
                            'sys.exit(7)\n')
            fake.chmod(0o755)
            for script in ('ldsc.py', 'munge_sumstats.py'):
                for overrides in ({}, {'OPENBLAS_NUM_THREADS': '2'}):
                    with self.subTest(script=script, overrides=overrides):
                        env = {**clean_environment(), 'CONDA_EXE': str(fake), **overrides}
                        result = subprocess.run([sys.executable, '-m', 'ldsc_gpca', script, '--example'],
                                                env=env, capture_output=True, text=True)
                        self.assertEqual(result.returncode, 7)
                        self.assertEqual(json.loads(result.stdout), {**dict.fromkeys(THREAD_ENV_VARS, '1'), **overrides})
                        self.assertIn('Numerical thread environment:', result.stderr)

    def test_gpca_r_child_inherits_environment_and_receives_worker_option(self):
        with tempfile.TemporaryDirectory() as folder:
            fake = Path(folder) / 'Rscript'
            fake.write_text('#!' + sys.executable + '\nimport json, os, sys\n'
                            'print("CHILD=" + json.dumps({"env": {k: os.environ[k] for k in ' +
                            repr(THREAD_ENV_VARS) + '}, "args": sys.argv[1:]}))\n')
            fake.chmod(0o755)
            env = {**clean_environment(), 'PATH': folder + os.pathsep + os.environ.get('PATH', '')}
            for command in (['gpca'], ['genomicsem', 'gpca']):
                for option in ('--n_cores',):
                    with self.subTest(command=command, option=option):
                        arguments = [*command, '--input', 'traits.csv', '--outdir', folder, '--validate_only', option, '22']
                        result = subprocess.run([sys.executable, '-m', 'ldsc_gpca', *arguments],
                                                env=env, capture_output=True, text=True)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        line = next(line for line in result.stdout.splitlines() if line.startswith('CHILD='))
                        child = json.loads(line.removeprefix('CHILD='))
                        self.assertEqual(child['env'], dict.fromkeys(THREAD_ENV_VARS, '1'))
                        self.assertEqual(child['args'][-2:], [option, '22'])


class PythonWorkerOptionsTests(unittest.TestCase):
    def test_genomicsem_ldsc_default_and_name(self):
        from ldsc_gpca.genomicsem_ldsc import build_parser
        base = ['--input', 'traits.csv', '--outdir', 'out', '--ld_ref', 'ld']
        parser = build_parser()
        self.assertEqual(parser.parse_args(base).n_cores, 1)
        for flag in ('--n_cores',):
            with self.subTest(flag=flag):
                self.assertEqual(parser.parse_args([*base, flag, '22']).n_cores, 22)
                for value in ('bad', None):
                    with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                        parser.parse_args([*base, flag, *([] if value is None else [value])])
                    self.assertEqual(error.exception.code, 2)

    def test_genomicsem_ldsc_rejects_zero_and_negative_before_reading_inputs(self):
        from ldsc_gpca.genomicsem_ldsc import main
        for flag in ('--n_cores',):
            for value in ('0', '-1'):
                output = io.StringIO()
                with self.subTest(flag=flag, value=value), contextlib.redirect_stderr(output), self.assertRaises(SystemExit) as error:
                    main(['--input', 'missing.csv', '--outdir', 'out', '--ld_ref', 'ld', flag, value])
                self.assertEqual(error.exception.code, 2)
                self.assertIn('--n_cores must be >=1', output.getvalue())

    def test_standalone_prepare_workers_and_embedded_stage_are_separate(self):
        from ldsc_gpca.prepare import main
        from ldsc_gpca.gpca import postprocess_parser
        for flag in ('--n_cores',):
            with self.subTest(flag=flag), patch('ldsc_gpca.prepare.prepare_inputs') as prepare:
                self.assertEqual(main(['--input', 'traits.csv', '--outdir', 'out', flag, '3']), 0)
                self.assertEqual(prepare.call_args.kwargs['prepare_workers'], 3)
        opts, remaining = postprocess_parser(include_prepare=True).parse_known_args(
            ['--prepare_workers', '4', '--n_cores', '22'])
        self.assertEqual(opts.prepare_workers, 4)
        self.assertEqual(remaining, ['--n_cores', '22'])
        for flag in ('--n_cores',):
            for value in ('0', '-1'):
                output = io.StringIO()
                with contextlib.redirect_stderr(output):
                    status = main(['--input', 'missing.csv', '--outdir', 'out', flag, value])
                self.assertEqual(status, 1)
                self.assertIn('must be a positive integer', output.getvalue())


class RWorkerOptionsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rscript = os.environ.get('LDSC_GPCA_TEST_RSCRIPT') or shutil.which('Rscript')
        if not cls.rscript:
            raise unittest.SkipTest('Rscript is not available')
        result = subprocess.run([cls.rscript, '-e', 'quit(status=if(requireNamespace("argparse", quietly=TRUE)) 0 else 1)'],
                                capture_output=True, text=True)
        if result.returncode:
            raise unittest.SkipTest('R argparse is not available: ' + result.stderr)

    def test_real_r_parsers_accept_single_name_and_reject_invalid_worker_counts(self):
        # Source only parser/workflow definitions; no GenomicSEM or GWAMA data
        # are needed to validate options. Negative values fail before analysis.
        for backend in ('python_ldsc', 'genomicsem'):
            with tempfile.TemporaryDirectory() as folder:
                manifest = Path(folder) / 'traits.csv'
                ldsc = Path(folder) / 'ldsc.csv'
                manifest.touch()
                ldsc.touch()
                module = ROOT / 'src' / 'ldsc_gpca' / 'r' / backend
                code = 'suppressPackageStartupMessages(library(argparse)); bundled_gwama_path <- "unused"; '
                code += f'source({json.dumps(str(module.parent / "shared" / "input.R"))}); '
                code += f'source({json.dumps(str(module / "cli.R"))}); '
                if backend == 'python_ldsc':
                    code += 'x <- parse_command_line(); x <- validate_cli_paths(x); cat("WORKERS=", x$n_cores, "\\n", sep="")'
                else:
                    code += f'source({json.dumps(str(module / "workflow.R"))}); '
                    code += 'x <- genomicsem_parser()$parse_args(); if(x$n_cores < 0L) genomicsem_main() else cat("WORKERS=", x$n_cores, "\\n", sep="")'
                ldsc_flag = '--ldsc_results'
                base = ['--input', str(manifest), ldsc_flag, str(ldsc), '--outdir', folder, '--validate_only']
                cases = [([], 0)]
                for flag in ('--n_cores',):
                    cases.extend(([flag, value], expected) for value, expected in
                                 [('0', 0), ('1', 1), ('22', 22), ('100', 100), ('-1', None), ('bad', None)])
                    cases.append(([flag], None))
                for arguments, expected in cases:
                    with self.subTest(backend=backend, arguments=arguments):
                        result = subprocess.run([self.rscript, '-e', code, *base, *arguments],
                                                env=clean_environment(), capture_output=True, text=True)
                        if expected is None:
                            self.assertNotEqual(result.returncode, 0)
                            if arguments[-1] == '-1':
                                self.assertIn('--n_cores must be', result.stderr)
                        else:
                            self.assertEqual(result.returncode, 0, result.stderr)
                            self.assertIn(f'WORKERS={expected}\n', result.stdout)


if __name__ == '__main__':
    unittest.main()
