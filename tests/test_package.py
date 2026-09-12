import contextlib
import io
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from importlib.resources import files
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from ldsc_gpca import cli, gpca, ldsc_cli, utils


class PackageTests(unittest.TestCase):
    def test_top_level_help_and_version(self):
        for args, expected in [([], 'gpca'), (['--help'], 'ldsc'), (['--version'], '0.1.0')]:
            result = subprocess.run([sys.executable, '-m', 'ldsc_gpca', *args], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(expected, result.stdout)

    def test_ldsc_help(self):
        result = subprocess.run([sys.executable, '-m', 'ldsc_gpca', 'ldsc', '--help'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('--ldsc-retries', result.stdout)

    def test_ldsc_dispatch_preserves_arguments(self):
        args = ['--input_file', '/a path/manifest.csv', '--ldsc-retries', '2']
        with patch.object(ldsc_cli, 'main', return_value=0) as run:
            self.assertEqual(cli.main(['ldsc', *args]), 0)
            run.assert_called_once_with(args)

    def test_gpca_arguments_exit_and_packaged_resource(self):
        with patch.object(gpca.shutil, 'which', return_value='/an R path/Rscript'), patch.object(gpca.subprocess, 'run') as run:
            run.return_value.returncode = 7
            self.assertEqual(cli.main(['gpca', '--input', '/a path/traits.csv']), 7)
            argv = run.call_args.args[0]
            self.assertEqual(argv[0], '/an R path/Rscript')
            self.assertTrue(Path(argv[1]).is_file())
            self.assertEqual(argv[2:], ['--input', '/a path/traits.csv'])
            self.assertNotIn('shell', run.call_args.kwargs)

    def test_missing_rscript(self):
        with patch.object(gpca.shutil, 'which', return_value=None), contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(gpca.main([]), 127)
            self.assertIn('Rscript was not found', err.getvalue())

    def test_r_script_is_unchanged(self):
        root = Path(__file__).resolve().parents[1]
        self.assertEqual(files('ldsc_gpca').joinpath('r/gpsca_gwama_python_ldsc.r').read_bytes(),
                         (root/'gpsca_gwama_python_ldsc.r').read_bytes())

    def test_bad_subcommand_and_missing_arguments(self):
        for argv in [['unknown'], ['ldsc']]:
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                cli.main(argv)
            self.assertEqual(error.exception.code, 2)

    def test_invalid_retry_cli(self):
        args = ['--input_file', 'missing.csv', '--output_folder', '/tmp/unused',
                '--ld_ref_snp_file', 'missing', '--ld_ref', 'missing', '--ldsc-retries', '-1']
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            ldsc_cli.main(args)
        self.assertEqual(error.exception.code, 2)

    def test_explicit_error_log_path(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(utils.subprocess, 'run') as run:
            run.return_value = subprocess.CompletedProcess('test', 1, '', 'simulated failure')
            with self.assertRaises(RuntimeError):
                utils.run_command('test', 'attempt 1/2', output_folder=folder)
            self.assertIn('simulated failure', (Path(folder)/'execution_errors.log').read_text())

    def test_spawn_import(self):
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context('spawn')) as pool:
            self.assertEqual(pool.submit(utils.optional_prevalence, 0.2, 'test').result(timeout=30), 0.2)


if __name__ == '__main__':
    unittest.main()
