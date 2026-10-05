"""Backend selection must preserve each workflow's arguments and input contract."""
import contextlib
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from ldsc_gpca import cli, genomicsem_ldsc, ldsc, ldsc_cli, pipeline, prepare
from ldsc_gpca.utils import DEFAULT_FILTERS


BASE = ['--input', 'traits.csv', '--outdir', 'results', '--ld_ref', 'ld']


class LDSCBackendTests(unittest.TestCase):
    def setUp(self):
        self.output = io.StringIO()
        stack = contextlib.ExitStack()
        self.addCleanup(stack.close)
        for redirect in (contextlib.redirect_stdout, contextlib.redirect_stderr):
            stack.enter_context(redirect(self.output))
        self.python = stack.enter_context(patch.object(ldsc_cli, 'main', return_value=7))
        self.genomicsem = stack.enter_context(patch.object(genomicsem_ldsc, 'main', return_value=0))

    def test_python_default_and_explicit_selection_preserve_arguments_and_status(self):
        arguments = [*BASE, '--ldsc_only', '--munged_dir', 'munged', '--chisq_max', 'auto']
        selections = [arguments, ['--ldsc_backend', 'python', *arguments],
                      [*arguments, '--ldsc_backend=python'],
                      [*BASE, '--ldsc_backend', 'python', *arguments[len(BASE):]]]
        for args in selections:
            with self.subTest(args=args):
                self.assertEqual(cli.main(['ldsc', *args]), 7)
                self.python.assert_called_with(arguments)
                self.genomicsem.assert_not_called()

    def test_genomicsem_raw_vcf_and_both_reuse_modes_select_the_same_workflow(self):
        for mode in (['--hm3', 'hm3.tsv'], ['--vcf_input', '--hm3', 'hm3.tsv'],
                     ['--munged_dir', 'munged'], ['--munged_input']):
            with self.subTest(mode=mode):
                args = [*BASE, *mode, '--chisq_max', '80']
                self.assertEqual(cli.main(['ldsc', *args, '--ldsc_backend=genomicsem']), 0)
                self.genomicsem.assert_called_with(args)
                self.python.assert_not_called()

    def test_existing_genomicsem_command_remains_available(self):
        args = [*BASE, '--vcf_input', '--hm3', 'hm3.tsv']
        self.assertEqual(cli.main(['genomicsem', 'ldsc', *args]), 0)
        self.genomicsem.assert_called_once_with(args)
        self.python.assert_not_called()

    def test_help_and_backend_only_commands_never_start_analysis(self):
        for backend, option, absent, output in (
            ('python', '--info_min', '--info_filter', 'ldsc_results.csv'),
            ('genomicsem', '--info_filter', '--info_min', 'genomicsem_LDSC.RData')):
            for help_flag in ([], ['--help'], ['-h']):
                with self.subTest(backend=backend, help_flag=help_flag):
                    self.output.seek(0); self.output.truncate()
                    args = ['ldsc', *help_flag, '--ldsc_backend', backend]
                    try:
                        code = cli.main(args)
                    except SystemExit as error:
                        code = error.code
                    self.assertEqual(code, 0)
                    text = self.output.getvalue()
                    self.assertIn('usage: ldsc-gpca ldsc ', text)
                    self.assertIn(f'SELECTED BACKEND: {backend}', text)
                    self.assertIn(option, text)
                    self.assertNotIn(absent, text)
                    self.assertIn(output, text)
                    self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_empty_command_keeps_python_help(self):
        self.assertEqual(cli.main(['ldsc']), 0)
        self.assertIn('SELECTED BACKEND: python', self.output.getvalue())
        self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_invalid_missing_abbreviated_and_conflicting_selectors_stop_before_analysis(self):
        for selector in (['--ldsc_backend', 'unknown'], ['--ldsc_backend'],
                         ['--ldsc_backend='], ['--ldsc_back', 'genomicsem'],
                         ['--ldsc-backend', 'python'],
                         ['--ldsc_backend', 'python', '--ldsc_backend=genomicsem']):
            with self.subTest(selector=selector), self.assertRaises(SystemExit) as error:
                cli.main(['ldsc', *BASE, *selector])
            self.assertEqual(error.exception.code, 2)
        self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_consistent_repeated_selector_is_removed_before_dispatch(self):
        self.assertEqual(cli.main(['ldsc', '--ldsc_backend', 'genomicsem', *BASE,
                                   '--ldsc_backend=genomicsem']), 0)
        self.genomicsem.assert_called_once_with(BASE)

    def test_backend_specific_options_cannot_cross_workflows(self):
        for backend, option in (
            ('python', ['--vcf_input']), ('python', ['--munged_input']),
            ('python', ['--info_filter', '.9']), ('python', ['--chisq_max', '1.5']),
            ('genomicsem', ['--ldsc_only']), ('genomicsem', ['--restart']),
            ('genomicsem', ['--info_min', '.9']), ('genomicsem', ['--chisq_max', 'auto'])):
            with self.subTest(backend=backend, option=option), self.assertRaises(SystemExit) as error:
                cli.main(['ldsc', '--ldsc_backend', backend, *BASE, *option])
            self.assertEqual(error.exception.code, 2)
        self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_no_separate_munging_override_on_ldsc_or_pipeline(self):
        for command in ('ldsc', 'pipeline'):
            for backend in ('python', 'genomicsem'):
                for munger in ('python', 'genomicsem'):
                    with self.subTest(command=command, backend=backend, munger=munger), self.assertRaises(SystemExit) as error:
                        cli.main([command, *BASE, '--ldsc_backend', backend, '--munge_backend', munger])
                    self.assertEqual(error.exception.code, 2)
        self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_missing_required_inputs_and_incompatible_modes_stop_before_analysis(self):
        for backend in ('python', 'genomicsem'):
            with self.subTest(backend=backend), self.assertRaises(SystemExit) as error:
                cli.main(['ldsc', '--ldsc_backend', backend, '--input', 'missing.csv'])
            self.assertEqual(error.exception.code, 2)
        with self.assertRaises(SystemExit) as error:
            cli.main(['ldsc', '--ldsc_backend', 'genomicsem', *BASE, '--vcf_input', '--munged_input'])
        self.assertEqual(error.exception.code, 2)
        self.python.assert_not_called(); self.genomicsem.assert_not_called()

    def test_help_copy_does_not_change_backend_or_pipeline_parsers(self):
        for backend in ('genomicsem', 'python', 'genomicsem', 'python'):
            parser = ldsc.build_parser(backend)
            opts = vars(parser.parse_args([*BASE, '--ldsc_backend', backend]))
            opts.pop('ldsc_backend')
            self.assertEqual(opts, vars(ldsc.ldsc_parser(backend).parse_args(BASE)))
        self.assertNotIn('--ldsc_backend', ldsc_cli.parser._option_string_actions)
        self.assertEqual(ldsc_cli.parser.prog, 'ldsc-gpca ldsc')
        for backend, workers in (('python', 5), ('genomicsem', 1)):
            opts = pipeline.build_parser(backend).parse_args(BASE)
            self.assertEqual(opts.ldsc_backend, backend)
            self.assertEqual(opts.n_cores, workers)

    def test_prepare_retains_munging_selection_and_omission_default(self):
        parser = prepare.build_parser()
        base = ['--input', 'traits.csv', '--outdir', 'out', '--mode', 'ldsc']
        self.assertEqual(parser.parse_args(base).munge_backend, 'python')
        self.assertEqual(parser.parse_args([*base, '--munge_backend', 'genomicsem']).munge_backend, 'genomicsem')
        with self.assertRaises(SystemExit) as error:
            parser.parse_args([*base, '--ldsc_backend', 'genomicsem'])
        self.assertEqual(error.exception.code, 2)


class MungedDirectorySafetyTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.supplied = self.root / 'supplied'; self.supplied.mkdir()
        self.out = self.root / 'out'
        self.manifest = self.root / 'traits.csv'
        self.vcf = self.root / 'A.vcf'
        self.vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        self.manifest.write_text('traitname,ref,population_prevalence,sample_prevalence,vcf_files\n'
                                 f'A,yes,,,{self.vcf}\nB,yes,,,{self.vcf}\n')
        self.ld = self.root / 'reference'; self.ld.mkdir()
        self.hm3 = self.root / 'hm3.tsv'; self.hm3.write_text('SNP\tA1\tA2\nrs1\tA\tG\n')
        for name in ('A', 'B'):
            target = self.supplied / f'{name}.sumstats.gz'
            target.write_bytes(gzip.compress(('SNP\tA1\tA2\tN\tZ\n' + ''.join(
                f'rs{i}\tA\tG\t10000\t{i/10}\n' for i in range(10))).encode(), mtime=0))
            (self.supplied / f'{name}.prevalence.json').write_text(json.dumps({
                'sample_size_column': 'NEF', 'sample_prevalence': None,
                'filters': DEFAULT_FILTERS, 'sha256': hashlib.sha256(target.read_bytes()).hexdigest()}))
        (self.supplied / 'keep.txt').write_text('Unrelated supplied file.\n')
        self.output = io.StringIO()
        stack = contextlib.ExitStack(); self.addCleanup(stack.close)
        for redirect in (contextlib.redirect_stdout, contextlib.redirect_stderr):
            stack.enter_context(redirect(self.output))

    def args(self):
        return ['--input', str(self.manifest), '--outdir', str(self.out),
                '--ld_ref', str(self.ld), '--n_cores', '1']

    def snapshot(self, folder):
        return {str(path.relative_to(folder)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in folder.rglob('*') if path.is_file()}

    def test_invalid_combination_preserves_all_hashes_and_stops_before_any_work(self):
        stages = ('read_manifest', 'check_runtime', 'prepare_compilation_outputs',
                  'run_vcf_to_table', 'parallel_munge_sumstats',
                  'parallel_ldsc_analysis', 'compile_results')
        for route in ('direct', 'default', 'explicit_python'):
            for with_hm3 in (True, False):
                for existing_output in (False, True):
                    with self.subTest(route=route, hm3=with_hm3, existing_output=existing_output):
                        self.out = self.root / f'{route}_{with_hm3}_{existing_output}'
                        if existing_output:
                            self.out.mkdir()
                            (self.out / 'ldsc_results.csv').write_text('Previous results must remain intact.\n')
                        args = self.args() + ['--munged_dir', str(self.supplied)]
                        if with_hm3: args += ['--hm3', str(self.hm3)]
                        if existing_output: args += ['--restart']
                        before = self.snapshot(self.root)
                        self.output.seek(0); self.output.truncate()
                        with contextlib.ExitStack() as stack:
                            mocks = [stack.enter_context(patch.object(ldsc_cli, stage,
                                side_effect=AssertionError('Unexpected stage: ' + stage))) for stage in stages]
                            with self.assertRaises(SystemExit) as error:
                                if route == 'direct': ldsc_cli.main(args)
                                else: cli.main(['ldsc', *(['--ldsc_backend', 'python']
                                                       if route == 'explicit_python' else []), *args])
                            self.assertEqual(error.exception.code, 2)
                            for mock in mocks: mock.assert_not_called()
                        self.assertIn('--munged_dir requires --ldsc_only', self.output.getvalue())
                        self.assertEqual(before, self.snapshot(self.root))
                        self.assertEqual(self.out.exists(), existing_output)

    def test_valid_reuse_preserves_inputs_in_explicit_and_default_directories(self):
        for explicit in (True, False):
            with self.subTest(explicit=explicit):
                self.out = self.root / f'reuse_{explicit}'
                directory = self.supplied
                if not explicit:
                    directory = self.out / 'ldsc_input'
                    shutil.copytree(self.supplied, directory)
                before = self.snapshot(directory)
                args = self.args() + ['--ldsc_only']
                if explicit: args += ['--munged_dir', str(directory)]
                with patch.object(ldsc_cli, 'check_runtime'), \
                     patch.object(ldsc_cli, 'run_vcf_to_table', side_effect=AssertionError('extraction')), \
                     patch.object(ldsc_cli, 'parallel_munge_sumstats', side_effect=AssertionError('munging')), \
                     patch.object(ldsc_cli, 'parallel_ldsc_analysis', return_value=[]) as regression, \
                     patch.object(ldsc_cli, 'compile_results', return_value='fixture.csv'):
                    cli.main(['ldsc', *args])
                self.assertEqual(regression.call_args.args[5], str(directory))
                self.assertEqual(before, self.snapshot(directory))

    def test_fresh_run_keeps_default_destination(self):
        before = self.snapshot(self.supplied)
        with patch.object(ldsc_cli, 'check_runtime'), \
             patch.object(ldsc_cli, 'run_vcf_to_table', side_effect=lambda n, frame, *args: frame) as extraction, \
             patch.object(ldsc_cli, 'parallel_munge_sumstats', return_value=['A', 'B']) as munging, \
             patch.object(ldsc_cli, 'parallel_ldsc_analysis', return_value=[]) as regression, \
             patch.object(ldsc_cli, 'compile_results', return_value='fixture.csv'):
            cli.main(['ldsc', *self.args(), '--hm3', str(self.hm3)])
        extraction.assert_called_once(); munging.assert_called_once(); regression.assert_called_once()
        self.assertEqual(munging.call_args.kwargs['ldsc_input_folder'], str(self.out / 'ldsc_input'))
        self.assertEqual(before, self.snapshot(self.supplied))


if __name__ == '__main__':
    unittest.main()
