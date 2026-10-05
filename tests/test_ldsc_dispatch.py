"""Backend selection must preserve each workflow's arguments and input contract."""
import contextlib
import io
import unittest
from unittest.mock import patch

from ldsc_gpca import cli, genomicsem_ldsc, ldsc, ldsc_cli, pipeline, prepare


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


if __name__ == '__main__':
    unittest.main()
