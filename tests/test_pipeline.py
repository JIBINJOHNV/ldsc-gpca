"""Pipeline contracts, routing, trait order, failure isolation and real preparation."""
import contextlib
import csv
import gzip
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

from ldsc_gpca import pipeline
from ldsc_gpca.interfaces import read_manifest


def csv_file(path, columns, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.manifest = self.root/'traits.csv'
        self.out = self.root/'output'
        self.ld = self.root/'ld'
        self.ld.mkdir()
        for chrom in range(1, 23):
            (self.ld/f'{chrom}.l2.ldscore.gz').write_bytes(b'reference fixture')
            (self.ld/f'{chrom}.l2.M_5_50').write_text('100\n')
        self.hm3 = self.root/'hm3.txt'
        self.hm3.write_text('SNP A1 A2\nrs1 A G\nrs2 C T\n')
        self.vcf = self.root/'input.vcf'
        self.vcf.write_text('VCF fixture; replaced for real-tool test\n')
        self.raw = self.root/'raw.tsv'
        self.raw.write_text('SNP\tA1\tA2\tN\tZ\tP\nrs1\tG\tA\t1000\t2\t0.0455\n')
        self.munged = self.root/'munged'
        self.munged.mkdir()
        for name in ('002', '001', '003'):
            with gzip.open(self.munged/f'{name}.sumstats.gz', 'wt') as stream:
                stream.write('SNP\tA1\tA2\tN\tZ\nrs1\tG\tA\t1000\t2\n')
        self.prepared = self.root/'prepared'
        self.prepared.mkdir()
        self.columns = ['traitname', 'vcf_files', 'sample_prevalence', 'population_prevalence']
        self.rows = [[name, 'input.vcf', 'NA', 'NA'] for name in ('002', '001', '003')]
        self.save_manifest()
        self.calls = []
        for stream in ('stdout', 'stderr'):
            manager = getattr(contextlib, 'redirect_' + stream)(io.StringIO())
            manager.__enter__()
            self.addCleanup(manager.__exit__, None, None, None)

    def save_manifest(self):
        csv_file(self.manifest, self.columns, self.rows)

    def args(self, backend='python', *extra):
        return ['--input', str(self.manifest), '--outdir', str(self.out), '--ld_ref', str(self.ld),
                '--ldsc_backend', backend, '--hm3', str(self.hm3), '--n_cores', '2',
                '--prepare_workers', '1', '--gwama_output_info', '.9', *extra]

    def fake_run(self, command, **kwargs):
        args = command[3:]
        self.calls.append(args)
        out = Path(args[args.index('--outdir') + 1])
        out.mkdir(parents=True)
        name = args[1] if args[0] == 'genomicsem' else args[0]
        manifest = Path(args[args.index('--input') + 1])
        columns, rows = read_manifest(manifest)
        if name == 'prepare':
            (out/'gpca_inputs').mkdir()
            if '--write_munge_inputs' in args:
                (out/'munge_inputs').mkdir()
                for row in rows:
                    (out/'munge_inputs'/f'{row["traitname"]}_munge_inputs.txt').write_text(self.raw.read_text())
        elif name == 'ldsc':
            if args[0] == 'genomicsem':
                (out/'genomicPCA_LDSC.RData').write_bytes(b'native result fixture')
                # Deliberately shuffled: pipeline must restore original row order.
                csv_file(out/'Selected_Traits.csv', ['traitname', 'traits', 'sampleprevalence'],
                         [['001', 'file1', 'NA'], ['002', 'file2', 'NA']])
            else:
                pairs = []
                for a in rows:
                    for b in rows:
                        rg = 1 if a['traitname'] == b['traitname'] else .2
                        pairs.append([a['traitname'], b['traitname'], rg, .1, rg/.1, .01, .2, .02, 1, .01, .05, .01])
                csv_file(out/'ldsc_results.csv', ['p1','p2','rg','se','z','p','h2_obs','h2_obs_se','h2_int','h2_int_se','gcov_int','gcov_int_se'], pairs)
        elif name == 'gpca':
            (out/'GenomicPCA_PC1_Weights_Used.csv').write_text('synthetic stage output\n')
            if '--validate_only' not in args:
                dataset = args[args.index('--dataset_id') + 1]
                (out/'GWAMA_Run_Status.csv').write_text('synthetic stage output\n')
                (out/f'{dataset}_postprocess.json').write_text('{}')
                (out/f'{dataset}_GWAMA_combined_results.txt.gz').touch()
                (out/'harmonisation_input').mkdir()
                (out/'harmonisation_input'/f'{dataset}_GPCA_inputs.txt.gz').touch()
        else:
            raise AssertionError(command)
        return subprocess.CompletedProcess(command, 0)

    def run_mock(self, args):
        with patch.object(pipeline.subprocess, 'run', side_effect=self.fake_run):
            return pipeline.main(args)

    def test_python_vcf_full_pipeline_and_cutoff_forwarding(self):
        self.assertEqual(self.run_mock(self.args('python', '--chisq_max', '80', '--gwama_output_n_eff', '33000')), 0)
        self.assertEqual([args[0] for args in self.calls], ['prepare', 'ldsc', 'gpca'])
        ldsc = self.calls[1]
        self.assertEqual(ldsc[ldsc.index('--chisq_max')+1], '80')
        self.assertNotIn('--chisq_max', self.calls[0]); self.assertNotIn('--chisq_max', self.calls[2])
        columns, rows = read_manifest(self.out/'manifests/ldsc_traits.csv')
        self.assertEqual([row['traitname'] for row in rows], ['002', '001', '003'])
        self.assertTrue(all(row['ref'] == 'yes' for row in rows))
        self.assertTrue(all(Path(row['vcf_files']).is_absolute() for row in rows))
        self.assertEqual(read_manifest(self.out/'manifests/gpca_traits.csv')[0], ['traitname'])
        report = json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['status'], 'completed')
        self.assertEqual(report['chisq_policy'], 'explicit_per_trait')
        self.assertTrue(all(stage['status'] == 'completed' for stage in report['stages']))
        self.assertEqual(self.calls[-1][self.calls[-1].index('--dataset_id')+1], 'output')

    def test_native_unmunged_handoff_and_cutoff(self):
        self.columns.append('sumstats_file')
        for row in self.rows: row.append('raw.tsv')
        self.save_manifest()
        self.assertEqual(self.run_mock(self.args('genomicsem', '--chisq_max', '80')), 0)
        self.assertEqual([args[:2] for args in self.calls[1:]], [['genomicsem','ldsc'],['genomicsem','gpca']])
        _, selected = read_manifest(self.out/'manifests/gpca_traits.csv')
        self.assertEqual(selected, [{'traitname':'002'}, {'traitname':'001'}])
        self.assertEqual(self.calls[1][self.calls[1].index('--chisq_max')+1], '80.0')
        self.assertNotIn('--chisq_max', self.calls[-1])
        report = json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['removed_before_gpca'], ['003'])

    def test_native_quantitative_vcf_uses_existing_preparation_exports(self):
        self.assertEqual(self.run_mock(self.args('genomicsem')), 0)
        self.assertIn('--write_munge_inputs', self.calls[0])
        selected_hm3 = Path(self.calls[0][self.calls[0].index('--hm3')+1])
        self.assertEqual(selected_hm3.read_text(), 'SNP\nrs1\nrs2\n')
        self.assertEqual(self.calls[1][self.calls[1].index('--hm3')+1], str(self.hm3))
        _, rows = read_manifest(self.out/'manifests/ldsc_traits.csv')
        for row in rows:
            self.assertEqual(row['sumstats_file'], str(self.out/'prepare/munge_inputs'/f'{row["traitname"]}_munge_inputs.txt'))
        self.assertNotIn('--chisq_max', self.calls[1])

    def test_existing_munged_and_prepared_files_both_backends(self):
        self.columns.remove('vcf_files')
        for row in self.rows:
            del row[1]
        self.save_manifest()
        for backend, flags in [('python',['--ldsc_only']), ('genomicsem',[])]:
            with self.subTest(backend=backend):
                self.out = self.root/backend
                self.calls.clear()
                self.assertEqual(self.run_mock(self.args(backend, '--munged_dir', str(self.munged),
                    '--gpca_input_folder', str(self.prepared), *flags)), 0)
                self.assertEqual(len(self.calls), 2)
                self.assertNotIn('--hm3', self.calls[0])
                self.assertIn('--munged_dir', self.calls[0])
                self.assertEqual(self.calls[-1][self.calls[-1].index('--gpca_input_folder')+1], str(self.prepared))

    def test_native_explicit_munged_paths(self):
        self.columns.append('munged_file')
        for row in self.rows: row.append(f'munged/{row[0]}.sumstats.gz')
        self.save_manifest()
        self.assertEqual(self.run_mock(self.args('genomicsem', '--munged_input', '--validate_only')), 0)
        self.assertEqual(len(self.calls), 2)
        self.assertIn('--munged_input', self.calls[0]); self.assertNotIn('--hm3', self.calls[0])

    def test_native_explicit_vcf_is_converted_in_ldsc_not_gwama_preparation(self):
        self.assertEqual(self.run_mock(self.args('genomicsem', '--vcf_input', '--validate_only')), 0)
        self.assertEqual([args[:2] for args in self.calls], [['genomicsem','ldsc'], ['genomicsem','gpca']])
        self.assertIn('--vcf_input', self.calls[0])
        self.assertIn('--hm3', self.calls[0])
        self.assertNotIn('--vcf_input', self.calls[1])
        columns, rows = read_manifest(self.out/'manifests/ldsc_traits.csv')
        self.assertIn('vcf_files', columns)
        self.assertNotIn('sumstats_file', columns)

    def test_native_binary_vcf_n_validated_before_running_stages(self):
        for row in self.rows:
            row[2:4] = ['.5', '.1']
        self.save_manifest()
        with patch.object(pipeline.subprocess, 'run') as run:
            self.assertEqual(pipeline.main(self.args('genomicsem','--vcf_input','--validate_only')), 1)
            run.assert_not_called()
        self.assertFalse(self.out.exists())
        self.columns.append('N')
        for row in self.rows:
            row.append('35000')
        self.save_manifest()
        self.assertEqual(self.run_mock(self.args('genomicsem','--vcf_input',
            '--gpca_input_folder',str(self.prepared))),0)
        self.assertEqual(len(self.calls),2)
        self.assertIn('--vcf_input',self.calls[0])

    def test_native_vcf_mode_cannot_be_combined_with_munged_modes(self):
        for flags in (['--munged_input'], ['--munged_dir', str(self.munged)]):
            with self.subTest(flags=flags), self.assertRaises(SystemExit):
                pipeline.main(self.args('genomicsem','--vcf_input',*flags))

    def test_validation_only_runs_ldsc_and_pca_without_gwama_preparation(self):
        args = self.args('python', '--validate_only')
        at = args.index('--gwama_output_info'); del args[at:at+2]
        self.assertEqual(self.run_mock(args), 0)
        self.assertEqual([args[0] for args in self.calls], ['ldsc', 'gpca'])
        self.assertIn('--validate_only', self.calls[-1])
        self.assertNotIn('--gpca_input_folder', self.calls[-1])

    def test_invalid_chisq_and_workers_fail_before_stages(self):
        for extra in (['--chisq_max','0'], ['--chisq_max','-1'], ['--chisq_max','nan'],
                      ['--chisq_max','inf'], ['--n_cores','0'], ['--prepare_workers','0']):
            with self.subTest(extra=extra), patch.object(pipeline.subprocess,'run') as run:
                args = self.args()
                if extra[0] in args:
                    at=args.index(extra[0]); del args[at:at+2]
                try:
                    code = pipeline.main(args+extra)
                except SystemExit as error:
                    code = error.code
                self.assertIn(code, (1, 2))
                run.assert_not_called(); self.assertFalse(self.out.exists())

    def test_python_auto_cutoff_and_trait_wide_normalization_forwarding(self):
        self.assertEqual(self.run_mock(self.args('python', '--chisq_max', 'auto',
            '--rg_normalization', 'trait_wide')), 0)
        ldsc, gpca = self.calls[1:]
        self.assertEqual(ldsc[ldsc.index('--chisq_max')+1], 'auto')
        self.assertEqual(gpca[gpca.index('--rg_normalization')+1], 'trait_wide')
        self.assertNotIn('--rg_normalization', ldsc)
        report = json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['chisq_policy'], 'automatic_per_trait')

    def test_python_dropped_traits_use_the_retained_manifest_in_original_order(self):
        def run(command, **kwargs):
            result = self.fake_run(command, **kwargs)
            if command[3] == 'ldsc':
                csv_file(self.out/'ldsc/LDSC_Retained_Traits.csv', ['traitname'], [['001'], ['002']])
            return result
        with patch.object(pipeline.subprocess, 'run', side_effect=run):
            self.assertEqual(pipeline.main(self.args('python', '--result_failure_action', 'drop_traits')), 0)
        _, rows = read_manifest(self.out/'manifests/gpca_traits.csv')
        self.assertEqual([row['traitname'] for row in rows], ['002', '001'])
        report = json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['removed_before_gpca'], ['003'])

    def test_missing_python_retained_manifest_stops_before_gpca(self):
        self.assertEqual(self.run_mock(self.args('python', '--result_failure_action', 'drop_traits')), 1)
        self.assertEqual([args[0] for args in self.calls], ['prepare', 'ldsc'])
        report = json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['status'], 'failed')
        self.assertIn('LDSC_Retained_Traits.csv', report['error'])

    def test_bad_manifest_duplicate_missing_and_incomplete_ref_fail(self):
        bad = [(['traitname'], [['A'],['B']]),
               (self.columns, [['A','input.vcf','NA','NA'],['A','input.vcf','NA','NA']]),
               (self.columns, [['A','missing','NA','NA'],['B','input.vcf','NA','NA']]),
               (self.columns+['ref'], [row+['no'] for row in self.rows]),
               (self.columns, [['A','input.vcf','bad','NA'],['B','input.vcf','NA','NA']])]
        for columns,rows in bad:
            with self.subTest(rows=rows), patch.object(pipeline.subprocess,'run') as run:
                csv_file(self.manifest,columns,rows)
                self.assertEqual(pipeline.main(self.args()),1)
                run.assert_not_called();self.assertFalse(self.out.exists())

    def test_binary_native_vcf_and_partial_prevalence_refused(self):
        for sample,population in [('.2','.05'), ('NA','.05'), ('.2','NA')]:
            for row in self.rows: row[2:4]=[sample,population]
            self.save_manifest()
            with patch.object(pipeline.subprocess,'run') as run:
                self.assertEqual(pipeline.main(self.args('genomicsem')),1)
                run.assert_not_called()

    def test_missing_info_and_reference_fail_before_work(self):
        args=self.args(); at=args.index('--gwama_output_info');del args[at:at+2]
        self.assertEqual(self.run_mock(args),1);self.assertFalse(self.calls)
        (self.ld/'22.l2.M_5_50').unlink()
        self.assertEqual(self.run_mock(self.args()),1);self.assertFalse(self.calls)

    def test_existing_output_is_never_overwritten(self):
        self.out.mkdir(); sentinel=self.out/'user_data';sentinel.write_text('keep')
        self.assertEqual(self.run_mock(self.args()),1)
        self.assertEqual(sentinel.read_text(),'keep');self.assertFalse(self.calls)

    def test_failed_stage_stops_following_stages_and_preserves_status(self):
        def fail(command,**kwargs):
            self.calls.append(command[3:]);return subprocess.CompletedProcess(command,7)
        with patch.object(pipeline.subprocess,'run',side_effect=fail):
            self.assertEqual(pipeline.main(self.args()),7)
        self.assertEqual(len(self.calls),1)
        report=json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['status'],'failed')
        self.assertEqual(report['stages'][0]['returncode'],7)
        self.assertEqual(report['stages'][0]['status'],'failed')
        self.assertFalse((self.out/'ldsc').exists())

    def test_success_exit_without_output_is_failure(self):
        with patch.object(pipeline.subprocess,'run',return_value=subprocess.CompletedProcess([],0)):
            self.assertEqual(pipeline.main(self.args()),1)
        report=json.loads((self.out/'Pipeline_Run_Status.json').read_text())
        self.assertEqual(report['status'],'failed')
        self.assertIn('without expected outputs',report['error'])

    def test_unknown_duplicate_or_empty_native_selection_fails(self):
        for names in (['002'],['002','002'],['002','unknown']):
            audit=self.root/'selected.csv';csv_file(audit,['traitname'],[[name] for name in names])
            with self.assertRaises(ValueError):pipeline.retained_manifest([dict(zip(self.columns,row)) for row in self.rows],audit)

    def test_numeric_strings_survive_manifest_roundtrip(self):
        self.assertEqual(self.run_mock(self.args()),0)
        _,rows=read_manifest(self.out/'manifests/gpca_traits.csv')
        self.assertEqual([row['traitname'] for row in rows],['002','001','003'])

    def test_backend_defaults_and_no_parent_parser_mutation(self):
        from ldsc_gpca.ldsc_cli import parser as original
        before=[(a.dest,a.default) for a in original._actions]
        for backend,workers in [('python',5),('genomicsem',1)]:
            parser=pipeline.build_parser(backend)
            opts=parser.parse_args(['--input','a','--outdir','b','--ld_ref','c'])
            self.assertEqual(opts.ldsc_backend,backend)
            self.assertEqual(opts.n_cores,workers)
            self.assertIsNone(opts.chisq_max)
            self.assertEqual(opts.prepare_workers,4)
            self.assertNotIn('--restart', parser._option_string_actions)
            self.assertIn('must not already exist', parser.format_help())
        self.assertEqual([(a.dest,a.default) for a in original._actions],before)
        self.assertIn('--restart', original._option_string_actions)

    def test_invalid_backend_specific_flags_and_mutually_exclusive_modes(self):
        for args in (self.args('genomicsem','--ldsc_only'),self.args('python','--munged_input'),
                     self.args('genomicsem','--munged_input','--munged_dir',str(self.munged)),
                     self.args('python','--n_cores','3'),self.args('python','--n_core','3'),
                     self.args('python','--chisq-max','80'), self.args('python','--restart'),
                     self.args('genomicsem','--rg_normalization','trait_wide'),
                     self.args('genomicsem','--chisq_max','auto')):
            with self.subTest(args=args),self.assertRaises(SystemExit):pipeline.main(args)

    def test_help_routes_do_not_run_analysis(self):
        from ldsc_gpca.cli import main
        for args in (['pipeline'],['pipeline','--help'],['pipeline','--ldsc_backend','genomicsem','--help'],
                     ['pipeline','--prepare_help'],['pipeline','--postprocess_help']):
            with patch.object(pipeline.subprocess,'run') as run:
                try: result=main(args)
                except SystemExit as error:result=error.code
                self.assertEqual(result,0);run.assert_not_called()

    def test_actual_chisq_filter_boundary_and_source_preservation(self):
        from ldsc_gpca.munging import _filter_one_munged_sumstats
        source=self.root/'source.gz'; original='SNP\tZ\nlo\t8\nat\t9\nneg\t-9\nhi\t9.01\n'
        with gzip.open(source,'wt') as stream:stream.write(original)
        filtered=self.root/'filtered.gz';excluded=self.root/'excluded.gz'
        result=_filter_one_munged_sumstats('A',str(source),str(filtered),str(excluded),81)
        self.assertEqual((result['variants_before'],result['variants_after'],result['variants_removed']),(4,3,1))
        with gzip.open(filtered,'rt') as stream:self.assertEqual(stream.read(),'SNP\tZ\nlo\t8\nat\t9\nneg\t-9\n')
        with gzip.open(source,'rt') as stream:self.assertEqual(stream.read(),original)

    def test_real_bcftools_preparation_in_native_vcf_pipeline(self):
        bcftools=shutil.which('bcftools')
        if not bcftools:self.skipTest('bcftools not installed')
        self.vcf.write_text('''##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=AF,Number=1,Type=Float,Description="ALT frequency">
##FORMAT=<ID=ES,Number=1,Type=Float,Description="ALT effect">
##FORMAT=<ID=SE,Number=1,Type=Float,Description="Standard error">
##FORMAT=<ID=LP,Number=1,Type=Float,Description="Negative log10 P">
##FORMAT=<ID=NEF,Number=1,Type=Float,Description="N">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS
1\t1000\trs1\tA\tG\t.\tPASS\t.\tAF:ES:SE:LP:NEF\t0.2:0.1:0.05:1.341986084:1000
''')
        original_run=subprocess.run
        def execute(command,**kwargs):
            if command[3]=='prepare':
                self.calls.append(command[3:])
                return original_run(command,capture_output=True,text=True,**kwargs)
            return self.fake_run(command,**kwargs)
        with patch.object(pipeline.subprocess,'run',side_effect=execute):
            self.assertEqual(pipeline.main(self.args('genomicsem','--splitby_chr','nosplit')),0)
        with (self.out/'prepare/gpca_inputs/002_GenomicPCA_inputs.tsv').open() as stream:
            rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual(rows[0]['SNPID'],'1_1000_A_G')
        self.assertEqual(float(rows[0]['Z']),2)
        self.assertEqual((rows[0]['EA'],rows[0]['OA']),('G','A'))
        with (self.out/'prepare/munge_inputs/002_munge_inputs.txt').open() as stream:
            header=stream.readline().split();row=stream.readline().split()
        values=dict(zip(header,row));self.assertEqual(values['SNP'],'rs1')
        self.assertEqual(float(values['N']),1000)


if __name__ == '__main__':
    unittest.main()
