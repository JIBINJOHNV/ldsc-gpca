"""Consistent parser help and generated R help, available without an R runtime."""
import argparse
from importlib.resources import files
import sys
import os


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _get_help_string(self, action):
        text = action.help or ''
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            return text
        if action.required:
            return text + ' [Required; no default.]'
        if action.option_strings and 'default' not in text.lower():
            text += ' [Default: %(default)s.]'
        if action.choices:
            text += ' [Choices: ' + ', '.join(map(str, action.choices)) + '.]'
        return text


class HelpParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('formatter_class', HelpFormatter)
        kwargs.setdefault('allow_abbrev', False)
        super().__init__(*args, **kwargs)

    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, f'\nERROR: {message}\n')


def print_analysis_help(r_script, arguments=(), file=None):
    backend = 'genomicsem' if r_script == 'gpsca_gwama_v2.r' else 'python_ldsc'
    text = files('ldsc_gpca').joinpath('r', backend, 'help.txt').read_text()
    file = sys.stdout if file is None else file
    mode = 'auto'
    for i, arg in enumerate(arguments):
        if arg.startswith('--color='):
            mode = arg.split('=',1)[1]
        elif arg == '--color' and i+1 < len(arguments):
            mode = arguments[i+1]
    enabled = mode == 'always' or (mode == 'auto' and 'NO_COLOR' not in os.environ and
              (file.isatty() or os.environ.get('CLICOLOR_FORCE','0') != '0'))
    if enabled:
        text = '\n'.join('\033[1;36m'+line+'\033[0m' if line and (line.isupper() or line.endswith(':')) else line
                         for line in text.split('\n'))
    print(text, end='', file=file)


PREPARE_INPUT_HELP = """INPUT FILE CONTRACT
  --input: comma-separated CSV with headers traitname,vcf_files.
  traitname: unique, non-empty identifier; vcf_files: one single-sample VCF per row.
  Relative VCF paths resolve beside the manifest. Plain VCF and .vcf.gz accepted.
  Required VCF FORMAT fields: AF, ES, SE, LP, NEF. LP means -log10(P).
  --hapmap-file: tab-separated text with SNP header; required only with
  --write-munge-inputs. Identifier selection only; no allele alignment.

OUTPUTS AND FIXED INPUT CONTRACT
  GPCA output is tab-separated: SNPID,CHR,BP,EA,OA,EAF,N,Z,P (in that order).
  Optional munging output is space-separated: SNP,CHR,POS,A1,A2,eaf_A1,beta,se,N,p.
  Chromosomes 1-22 only; split mode needs all 22 chromosomes for every trait.
  Invalid rows are removed and audited; source VCF files are never changed.
  This module does not run LDSC and assumes prior reference/allele harmonisation.
"""

LDSC_INPUT_HELP = """INPUT FILE CONTRACT
  --input_file: comma-separated CSV. VCF workflow headers:
    gwas_name,vcf_files,ref,pop_prevalence,sample_prevalence
  With --ldsc_only, vcf_files is optional; all other headers remain required.
  Names must be unique/non-empty. Use absolute VCF paths when VCFs are processed.
  ref=yes selects an LDSC reference trait; ref=no leaves it as a target.
  Use ref=yes for every trait to generate complete GPCA pairwise coverage.
  Prevalence columns must exist; blank/NA values are allowed.
  Population prevalence absent: use NEF; sample prevalence is not used.
  Population prevalence supplied: use NC+NCO; derive missing sample prevalence.
  --ld_ref_snp_file: whitespace-separated HapMap allele table with headers
    SNP,A1,A2 (tabs or spaces). Passed to standard LDSC --merge-alleles;
    required for VCF/munging runs and not required with --ldsc_only.
  --ld_ref: directory of chromosome reference files; see FIXED CONTRACT below.
  --ldsc_input_folder: .sumstats.gz directory used with --ldsc_only;
    filenames: {gwas_name}.sumstats.gz. Required sidecars must remain alongside.
  --chisq-max: optional positive threshold applied independently to each munged
    trait as Z^2 <= threshold before pairwise LDSC. Original files are unchanged.
    This is not forwarded to native LDSC's cross-product --chisq-max behavior.

FIXED CONTRACT (NOT CONFIGURABLE BY THESE OPTIONS)
  Local bcftools/Bash/awk extract variants; conda run launches isolated CBIIT LDSC.
  Setup: bash scripts/setup_environments.sh; Docker is not used.
  Extraction expects INFO/AF, INFO/EUR and FORMAT/SI, AF, EZ, LP, NEF;
  population-prevalence traits additionally require FORMAT/NC and FORMAT/NCO.
  LD reference prefix: <ld_ref>/<CHR>.l2.ldscore.gz and associated M files;
  the same directory is currently used for reference and regression weights.
  This is an EUR-field-specific extraction schema, not a generic GWAS-VCF reader.
  Threshold options below do not change field names or ancestry.
"""
