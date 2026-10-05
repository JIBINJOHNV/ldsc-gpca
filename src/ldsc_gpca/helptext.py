"""Consistent parser help and generated R help, available without an R runtime."""
import argparse
from importlib.resources import files
import sys
import os
from .interfaces import validate_option_spelling, reject_conflicting_options
from .help_formatter import HelpFormatter


class HelpParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('formatter_class', HelpFormatter)
        kwargs.setdefault('allow_abbrev', False)
        super().__init__(*args, **kwargs)

    def parse_known_args(self, args=None, namespace=None):
        validate_option_spelling(self, args)
        reject_conflicting_options(self, args)
        return super().parse_known_args(args, namespace)

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
  --hm3: tab-separated text with SNP header; required only with
  --write_munge_inputs. Identifier selection only; no allele alignment.
  These options prepare GPCA/raw tables. For GenomicSEM LDSC munged-file reuse, run
  standalone prepare --mode ldsc (or --mode both) --munge_backend genomicsem
  with the same --hm3 allele reference for all traits. Do not directly reuse
  Python-munged files in GenomicSEM.

OUTPUTS AND FIXED INPUT CONTRACT
  GPCA output is tab-separated: SNPID,CHR,BP,EA,OA,EAF,N,Z,P (in that order).
  Optional munging output is space-separated: SNP,CHR,POS,A1,A2,eaf_A1,beta,se,N,p.
  Chromosomes 1-22 only; split mode needs all 22 chromosomes for every trait.
  Invalid rows are removed and audited; source VCF files are never changed.
  This module does not run LDSC and assumes prior reference/allele harmonisation.
"""

LDSC_INPUT_HELP = """INPUT FILE CONTRACT
  This command runs Python LDSC and writes ldsc_results.csv.
  For GenomicSEM LDSC and genomicsem_LDSC.RData, use ldsc --ldsc_backend genomicsem
  (the existing genomicsem ldsc command also remains available).
  --input: comma-separated CSV. VCF workflow headers:
    traitname,vcf_files,ref,population_prevalence,sample_prevalence
  With --ldsc_only, vcf_files is optional; all other headers remain required.
  Extra columns and any column order are allowed; required headers remain mandatory.
  Names must be unique/non-empty. Use absolute VCF paths when VCFs are processed.
  ref=yes selects an LDSC reference trait; ref=no leaves it as a target.
  Use ref=yes for every trait to generate complete GPCA pairwise coverage.
  Prevalence columns must exist; blank/NA values are allowed.
  Population prevalence absent: use NEF; sample prevalence is not used.
  Population prevalence supplied: use NC+NCO; derive missing sample prevalence.
  --hm3: whitespace-separated HapMap allele table with headers
    SNP,A1,A2 (tabs or spaces). Passed to standard LDSC --merge-alleles;
    required for VCF/munging runs and not required with --ldsc_only.
  --ld_ref: directory of chromosome reference files; see FIXED CONTRACT below.
  --munged_dir: existing .sumstats.gz input directory; requires --ldsc_only.
    Without --ldsc_only this combination is rejected before outputs are created;
    filenames: {traitname}.sumstats.gz. Required sidecars must remain alongside.
  Managed munging writes Python LDSC computed Z/N with 17 significant digits.
    --ldsc_only preserves existing precision; rerun munging to replace rounded files.
    Python-munged files are for Python LDSC. For GenomicSEM reuse, prepare
    files separately with --mode ldsc (or --mode both) --munge_backend genomicsem.
    Matching the same HapMap reference does not guarantee GenomicSEM-compatible allele orientation.
  --chisq_max INTEGER|auto: optional per-trait Z^2 <= cutoff filtering before LDSC.
    Omitted: disabled. INTEGER must be >0; auto uses max(80, 0.001 * max(N))
    separately per trait, after complete-row matching to reference AND weight
    LD-score SNPs. N is the munged file's N; auto cutoffs are not rounded.
    Each cutoff and matched maximum N are saved in LDSC_ChiSquare_Filter_Summary.csv.
    Original files are unchanged. Applies to both VCF and --ldsc_only inputs.
    Missing Z placeholders are preserved and counted separately for LDSC to discard.
    This is not forwarded to Python LDSC's cross-product --chisq-max behavior.

FIXED CONTRACT (NOT CONFIGURABLE BY THESE OPTIONS)
  Failed extraction/munging/filter jobs get one retry (2 total attempts).
  Exhausted worker failures stop the analysis, including in drop_traits mode.
  LDSC command retries default to 1 and remain configurable with --ldsc_retries.
  Local bcftools/Bash/awk extract variants; conda run launches isolated CBIIT LDSC.
  Setup: bash scripts/setup_environments.sh; Docker is not used.
  Extraction expects INFO/AF, INFO/EUR and FORMAT/SI, AF, EZ, LP, NEF;
  population-prevalence traits additionally require FORMAT/NC and FORMAT/NCO.
  Missing INFO/AF or INFO/EUR records are excluded and counted in per-trait
  munge_input/*_AF_Filter_QC.csv files; FORMAT/AF is not a fallback. Missing
  INFO header declarations are fatal. --ldsc_only reuses prior munged data.
  LD reference prefix: <ld_ref>/<CHR>.l2.ldscore.gz and associated M files;
  weights default to --ld_ref; --ld_weights can select a separate directory.
  Each LDSC batch exports .results.csv directly from its in-memory estimates.
  Compilation reads each CSV once; human-readable logs are not parsed.
  This is an EUR-field-specific extraction schema, not a generic GWAS-VCF reader.
  Threshold options below do not change field names or ancestry.
"""
