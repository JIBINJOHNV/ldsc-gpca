"""Validate and launch native GenomicSEM munging followed by two-pass LDSC."""
import csv
import math
from pathlib import Path
import shutil
import subprocess
import sys
from importlib.resources import files

from .helptext import HelpParser


def build_parser():
    p = HelpParser(prog='ldsc-gpca genomicsem ldsc', description=
        'Run native GenomicSEM munge then LDSC, or LDSC alone from supplied munged files.', epilog='''INPUT FILE CONTRACT
  --input: comma-separated CSV with unique traitname, sampleprevalence,
    populationprevalence. Both prevalence values must be blank/NA for quantitative
    traits, or both strictly between 0 and 1 for liability-scale binary traits.
    Mixed quantitative/binary rows are supported; no prevalence is inferred.
  Default mode additionally requires munge_inputs: paths to whitespace-delimited
    GWAS tables with headers SNP,A1,A2,P and a signed effect (Z or BETA), plus N.
    Optional manifest N supplies a positive constant per trait instead of file N.
    INFO and MAF/effect-allele frequency are used when recognized by GenomicSEM;
    inspect its logs for column interpretation and unavailable filtering fields.
  --hm3: reference table with SNP,A1,A2 headers (tabs/spaces), required in default mode.
  --munge-output: directory of {traitname}.sumstats.gz or .sumstats files.
    Alternatively --munged-input uses a traits column of file paths in the CSV.
    Munged files are tab-separated with SNP,A1,A2,N,Z headers.
  Relative manifest paths resolve beside the manifest. Row order is preserved.
  --ld: directory with <CHR>.l2.ldscore.gz and <CHR>.l2.M_5_50 files.
  --wld: optional separate directory with <CHR>.l2.ldscore.gz weights.

OUTPUTS / DEPENDENCIES
  Fresh output directory required; existing non-empty directories are refused.
  Writes genomicPCA_LDSC_raw.RData, genomicPCA_LDSC.RData (LDSCoutput with
  S,V,I,S_Stand,V_Stand), GenomicSEM_LDSC_Trait_QC.csv, Selected_Traits.csv,
  GenomicSEM_LDSC_Events.csv, resolved manifest, logs and sessionInfo.txt.
  Newly munged files are under outdir/munge_output. GWAMA is NOT run here.
  Requires local Rscript and GenomicSEM; no Docker or GWAMA source is used.
  stand=TRUE is required for GPCA output; matrices are never fabricated/repaired.
''')
    p.add_argument('--input', required=True, metavar='MANIFEST.csv', help='Trait manifest; columns and separator below.')
    p.add_argument('--outdir', required=True, metavar='DIRECTORY', help='Fresh output directory.')
    p.add_argument('--ld', required=True, metavar='DIRECTORY', help='LD scores and M reference files.')
    p.add_argument('--wld', metavar='DIRECTORY', help='Separate regression weights; unset uses --ld.')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--munge-output', metavar='DIRECTORY', help='Use existing per-trait munged files; skip munge.')
    group.add_argument('--munged-input', action='store_true', help='Use manifest traits paths; skip munge.')
    p.add_argument('--hm3', metavar='REFERENCE.tsv', help='HapMap reference; required unless skipping munge.')
    p.add_argument('--cores', type=int, default=1, help='Munging workers; 1 is sequential. LDSC is not parallelized by this option.')
    p.add_argument('--info-filter', type=float, default=.9, help='GenomicSEM munging INFO threshold.')
    p.add_argument('--maf-filter', type=float, default=.01, help='GenomicSEM munging MAF threshold.')
    p.add_argument('--chromosomes', type=int, default=22, help='Use chromosome files 1 through this number (1–22).')
    p.add_argument('--n-blocks', type=int, default=200, help='Requested jackknife blocks; GenomicSEM overrides this for more than 18 traits; see its log.')
    p.add_argument('--chisq-max', type=float, help='Positive chi-square exclusion threshold; unset uses GenomicSEM automatic rule.')
    p.add_argument('--invalid-h2-action', choices=['drop','error'], default='drop', help='Non-finite/non-positive raw h2: audit and drop before second pass, or stop.')
    p.add_argument('--rscript', default='Rscript', help='Rscript executable name or path.')
    return p


def resolve_manifest(opts):
    manifest = Path(opts.input).resolve()
    with manifest.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        mode = 'existing' if opts.munge_output or opts.munged_input else 'munge'
        path_column = 'traits' if opts.munged_input else 'munge_inputs'
        required = {'traitname','sampleprevalence','populationprevalence'}
        if not opts.munge_output:
            required.add(path_column)
        if len(columns) != len(set(columns)) or not required.issubset(columns):
            raise ValueError('Manifest needs unique headers including: ' + ','.join(sorted(required)))
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError('At least two traits are required.')
    names = set()
    for row in rows:
        name = row['traitname'] or ''
        if not name or any(c.isspace() for c in name) or name in names or any(c in name for c in '/\\') or name in ('.','..'):
            raise ValueError('traitname must be unique, non-empty, and safe as a filename: ' + repr(name))
        names.add(name)
        prev = []
        for column in ('sampleprevalence','populationprevalence','N'):
            value = (row.get(column) or '').strip()
            value = None if value.upper() in ('','NA','NAN','.') else float(value)
            if value is not None and (not math.isfinite(value) or value <= 0 or (column != 'N' and value >= 1)):
                raise ValueError(f'{name}: invalid {column}.')
            row[column] = 'NA' if value is None else value
            if column != 'N':
                prev.append(value)
        if (prev[0] is None) != (prev[1] is None):
            raise ValueError(f'{name}: supply both prevalences or neither; no sample prevalence is inferred.')
        if opts.munge_output:
            candidates = [Path(opts.munge_output).resolve() / (name + suffix) for suffix in ('.sumstats.gz','.sumstats')]
            candidates = [path for path in candidates if path.is_file()]
            if len(candidates) != 1:
                raise ValueError(f'{name}: need exactly one .sumstats.gz or .sumstats file in --munge-output.')
            path = candidates[0]
        else:
            value = row[path_column] or ''
            path = Path(value)
            if not path.is_absolute():
                path = manifest.parent / path
            if not value or not path.is_file():
                raise ValueError(f'{name}: input file not found: {path}')
        if path.stat().st_size == 0:
            raise ValueError(f'{name}: empty input file: {path}')
        row['source_file'] = str(path.resolve())
    return rows, mode


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    if not argv:
        parser.print_help()
        return 0
    opts = parser.parse_args(argv)
    try:
        if opts.cores < 1 or not 1 <= opts.chromosomes <= 22 or opts.n_blocks < 2:
            raise ValueError('cores must be >=1, chromosomes 1–22, n-blocks >=2.')
        if not 0 <= opts.info_filter <= 1 or not 0 <= opts.maf_filter <= .5:
            raise ValueError('info-filter must be in [0,1]; maf-filter in [0,0.5].')
        if opts.chisq_max is not None and (not math.isfinite(opts.chisq_max) or opts.chisq_max <= 0):
            raise ValueError('chisq-max must be finite and positive.')
        rows, mode = resolve_manifest(opts)
        if mode == 'munge' and (not opts.hm3 or not Path(opts.hm3).is_file()):
            raise ValueError('--hm3 reference file is required when running munge.')
        if mode == 'munge':
            with Path(opts.hm3).open() as stream:
                if not {'SNP','A1','A2'}.issubset(stream.readline().split()):
                    raise ValueError('--hm3 must have SNP,A1,A2 headers for native allele alignment.')
        if mode == 'existing' and opts.hm3:
            raise ValueError('--hm3 is unused when munging is skipped; omit it.')
        ld = Path(opts.ld).resolve()
        wld = Path(opts.wld).resolve() if opts.wld else ld
        for chrom in range(1, opts.chromosomes + 1):
            for path in (ld/f'{chrom}.l2.ldscore.gz',ld/f'{chrom}.l2.M_5_50',wld/f'{chrom}.l2.ldscore.gz'):
                if not path.is_file():
                    raise ValueError(f'Reference file not found: {path}')
        out = Path(opts.outdir).resolve()
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            raise ValueError('Use a fresh or empty --outdir; existing results are not overwritten.')
        rscript = shutil.which(opts.rscript)
        if not rscript:
            raise ValueError('Rscript not found; install R and GenomicSEM, or supply --rscript.')
    except (ValueError, OSError) as error:
        parser.error(str(error))
    out.mkdir(parents=True, exist_ok=True)
    resolved = out/'Resolved_Manifest.csv'
    with resolved.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['traitname','source_file','sampleprevalence','populationprevalence','N'])
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    script = files('ldsc_gpca').joinpath('r','genomicsem','ldsc_entry.R')
    command = [rscript, str(script), str(resolved), str(out), str(ld), str(wld),
               str(Path(opts.hm3).resolve()) if opts.hm3 else '', mode, str(opts.cores),
               str(opts.info_filter), str(opts.maf_filter), str(opts.chromosomes),
               str(opts.n_blocks), str(opts.chisq_max) if opts.chisq_max is not None else 'NA', opts.invalid_h2_action]
    return subprocess.run(command, check=False).returncode
