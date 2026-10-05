"""Validate and launch GenomicSEM munging followed by two-pass LDSC."""
import csv
import math
from pathlib import Path
import shutil
import subprocess
import sys
from importlib.resources import files

from .helptext import HelpParser
from .interfaces import read_manifest
from .genomicsem_inputs import validate_reuse, write_input_audit, InputCompatibilityError, AUDIT_FILE


def build_parser():
    p = HelpParser(prog='ldsc-gpca genomicsem ldsc', description=
        'Run GenomicSEM LDSC from raw GWAS tables, GWAS VCFs or munged files.',
        usage='%(prog)s --input MANIFEST.csv --outdir DIRECTORY --ld_ref DIRECTORY [options]', epilog='''INPUT FILE CONTRACT
  --input: comma-separated CSV with unique traitname, sample_prevalence,
    population_prevalence. Both prevalence values must be blank/NA for quantitative
    traits, or both strictly between 0 and 1 for liability-scale binary traits.
    Mixed quantitative/binary rows are supported; no prevalence is inferred.
  Default mode additionally requires sumstats_file: paths to whitespace-delimited
    GWAS tables with headers SNP,A1,A2,P and a signed effect (Z or BETA), plus N.
    Optional manifest N supplies a positive constant per trait instead of file N.
    INFO and MAF/effect-allele frequency are used when recognized by GenomicSEM;
    inspect its logs for column interpretation and unavailable filtering fields.
  --hm3: SNP,A1,A2 reference table (tabs/spaces), required for raw and VCF modes.
  --vcf_input: read manifest vcf_files, convert GWAS VCFs to raw tables, then munge.
    Each VCF needs exactly one GWAS sample, variant IDs and FORMAT AF,ES,SE,LP,SI.
    Zero-sample and multi-sample VCFs are rejected before extraction.
    Quantitative N comes from FORMAT/NEF unless manifest N overrides it.
    Binary VCFs require an explicit appropriate manifest N and both prevalences;
    no case-count or effective-N convention is inferred. VCF IDs must match hm3.
    GenomicSEM INFO/MAF filters apply during munging; Python extraction filters do not.
  --munged_dir: directory of {traitname}.sumstats.gz or .sumstats files.
    Alternatively --munged_input uses a munged_file column of file paths in the CSV.
    Munged files are tab-separated with SNP,A1,A2,N,Z headers.
    Prepare reused files with prepare --mode ldsc (or --mode both) and
    --munge_backend genomicsem, using the same --hm3 allele reference for all traits.
    Do not directly reuse Python-munged files: HapMap matching can retain strand
    complements that GenomicSEM LDSC does not align correctly. Keep the entire
    completed preparation bundle: manifest, settings, sidecars and allele reference.
    Reuse verifies backend/runtime provenance, checksums, N/prevalences, unique SNPs
    and exact reference A1/A2 order. Swapped/complemented/incompatible pairs fail;
    no alleles or Z values are changed. Legacy files without this provenance fail.
    GenomicSEM_Input_QC.csv records acceptance or the first failing trait/SNP.
    Raw-table and --vcf_input modes run GenomicSEM munging themselves.
  Relative manifest paths resolve beside the manifest. Row order is preserved.
  --ld_ref: directory with <CHR>.l2.ldscore.gz and <CHR>.l2.M_5_50 files.
  --ld_weights: optional separate directory with <CHR>.l2.ldscore.gz weights.

OUTPUTS / DEPENDENCIES
  Fresh output directory required; existing non-empty directories are refused.
  Writes genomicsem_LDSC_raw.RData, genomicsem_LDSC.RData (LDSCoutput with
  S,V,I,S_Stand,V_Stand), GenomicSEM_LDSC_Trait_QC.csv, Selected_Traits.csv,
  GenomicSEM_LDSC_Events.csv, resolved manifest, logs and sessionInfo.txt.
  Newly munged files are under outdir/munge_output. GWAMA is NOT run here.
  VCF mode also writes vcf_input/*_munge_inputs.tsv and GenomicSEM_VCF_* QC/attempts.
  Requires local Rscript and GenomicSEM; no Docker or GWAMA source is used.
  stand=TRUE is required for GPCA output; matrices are never fabricated/repaired.
''')
    required = p.add_argument_group('Required inputs')
    inputs = p.add_argument_group('Input mode and references')
    filters = p.add_argument_group('Munging and chi-square filters')
    regression = p.add_argument_group('LDSC settings and failed heritabilities')
    execution = p.add_argument_group('Execution and local tools')
    required.add_argument('--input', required=True, metavar='MANIFEST.csv', help='Trait manifest; columns and separator below.')
    required.add_argument('--outdir', required=True, metavar='DIRECTORY', help='Fresh output directory.')
    required.add_argument('--ld_ref', required=True, metavar='DIRECTORY', help='LD scores and M reference files.')
    inputs.add_argument('--ld_weights', metavar='DIRECTORY', help='Separate regression weights. Default: use --ld_ref.')
    group = inputs.add_mutually_exclusive_group()
    group.add_argument('--munged_dir', metavar='DIRECTORY', help='Reuse files from prepare --munge_backend genomicsem; skip munge. Default: unset. Mutually exclusive with --munged_input and --vcf_input.')
    group.add_argument('--munged_input', action='store_true', help='Reuse manifest munged_file paths from prepare --munge_backend genomicsem; skip munge.')
    group.add_argument('--vcf_input', action='store_true', help='Read manifest vcf_files; extract raw tables with INFO, then run GenomicSEM munging and LDSC. Default: off; read raw sumstats_file unless using a munged mode.')
    inputs.add_argument('--hm3', metavar='REFERENCE.tsv', help='HapMap reference. Required for raw/VCF munging; omit for reuse, which verifies the bundled reference.')
    execution.add_argument('--n_cores', dest='n_cores', type=int, default=1,
                   help='VCF extraction and munging workers; 1 is sequential. Each failed trait gets one retry (2 total attempts); exhausted failures stop before LDSC. LDSC regression stays sequential.')
    vcf = p.add_argument_group('VCF conversion (--vcf_input only)')
    vcf.add_argument('--bcftools', default='bcftools', help='VCF query executable name/path. Default: bcftools on PATH; unused for raw or munged tables.')
    vcf.add_argument('--p_min', type=float, default=1e-300, help='Floor for P calculated from FORMAT/LP; range (0,1). Default: 1e-300. Adjusted records are audited; GenomicSEM munging derives Z from P and effect direction.')
    filters.add_argument('--info_filter', type=float, default=.9, help='GenomicSEM munging INFO threshold.')
    filters.add_argument('--maf_filter', type=float, default=.01, help='GenomicSEM munging MAF threshold.')
    regression.add_argument('--chromosomes', type=int, default=22, help='Use chromosome files 1 through this number (1–22).')
    regression.add_argument('--n_blocks', type=int, default=200, help='Requested jackknife blocks; GenomicSEM overrides this for more than 18 traits; see its log.')
    filters.add_argument('--chisq_max', type=float, help='Positive chi-square exclusion threshold. Default: unset; use GenomicSEM automatic rule.')
    regression.add_argument('--invalid_h2_action', choices=['drop','error'], default='drop', help='Non-finite/non-positive raw h2: audit and drop before second pass, or stop.')
    execution.add_argument('--rscript', default='Rscript', help='Rscript executable name or path.')
    return p


def resolve_manifest(opts):
    manifest = Path(opts.input).resolve()
    columns, rows = read_manifest(manifest)
    mode = 'existing' if opts.munged_dir or opts.munged_input else 'vcf' if opts.vcf_input else 'munge'
    path_column = 'vcf_files' if opts.vcf_input else 'munged_file' if opts.munged_input else 'sumstats_file'
    required = {'traitname', 'sample_prevalence', 'population_prevalence'}
    if not opts.munged_dir:
        required.add(path_column)
    if not required.issubset(columns):
        raise ValueError('Manifest needs headers including: ' + ','.join(sorted(required)))
    # Keep the GenomicSEM R runner and existing audit output schema compatible.
    for row in rows:
        row['sampleprevalence'] = row.pop('sample_prevalence')
        row['populationprevalence'] = row.pop('population_prevalence')
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
        if mode == 'vcf' and prev[1] is not None and row['N'] == 'NA':
            raise ValueError(f'{name}: binary VCF input requires an explicit appropriate N in the manifest; '
                             'do not assume FORMAT/NEF or NC+NCO has the required convention.')
        if opts.munged_dir:
            candidates = [Path(opts.munged_dir).resolve() / (name + suffix) for suffix in ('.sumstats.gz','.sumstats')]
            candidates = [path for path in candidates if path.is_file()]
            if len(candidates) != 1:
                raise ValueError(f'{name}: need exactly one .sumstats.gz or .sumstats file in --munged_dir.')
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
    if mode == 'existing':
        validate_reuse(rows)
    return rows, mode


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    if not argv:
        parser.print_help()
        return 0
    opts = parser.parse_args(argv)
    try:
        out = Path(opts.outdir).resolve()
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            raise ValueError('Use a fresh or empty --outdir; existing results are not overwritten.')
        if opts.n_cores < 1 or not 1 <= opts.chromosomes <= 22 or opts.n_blocks < 2:
            raise ValueError('--n_cores must be >=1, chromosomes 1–22, n-blocks >=2.')
        if not 0 <= opts.info_filter <= 1 or not 0 <= opts.maf_filter <= .5:
            raise ValueError('info-filter must be in [0,1]; maf-filter in [0,0.5].')
        if opts.chisq_max is not None and (not math.isfinite(opts.chisq_max) or opts.chisq_max <= 0):
            raise ValueError('chisq-max must be finite and positive.')
        if not math.isfinite(opts.p_min) or not 0 < opts.p_min < 1:
            raise ValueError('--p_min must be finite and strictly between 0 and 1.')
        rows, mode = resolve_manifest(opts)
        if mode != 'existing' and (not opts.hm3 or not Path(opts.hm3).is_file()):
            raise ValueError('--hm3 reference file is required when running munge.')
        if mode != 'existing':
            with Path(opts.hm3).open() as stream:
                if not {'SNP','A1','A2'}.issubset(stream.readline().split()):
                    raise ValueError('--hm3 must have SNP,A1,A2 headers for GenomicSEM allele alignment.')
        if mode == 'existing' and opts.hm3:
            raise ValueError('--hm3 is unused when munging is skipped; omit it.')
        ld = Path(opts.ld_ref).resolve()
        wld = Path(opts.ld_weights).resolve() if opts.ld_weights else ld
        for chrom in range(1, opts.chromosomes + 1):
            for path in (ld/f'{chrom}.l2.ldscore.gz',ld/f'{chrom}.l2.M_5_50',wld/f'{chrom}.l2.ldscore.gz'):
                if not path.is_file():
                    raise ValueError(f'Reference file not found: {path}')
        rscript = shutil.which(opts.rscript)
        if not rscript:
            raise ValueError('Rscript not found; install R and GenomicSEM, or supply --rscript.')
        if mode == 'vcf':
            executable = shutil.which(opts.bcftools)
            if not executable:
                raise ValueError('VCF input requires bcftools; install it or supply --bcftools.')
    except (ValueError, OSError) as error:
        if isinstance(error, InputCompatibilityError):
            out.mkdir(parents=True, exist_ok=True)
            write_input_audit(out/AUDIT_FILE, [error.audit])
        parser.error(str(error))
    out.mkdir(parents=True, exist_ok=True)
    if mode == 'existing':
        write_input_audit(out/AUDIT_FILE, [row['_reuse_audit'] for row in rows])
    if mode == 'vcf':
        from .genomicsem_vcf import prepare_vcf_inputs
        try:
            rows = prepare_vcf_inputs(rows, out, executable, opts.n_cores, opts.p_min)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        mode = 'munge'
    resolved = out/'Resolved_Manifest.csv'
    with resolved.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['traitname','source_file','sampleprevalence','populationprevalence','N'])
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    script = files('ldsc_gpca').joinpath('r','genomicsem','ldsc_entry.R')
    command = [rscript, str(script), str(resolved), str(out), str(ld), str(wld),
               str(Path(opts.hm3).resolve()) if opts.hm3 else '', mode, str(opts.n_cores),
               str(opts.info_filter), str(opts.maf_filter), str(opts.chromosomes),
               str(opts.n_blocks), str(opts.chisq_max) if opts.chisq_max is not None else 'NA', opts.invalid_h2_action]
    return subprocess.run(command, check=False).returncode
