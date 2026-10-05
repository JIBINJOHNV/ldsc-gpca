"""LDSC summary-statistic munging and provenance."""
import gzip
import io
import os
import shlex
import json
import math
import shutil
import tempfile
import pandas as pd
import polars as pl
from .compression import compressed_writer, compression_settings
from .utils import DEFAULT_FILTERS, is_valid_gz, run_command
from .ldsc_runtime import ldsc_munging_command
from .ldsc_munge import FLOAT_FORMAT as MUNGE_FLOAT_FORMAT
from .workers import run_parallel_jobs
from .pairwise import _result_stamp
from .restart import file_digest


def _split_sumstats_line(line, tab_separated):
    """Split an LDSC row without changing the text written to filtered output."""
    stripped = line.rstrip('\r\n')
    return stripped.split('\t') if tab_separated else stripped.split()


def _temporary_path(destination):
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    descriptor, path = tempfile.mkstemp(
        prefix=f'.{os.path.basename(destination)}.', suffix='.tmp',
        dir=os.path.dirname(destination))
    os.close(descriptor)
    return path


def _read_ld_snps(path):
    """Project standard TSVs with Polars; retain pandas' whitespace/NA rules."""
    with gzip.open(path, 'rb') as source:
        data = source.read()
    # Only plain TSV is unambiguous under both parsers. Quotes, padding, mixed
    # whitespace and empty tab fields use the existing whitespace reader.
    if (b'\t' not in data.partition(b'\n')[0] or data.startswith(b'\t') or data.endswith(b'\t')
            or any(token in data for token in (b' ', b'"', b'\r', b'\v', b'\f', b'\t\t', b'\t\n', b'\n\t'))):
        return pd.read_csv(io.BytesIO(data), sep=r'\s+', usecols=['SNP'], dtype=str)['SNP'].dropna()
    # pandas 2.x default missing-value literals; SNP IDs stay strings.
    missing = ['', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN',
               '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA', 'NULL', 'NaN',
               'None', 'n/a', 'nan', 'null']
    return pl.read_csv(data, separator='\t', columns=['SNP'], infer_schema_length=0,
                       null_values=missing, n_threads=1)['SNP'].drop_nulls().to_list()


def _matched_ld_snps(ld_ref_dir, ld_weights_dir):
    """Load the reference/weight SNP intersection once, shared by trait workers."""
    if not ld_ref_dir:
        raise ValueError('--chisq_max auto requires an LD-score reference directory')
    matched = None
    for directory in dict.fromkeys(map(os.path.realpath, (ld_ref_dir, ld_weights_dir or ld_ref_dir))):
        snps = set()
        for chromosome in range(1, 23):
            path = os.path.join(directory, f'{chromosome}.l2.ldscore.gz')
            snps.update(_read_ld_snps(path))
        matched = snps if matched is None else matched.intersection(snps)
    if not matched:
        raise ValueError('--chisq_max auto: reference and weight LD scores have no shared SNPs')
    return frozenset(matched)


def _automatic_chisq_max(trait, source_file, matched_snps):
    """Match GenomicSEM's complete-row, per-trait maximum N before filtering."""
    maximum_n, matched_rows = 0., 0
    missing = frozenset(('', 'NA', 'nan', 'NaN', '.'))
    previous_n = None
    # A preliminary pass is needed: a later SNP can raise the cutoff for earlier SNPs.
    # Stream rows to bound memory; match readr's trim_ws before missing-value checks.
    with gzip.open(source_file, 'rt', encoding='utf-8') as source:
        header = source.readline()
        tab_separated = '\t' in header
        columns = _split_sumstats_line(header, tab_separated)
        if any(columns.count(name) != 1 for name in ('SNP', 'A1', 'A2', 'N', 'Z')):
            raise ValueError(f'{trait}: --chisq_max auto requires unique SNP,A1,A2,N,Z columns')
        snp_index, n_index = columns.index('SNP'), columns.index('N')
        for line_number, line in enumerate(source, start=2):
            values = _split_sumstats_line(line, tab_separated)
            if len(values) != len(columns):
                raise ValueError(f'{trait}: malformed munged LDSC row {line_number}')
            if values[snp_index].strip() not in matched_snps or not missing.isdisjoint(map(str.strip, values)):
                continue
            # Constant N is common. Reuse only the immediately preceding validated
            # text within this pass; no cross-file/reference/settings cache is used.
            raw_n = values[n_index]
            if raw_n != previous_n:
                try:
                    n = float(raw_n)
                except ValueError as error:
                    raise ValueError(f'{trait}: non-numeric N at LD-matched row {line_number}') from error
                if not math.isfinite(n) or n <= 0:
                    raise ValueError(f'{trait}: LD-matched N must be finite and positive (row {line_number})')
                maximum_n = max(maximum_n, n)
                previous_n = raw_n
            matched_rows += 1
    if not matched_rows:
        raise ValueError(f'{trait}: --chisq_max auto has no complete rows matching both LD-score files')
    return max(80., .001 * maximum_n), maximum_n, matched_rows


def _filter_one_munged_sumstats(trait, source_file, filtered_file, excluded_file, chisq_max,
                               matched_snps=None, compression=None):
    """Write one independently chi-square-filtered LDSC input and its exclusions."""
    automatic = chisq_max == 'auto'
    compression = compression or compression_settings(1, 9)
    maximum_n = matched_rows = None
    if automatic:
        chisq_max, maximum_n, matched_rows = _automatic_chisq_max(trait, source_file, matched_snps)
    filtered_temp = _temporary_path(filtered_file)
    excluded_temp = _temporary_path(excluded_file)
    before = kept = removed = missing_z = 0
    maximum_chisq = None
    try:
        with gzip.open(source_file, 'rt', encoding='utf-8', newline='') as source, \
                compressed_writer(filtered_temp, compression, text=True) as destination, \
                compressed_writer(excluded_temp, compression, text=True) as excluded:
            header = source.readline()
            if not header:
                raise ValueError(f'{trait}: munged LDSC input is empty')
            tab_separated = '\t' in header
            columns = _split_sumstats_line(header, tab_separated)
            if columns.count('Z') != 1 or columns.count('SNP') != 1:
                raise ValueError(f'{trait}: munged LDSC input must contain exactly one SNP and one Z column')
            z_index = columns.index('Z')
            snp_index = columns.index('SNP')
            destination.write(header if header.endswith(('\n', '\r')) else header + '\n')
            excluded.write('gwas_name\tSNP\tZ\tCHISQ\n')

            for line_number, line in enumerate(source, start=2):
                values = _split_sumstats_line(line, tab_separated)
                if len(values) != len(columns):
                    raise ValueError(
                        f'{trait}: malformed munged LDSC row {line_number}; '
                        f'expected {len(columns)} fields, observed {len(values)}')
                before += 1
                raw_z = values[z_index].strip()
                # Allele merging can leave reference SNP placeholders. Native
                # LDSC drops missing rows; preserve their text rather than impute Z.
                if raw_z in ('', 'NA', 'NaN', 'nan', '.'):
                    destination.write(line if line.endswith(('\n', '\r')) else line + '\n')
                    missing_z += 1
                    continue
                try:
                    z_value = float(raw_z)
                except ValueError as error:
                    raise ValueError(
                        f'{trait}: non-numeric Z at munged LDSC row {line_number}') from error
                if not math.isfinite(z_value):
                    raise ValueError(f'{trait}: non-finite Z at munged LDSC row {line_number}')
                chisq = z_value * z_value
                if not math.isfinite(chisq):
                    raise ValueError(f'{trait}: Z^2 overflow at munged LDSC row {line_number}')
                maximum_chisq = chisq if maximum_chisq is None else max(maximum_chisq, chisq)
                if chisq <= chisq_max:
                    destination.write(line if line.endswith(('\n', '\r')) else line + '\n')
                    kept += 1
                else:
                    excluded.write(
                        f'{trait}\t{values[snp_index]}\t{z_value:.17g}\t{chisq:.17g}\n')
                    removed += 1

        if before == 0:
            raise ValueError(f'{trait}: munged LDSC input contains no variants')
        if kept == 0:
            raise ValueError(
                f'{trait}: no variants with finite Z remain after --chisq_max {chisq_max} '
                f'({missing_z:,} missing Z; {removed:,} above threshold)')
        os.replace(excluded_temp, excluded_file)
        os.replace(filtered_temp, filtered_file)
    except Exception:
        for path in (filtered_temp, excluded_temp):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
        raise

    return {
        'gwas_name': trait,
        'chisq_max': chisq_max,
        'threshold_mode': 'auto' if automatic else 'fixed',
        'maximum_matched_n': maximum_n,
        'complete_ld_matched_rows': matched_rows,
        'variants_before': before,
        'variants_removed': removed,
        'variants_after': kept + missing_z,
        'variants_missing_z': missing_z,
        'maximum_chisq': maximum_chisq,
        'source_file': os.path.abspath(source_file),
        'filtered_file': os.path.abspath(filtered_file),
    }


def filter_munged_sumstats(n_parallel, input_df, input_folder, output_folder, chisq_max,
                          *, ld_ref_dir=None, ld_weights_dir=None):
    """Apply GenomicSEM-style per-trait Z^2 filtering to munged LDSC inputs."""
    if chisq_max != 'auto' and (isinstance(chisq_max, bool) or not isinstance(chisq_max, (int, float))
            or (isinstance(chisq_max, float) and not math.isfinite(chisq_max)) or chisq_max <= 0):
        raise ValueError('chisq_max must be a positive finite number or auto')
    if not isinstance(n_parallel, int) or isinstance(n_parallel, bool) or n_parallel < 1:
        raise ValueError('n_parallel must be a positive integer')

    filtered_folder = os.path.join(output_folder, 'ldsc_input_chisq_filtered')
    if os.path.realpath(input_folder) == os.path.realpath(filtered_folder):
        raise ValueError('ldsc_input_folder must not be the chi-square filtered output directory')
    matched_snps = _matched_ld_snps(ld_ref_dir, ld_weights_dir) if chisq_max == 'auto' else None
    # Two output streams per active trait; never give every trait the full pool.
    active = max(1, min(n_parallel, len(input_df)))
    compression = compression_settings(max(1, n_parallel // (2 * active)), 9)
    print(f"[Chi-square filtering] Compression: {compression['backend']}, level 9, "
          f"{compression['workers']} worker(s) per output stream")
    os.makedirs(filtered_folder, exist_ok=True)
    tasks = []
    for trait in input_df['gwas_name']:
        source_file = os.path.join(input_folder, f'{trait}.sumstats.gz')
        if not is_valid_gz(source_file):
            raise ValueError(f'{trait}: missing or invalid munged LDSC input {source_file}')
        tasks.append((
            trait,
            source_file,
            os.path.join(filtered_folder, f'{trait}.sumstats.gz'),
            os.path.join(filtered_folder, f'{trait}.chisq_excluded.tsv.gz'),
            chisq_max,
            matched_snps,
            compression,
        ))

    summaries = run_parallel_jobs([(task[0], _filter_one_munged_sumstats, task) for task in tasks],
        n_parallel, stage='Chi-square filtering',
        status_path=os.path.join(output_folder, 'LDSC_ChiSquare_Worker_Attempts.csv'))

    order = {trait: index for index, trait in enumerate(input_df['gwas_name'])}
    summaries.sort(key=lambda row: order[row['gwas_name']])
    summary_path = os.path.join(output_folder, 'LDSC_ChiSquare_Filter_Summary.csv')
    pd.DataFrame(summaries).to_csv(summary_path, index=False)

    excluded_path = os.path.join(output_folder, 'LDSC_ChiSquare_Excluded_Variants.tsv.gz')
    excluded_temp = _temporary_path(excluded_path)
    try:
        with compressed_writer(excluded_temp, compression_settings(n_parallel, 9), text=True) as combined:
            combined.write('gwas_name\tSNP\tZ\tCHISQ\n')
            for trait in input_df['gwas_name']:
                per_trait = os.path.join(filtered_folder, f'{trait}.chisq_excluded.tsv.gz')
                with gzip.open(per_trait, 'rt', encoding='utf-8', newline='') as source:
                    next(source)
                    shutil.copyfileobj(source, combined)
        os.replace(excluded_temp, excluded_path)
    except Exception:
        try:
            os.remove(excluded_temp)
        except FileNotFoundError:
            pass
        raise

    for row in summaries:
        print(f"{row['gwas_name']}: removed {row['variants_removed']:,} of "
              f"{row['variants_before']:,} variants with Z^2 > {row['chisq_max']}; "
              f"preserved {row['variants_missing_z']:,} missing-Z rows for LDSC to discard")
    return filtered_folder


def parallel_munge_sumstats(n_parallel, input_df, output_folder, snp_include_file, filters=None, *, ldsc_input_folder=None, munge_input_folder=None, runtime=None):
    ldsc_input_folder = ldsc_input_folder or os.path.join(output_folder, "ldsc_input")
    munge_input_folder = munge_input_folder or os.path.join(output_folder, "munge_input")
    filters = {**DEFAULT_FILTERS, **(filters or {})}
    os.makedirs(ldsc_input_folder, exist_ok=True)

    def munge_sumstats(row):
        sample_name = row['gwas_name']
        n_column = row['sample_size_column']
        in_file = os.path.join(munge_input_folder, f"{sample_name}_mungeinput.tsv")
        out_prefix = os.path.join(ldsc_input_folder, sample_name)
        final_out_file = f"{out_prefix}.sumstats.gz"

        # Re-munge full runs: an existing file may have used a different N convention.
        metadata_file = out_prefix + '.prevalence.json'
        if os.path.exists(metadata_file):
            os.remove(metadata_file)
        ignore = ['--ignore', 'N_CASES,N_CONTROLS,NEF'] if n_column == 'N_TOTAL' else []
        command = shlex.join(ldsc_munging_command(**(runtime or {})) + [
            '--sumstats', in_file, '--N-col', n_column, *ignore, '--snp', 'ID', '--a1', 'ALT',
            '--a2', 'REF', '--p', 'P', '--frq', 'AF', '--maf-min', str(filters['munge_maf_min']),
            '--signed-sumstats', 'EZ,0', '--merge-alleles', snp_include_file, '--out', out_prefix])

        previous = _result_stamp(final_out_file)
        run_command(command, f"Munge_{sample_name}", output_folder=output_folder)
        if not is_valid_gz(final_out_file) or _result_stamp(final_out_file) == previous:
            raise RuntimeError(f'{sample_name}: munging produced no fresh valid gzip output')
        digest = file_digest(final_out_file)
        with open(metadata_file, 'w') as handle:
            json.dump({'sample_size_column': n_column, 'sha256': digest, 'filters': filters,
                       'float_format': MUNGE_FLOAT_FORMAT,
                       'sample_prevalence': None if pd.isna(row['sample_prevalence']) else float(row['sample_prevalence'])}, handle)
        return sample_name

    return run_parallel_jobs([(row['gwas_name'], munge_sumstats, (row,)) for _, row in input_df.iterrows()],
        n_parallel, stage='LDSC munging',
        status_path=os.path.join(output_folder, 'LDSC_Munging_Worker_Attempts.csv'))
