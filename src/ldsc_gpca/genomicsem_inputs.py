"""Strict, read-only validation of prepared GenomicSEM regression inputs.

The strand rule is exact reference A1/A2 order. No input row is dropped, flipped
or rewritten. Checksums describe provenance; the SNP scan verifies compatibility.
"""
import csv
import gzip
import hashlib
import math
from pathlib import Path

from .restart import file_digest
from .preparation_provenance import read_json_object, validate_completed_bundle

GENOMICSEM_COMMIT = '6b65ca5db39fdade08b0d811477be1cdd57b5039'
REFERENCE_FILE = 'Allele_Reference.tsv'
AUDIT_FILE = 'GenomicSEM_Input_QC.csv'
AUDIT_COLUMNS = ['traitname', 'source_file', 'status', 'reason', 'row', 'SNP',
                 'A1', 'A2', 'reference_A1', 'reference_A2', 'n_snps',
                 'snp_order_sha256', 'sumstats_sha256', 'reference_sha256']


class InputCompatibilityError(ValueError):
    def __init__(self, reason, audit):
        self.audit = {**audit, 'status': 'error', 'reason': reason}
        super().__init__(f"{audit.get('traitname', 'reference')}: {reason}. "
                         'Re-run prepare --mode ldsc --munge_backend genomicsem from the raw inputs; '
                         'do not relabel or edit munged files.')


def write_input_audit(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=AUDIT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def _records(path, columns, *, tabs=False):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        header = stream.readline().rstrip('\r\n')
        separated = tabs or '\t' in header
        fields = next(csv.reader([header], delimiter='\t')) if separated else header.split()
        if any(fields.count(column) != 1 for column in columns):
            raise ValueError(f'{path}: expected unique {",".join(columns)} headers')
        for number, line in enumerate(stream, 2):
            values = next(csv.reader([line], delimiter='\t')) if separated else line.split()
            if len(values) != len(fields):
                raise ValueError(f'{path}: malformed row {number}')
            yield number, dict(zip(fields, values))


def read_allele_reference(path):
    reference = {}
    for number, row in _records(path, ('SNP', 'A1', 'A2')):
        snp, pair = row['SNP'], (row['A1'], row['A2'])
        if not snp or snp != snp.strip() or snp in reference:
            raise ValueError(f'{path}: empty, padded or duplicate reference SNP at row {number}: {snp!r}')
        if pair[0] not in 'ACGT' or len(pair[0]) != 1 or pair[1] not in 'ACGT' or len(pair[1]) != 1 or pair[0] == pair[1]:
            raise ValueError(f'{path}: invalid reference allele pair at row {number}: {pair}')
        reference[snp] = pair
    if not reference:
        raise ValueError(f'{path}: empty allele reference')
    return reference


def validate_alleles(path, trait, reference, *, n_override=None):
    """Scan every row, retain exact signed Z, and report the first invalid SNP."""
    audit = {'traitname': trait, 'source_file': str(path)}
    seen, snps = set(), hashlib.sha256()
    complement = str.maketrans('ACGT', 'TGCA')
    try:
        for number, row in _records(path, ('SNP', 'N', 'Z', 'A1', 'A2'), tabs=True):
            snp, pair = row['SNP'], (row['A1'], row['A2'])
            expected = reference.get(snp)
            location = {**audit, 'row': number, 'SNP': snp, 'A1': pair[0], 'A2': pair[1],
                        'reference_A1': expected[0] if expected else '',
                        'reference_A2': expected[1] if expected else ''}
            reason = None
            if snp in seen:
                reason = 'duplicate SNP ID'
            elif expected is None:
                reason = 'SNP absent from allele reference'
            elif pair != expected:
                if pair == expected[::-1]:
                    reason = 'swapped allele order'
                elif tuple(a.translate(complement) for a in pair) == expected:
                    reason = 'strand-complemented allele pair'
                elif tuple(a.translate(complement) for a in pair) == expected[::-1]:
                    reason = 'swapped strand-complemented allele pair'
                else:
                    reason = 'incompatible allele pair'
            if reason:
                raise InputCompatibilityError(f'{reason} at row {number}, SNP {snp!r}', location)
            try:
                n, z = float(row['N']), float(row['Z'])
            except ValueError:
                raise InputCompatibilityError(f'nonnumeric N/Z at row {number}, SNP {snp!r}', location) from None
            if not math.isfinite(n) or n <= 0 or not math.isfinite(z):
                raise InputCompatibilityError(f'invalid N/Z at row {number}, SNP {snp!r}', location)
            if n_override is not None and not math.isclose(n, n_override, rel_tol=1e-12, abs_tol=0):
                raise InputCompatibilityError(f'N differs from preparation override at row {number}, SNP {snp!r}', location)
            seen.add(snp)
            snps.update((snp+'\n').encode())
        if not seen:
            raise ValueError('no usable SNP rows')
    except InputCompatibilityError:
        raise
    except (ValueError, OSError, EOFError) as error:
        raise InputCompatibilityError(str(error), audit) from error
    return {**audit, 'status': 'accepted', 'reason': 'exact reference allele order; no edits',
            'n_snps': len(seen), 'snp_order_sha256': snps.hexdigest()}


def _check_metadata(metadata, row):
    if metadata.get('schema_version') != 1 or metadata.get('preparation_method') != 'shared':
        raise ValueError('missing supported preparation provenance')
    if metadata.get('munge_backend') != 'genomicsem':
        raise ValueError('munging backend must be genomicsem')
    if metadata.get('traitname') != row['traitname']:
        raise ValueError('sidecar traitname differs from manifest')
    for field in ('sample_prevalence', 'population_prevalence'):
        expected = row[field.replace('_', '')]
        expected = None if expected == 'NA' else expected
        if field not in metadata or metadata[field] != expected:
            raise ValueError(f'{field} differs from preparation metadata')
    if 'N_override' not in metadata:
        raise ValueError('missing preparation N_override')
    override = metadata['N_override']
    if override is not None and (isinstance(override, bool) or not isinstance(override, (int, float)) or not math.isfinite(override) or override <= 0):
        raise ValueError('invalid preparation N_override')
    binary = metadata['population_prevalence'] is not None
    convention = 'total' if binary else 'manifest_N' if override is not None else 'NEF'
    if (metadata.get('N_convention') != convention or
            metadata.get('sample_size_column') != ('N_TOTAL' if binary else 'NEF') or
            metadata.get('N_source') != ('manifest_N' if override is not None else 'FORMAT/NEF') or
            (binary and override is None)):
        raise ValueError('inconsistent preparation N convention/source')
    if row['N'] != 'NA' and row['N'] != override:
        raise ValueError('manifest N differs from preparation N_override; reused N cannot be replaced')
    runtime = metadata.get('runtime')
    if (not isinstance(runtime, dict) or runtime.get('source_commit') != GENOMICSEM_COMMIT or
            not runtime.get('package_version') or not runtime.get('r_version')):
        raise ValueError('missing pinned GenomicSEM runtime provenance')
    return override


def validate_reuse(rows):
    """Verify selected files against completed bundles and one common reference."""
    references, bundles, common_reference = {}, {}, None
    for row in rows:
        path = Path(row['source_file'])
        audit = {'traitname': row['traitname'], 'source_file': str(path)}
        try:
            sidecar = path.parent/f"{row['traitname']}.prevalence.json"
            metadata = read_json_object(sidecar)
            override = _check_metadata(metadata, row)
            actual, artifacts = validate_completed_bundle(path, sidecar, metadata, cache=bundles)
            if metadata.get('reference_file') != REFERENCE_FILE:
                raise ValueError('missing bundled allele reference')
            reference_path = path.parent/REFERENCE_FILE
            if reference_path not in references:
                digest = file_digest(reference_path)
                if metadata.get('reference_sha256') != digest or artifacts.get(REFERENCE_FILE) != digest:
                    raise ValueError('allele reference checksum mismatch')
                references[reference_path] = (digest, read_allele_reference(reference_path))
            digest, reference = references[reference_path]
            if metadata.get('reference_sha256') != digest or artifacts.get(REFERENCE_FILE) != digest:
                raise ValueError('allele reference checksum mismatch')
            if common_reference is not None and common_reference != digest:
                raise ValueError('traits were prepared against different allele references')
            common_reference = digest
            audit = validate_alleles(path, row['traitname'], reference, n_override=override)
            audit.update(sumstats_sha256=actual, reference_sha256=digest)
            row['_reuse_audit'] = audit
        except InputCompatibilityError:
            raise
        except (ValueError, OSError, EOFError, RuntimeError) as error:
            raise InputCompatibilityError(str(error), audit) from error
