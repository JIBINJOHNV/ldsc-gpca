"""Annotate completed LDSC estimates; never replace native estimates or SEs."""
import argparse
import csv
import gzip
from pathlib import Path
import polars as pl
from .ldsc_export import NATIVE_COLUMNS, write_results_csv


DERIVED_COLUMNS = ('gcov_pair', 'rg_trait_wide', 'normalization_status')


def add_trait_wide_columns(frame):
    """Prefer native covariance and observed self h2; reconstruct older exports.

    Exactly one populated scale is required per row, matching that trait's
    self-pair. Mixed traits are supported without converting their scales.
    Missing/ambiguous information produces null annotations, not imputation.
    Original columns, values and row order are preserved, including duplicates.
    """
    missing = {'p1', 'p2', 'rg'} - set(frame.columns)
    if missing:
        raise ValueError(f'Missing LDSC columns: {sorted(missing)}')

    def numeric(name):
        return (pl.col(name).cast(pl.Float64, strict=False).fill_nan(None)
                if name in frame.columns else pl.lit(None, dtype=pl.Float64))

    rows = frame.select('p1', 'p2', numeric('rg').alias('_rg'),
        numeric('h2_obs').alias('_obs'), numeric('h2_liab').alias('_liab'),
        *[numeric(name).alias(alias) for name, alias in zip(NATIVE_COLUMNS, ('_ncov', '_n1', '_n2'))])
    rows = rows.with_row_index('_row').with_columns(
        pl.any_horizontal(pl.col(name).is_not_null() for name in ('_ncov', '_n1', '_n2')).alias('_native'),
        pl.coalesce('_obs', '_liab').alias('_h2'),
        pl.when(pl.col('_obs').is_not_null() & pl.col('_liab').is_null()).then(pl.lit('observed'))
        .when(pl.col('_liab').is_not_null() & pl.col('_obs').is_null()).then(pl.lit('liability'))
        .otherwise(None).alias('_scale'),
    ).with_columns(
        pl.struct('_h2', '_scale').n_unique().over('p1', 'p2').alias('_count'),
        pl.struct('_ncov', '_n1', '_n2').n_unique().over('p1', 'p2').alias('_native_count'),
        pl.when(pl.col('_native')).then(pl.col('_n2')).otherwise(pl.col('_obs')).alias('_native_h2'))
    lookup = rows.unique(subset=['p1', 'p2'], keep='first')
    rows = rows.join(lookup.select(pl.col('p2').alias('p1'), pl.col('p1').alias('p2'),
        pl.col('_h2').alias('_other_h2'), pl.col('_scale').alias('_other_scale'),
        pl.col('_count').alias('_other_count')), on=['p1', 'p2'], how='left', validate='m:1', coalesce=True)
    self_pairs = lookup.filter(pl.col('p1') == pl.col('p2'))
    for key, suffix in [('p1', 'a'), ('p2', 'b')]:
        rows = rows.join(self_pairs.select(pl.col('p1').alias(key),
            *[pl.col('_' + field).alias(f'_self_{suffix}_{field}')
              for field in ('h2', 'scale', 'count', 'native_h2', 'native_count')]),
            on=key, how='left', validate='m:1', coalesce=True)

    native = pl.col('_native')
    # Keep gcov_pair on its historical reported h2 scales. Native normalization
    # uses only observed-scale quantities, so liability conversion cannot leak in.
    rows = rows.with_columns(
        pl.when(native).then(pl.col('_n1') * pl.col('_self_a_h2') / pl.col('_self_a_native_h2'))
        .otherwise(pl.col('_other_h2')).alias('_other_h2'),
        pl.when(native).then(pl.col('_self_a_scale')).otherwise(pl.col('_other_scale')).alias('_other_scale'),
        pl.when(native).then(pl.lit(1)).otherwise(pl.col('_other_count')).alias('_other_count'))
    native_unique = pl.all_horizontal(pl.col(name) == 1 for name in
        ('_native_count', '_self_a_native_count', '_self_b_native_count'))
    native_valid = pl.col('_ncov').is_finite() & pl.all_horizontal(
        pl.col(name).is_finite() & (pl.col(name) > 0)
        for name in ('_n1', '_n2', '_self_a_native_h2', '_self_b_native_h2'))
    native_rg = pl.col('_ncov') / pl.col('_n1').sqrt() / pl.col('_n2').sqrt()
    scale_ratio = (pl.col('_h2') / pl.col('_n2')) / (pl.col('_self_b_h2') / pl.col('_self_b_native_h2'))
    consistent = ((native_rg - pl.col('_rg')).abs() <= 1e-10 + 1e-8 * pl.col('_rg').abs()) & \
                 ((scale_ratio - 1).abs() <= 1e-8)

    h2_columns = ('_h2', '_other_h2', '_self_a_h2', '_self_b_h2')
    unique = pl.all_horizontal(pl.col(name) == 1 for name in
                              ('_count', '_other_count', '_self_a_count', '_self_b_count'))
    positive = pl.all_horizontal(pl.col(name).is_finite() & (pl.col(name) > 0)
                                for name in h2_columns)
    same_scale = ((pl.col('_scale') == pl.col('_self_b_scale')) &
                  (pl.col('_other_scale') == pl.col('_self_a_scale')))
    status = (pl.when(pl.col('_other_count').is_null()).then(pl.lit('missing_reverse_pair'))
        .when(pl.col('_self_a_count').is_null() | pl.col('_self_b_count').is_null())
        .then(pl.lit('missing_self_pair'))
        .when(~unique).then(pl.lit('ambiguous_heritability'))
        .when(native & ~native_unique.fill_null(False)).then(pl.lit('ambiguous_native_estimate'))
        .when(native & ~native_valid.fill_null(False)).then(pl.lit('invalid_native_estimate'))
        .when(~positive.fill_null(False) | ~pl.col('_rg').is_finite().fill_null(False))
        .then(pl.lit('invalid_estimate'))
        .when(~same_scale.fill_null(False)).then(pl.lit('incompatible_or_ambiguous_scale'))
        .when(native & ~consistent.fill_null(False)).then(pl.lit('inconsistent_native_estimate'))
        .otherwise(pl.lit('calculated')))
    rows = rows.with_columns(status.alias('normalization_status'),
        (pl.col('_rg') * pl.col('_h2').sqrt() * pl.col('_other_h2').sqrt()).alias('gcov_pair'))
    rows = rows.with_columns(
        pl.when(pl.col('p1') == pl.col('p2')).then(1.0)
        .when(native).then(pl.col('_ncov') / pl.col('_self_a_native_h2').sqrt() / pl.col('_self_b_native_h2').sqrt())
        .otherwise(pl.col('gcov_pair') / pl.col('_self_a_h2').sqrt() / pl.col('_self_b_h2').sqrt())
        .alias('rg_trait_wide'))
    rows = rows.with_columns(pl.when(
        (pl.col('normalization_status') == 'calculated') &
        ~(pl.col('gcov_pair').is_finite() & pl.col('rg_trait_wide').is_finite()).fill_null(False))
        .then(pl.lit('non_finite_normalization')).otherwise(pl.col('normalization_status'))
        .alias('normalization_status'))
    annotations = rows.sort('_row').select(
        *[pl.when(pl.col('normalization_status') == 'calculated').then(pl.col(name))
          .otherwise(None).alias(name) for name in DERIVED_COLUMNS[:2]], 'normalization_status')
    return frame.drop([name for name in DERIVED_COLUMNS if name in frame.columns]).hstack(annotations)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Append trait-wide annotations to an existing LDSC CSV.')
    parser.add_argument('--input', required=True, help='Completed LDSC CSV (gzip accepted).')
    parser.add_argument('--out', required=True, help='Output CSV; may equal --input for atomic replacement.')
    args = parser.parse_args(argv)
    if Path(args.out).suffix == '.gz':
        parser.error('--out must be an uncompressed CSV path.')
    opener = gzip.open if str(args.input).endswith('.gz') else open
    with opener(args.input, 'rt', encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle), [])
    if not header or any(not name.strip() for name in header) or len(set(header)) != len(header):
        raise ValueError('Empty or duplicate column headers.')
    # Strings preserve original numeric text, exact trait names (including NA),
    # and unrelated metadata. Only temporary calculation columns become Float64.
    # A stream avoids Polars 0.20's glob-path handling of compressed file names.
    with opener(args.input, 'rb') as handle:
        frame = pl.read_csv(handle.read(), infer_schema_length=0, missing_utf8_is_empty_string=True)
    annotated = add_trait_wide_columns(frame)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    write_results_csv(annotated, args.out)
    print(annotated.group_by('normalization_status').len().sort('normalization_status'))


if __name__ == '__main__':
    main()
