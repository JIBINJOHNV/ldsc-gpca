"""Report LDSC estimation failures without changing any reported estimate."""
import math
import pandas as pd

from .ldsc_export import BASE_RESULT_COLUMNS, HERITABILITY_COLUMNS


def _missing_value(series):
    return series.isna() | series.astype(str).str.strip().isin(('', 'NA', 'NaN', 'nan'))


def result_status(frame, source):
    """Classify rows; missing estimates are reportable, malformed text is fatal."""
    h2_pairs = [pair for pair in HERITABILITY_COLUMNS if set(pair) <= set(frame)]
    active_pairs = {pair: ~_missing_value(frame[pair[0]]) | ~_missing_value(frame[pair[1]])
                    for pair in h2_pairs}
    no_scale = ~pd.DataFrame(active_pairs).any(axis=1)
    numeric = [c for c in BASE_RESULT_COLUMNS if c not in ('p1', 'p2')]
    numeric += [c for pair in h2_pairs for c in pair]
    reasons = [[] for _ in range(len(frame))]
    warnings = [[] for _ in range(len(frame))]

    def mark(mask, message, destination=reasons):
        for i, applies in enumerate(mask):
            if applies:
                destination[i].append(message)

    for column in numeric:
        original = frame[column]
        missing = _missing_value(original)
        try:
            converted = original.where(~missing, float('nan')).map(float)
        except (ValueError, TypeError) as error:
            raise RuntimeError(f'Non-numeric LDSC {column} in {source}') from error
        frame[column] = converted
        # An empty alternative scale for this row is allowed. If neither scale
        # has an estimate, flag the row rather than treating both as unused.
        active = next((mask | no_scale for pair, mask in active_pairs.items() if column in pair),
                      pd.Series(True, index=frame.index))
        mark(active & ~converted.map(math.isfinite), f'Non-finite LDSC {column}')
        if column.endswith('_se') or column == 'se':
            mark(converted.le(0), f'Non-positive LDSC {column}')
        if column in ('h2_obs', 'h2_liab'):
            mark(converted.le(0), f'Non-positive LDSC {column}')
    mark(frame.p.notna() & ~frame.p.between(0, 1), 'LDSC p-values outside [0,1]')
    self_pair = frame.p1 == frame.p2
    mark(self_pair & frame.rg.map(math.isfinite) & (frame.rg - 1).abs().gt(0.01 + 1e-12),
         'Self-pair rg differs from 1 beyond tolerance 0.01')
    mark(~self_pair & frame.rg.map(math.isfinite) & frame.rg.abs().gt(1),
         'Estimated rg outside [-1,1]', warnings)
    for h2, se in h2_pairs:
        mark(self_pair & frame[h2].gt(0) & frame[se].gt(0) & (frame[h2] / frame[se]).lt(2),
             f'Low self {h2}/SE (<2); diagnostic only', warnings)
    return pd.DataFrame({
        'p1': frame.p1, 'p2': frame.p2,
        'Source_File': str(source), 'Source_Row': range(2, len(frame) + 2),
        'Status': ['failed_estimate' if r else 'valid' for r in reasons],
        'Reason': ['; '.join(r) for r in reasons],
        'Warning': ['; '.join(w) for w in warnings],
    })


def trait_status(frame, pairs, trait_order):
    """Manifest-ordered self estimates and counts; pair h2 belongs to p2."""
    records = []
    for position, trait in enumerate(trait_order, 1):
        self_rows = frame[(frame.p1 == trait) & (frame.p2 == trait)]
        self_qc = pairs[(pairs.p1 == trait) & (pairs.p2 == trait)]
        involved = pairs[(pairs.p1 == trait) | (pairs.p2 == trait)]
        record = {'Manifest_Order': position, 'Trait': trait,
                  'Self_Status': ('not_requested' if self_qc.empty else
                                  'failed_estimate' if self_qc.Status.ne('valid').any() else 'valid'),
                  'Self_Reason': '; '.join(dict.fromkeys(self_qc.Reason[self_qc.Reason.ne('')])),
                  'Failed_Pair_Rows': int(involved.Status.ne('valid').sum())}
        for col in ('h2_obs', 'h2_obs_se', 'h2_liab', 'h2_liab_se', 'h2_int', 'h2_int_se'):
            record[col] = self_rows[col].iloc[0] if len(self_rows) else float('nan')
        records.append(record)
    return pd.DataFrame(records)
