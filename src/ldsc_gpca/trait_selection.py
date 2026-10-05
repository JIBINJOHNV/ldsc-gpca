"""Select an audited, complete LDSC subset without an R dependency.

Matches the existing R selection order: failed self-pairs, highest failed-pair
degree, lower self-h2 Z, then later manifest position. This greedy procedure
does not guarantee the largest subset. It never changes an estimate.
"""
from itertools import combinations_with_replacement
import math

import pandas as pd

from .ldsc_export import HERITABILITY_COLUMNS


def _pair_key(a, b, order):
    return tuple(sorted((a, b), key=order.__getitem__))


def validate_duplicate_estimates(frame, order, *, sources=None):
    """Use the R duplicate tolerances before exclusions can hide a conflict."""
    groups = {}
    self_columns = [c for pair in HERITABILITY_COLUMNS for c in pair if c in frame]
    self_columns += ['h2_int', 'h2_int_se']
    for row in frame.itertuples():
        groups.setdefault(_pair_key(row.p1, row.p2, order), []).append(row)
    for pair, rows in groups.items():
        if len(rows) < 2:
            continue
        conflicts = []
        columns = ['rg', 'se', 'z', 'p', 'gcov_int', 'gcov_int_se']
        # Off-diagonal h2/intercept fields belong to p2, so reverse rows
        # describe different traits. Only repeated self-pairs must agree here.
        if pair[0] == pair[1]:
            columns += self_columns
        for column in columns:
            values = [getattr(row, column) for row in rows]
            finite = [v for v in values if math.isfinite(v)]
            tolerance = .01 if column == 'z' else .001
            kinds = {'finite' if math.isfinite(v) else 'missing' if math.isnan(v)
                     else 'positive_inf' if v > 0 else 'negative_inf' for v in values}
            disagreement = len(kinds) > 1
            if len(finite) > 1:
                low, high = min(finite), max(finite)
                disagreement |= abs(high - low) > tolerance + 1e-12
            if disagreement:
                conflicts.append(column)
        if conflicts:
            locations = '' if sources is None else ' Sources: ' + '; '.join(
                f'{sources.at[row.Index, "Source_File"]}:row {sources.at[row.Index, "Source_Row"]}'
                for row in rows[:10])
            raise RuntimeError(f'Conflicting duplicate LDSC estimates: {pair[0]} <-> {pair[1]} '
                               f'[{", ".join(conflicts)}]. Trait removal cannot hide this conflict.'
                               + locations)


def _self_estimates(frame, traits):
    self_rows = {trait: rows for trait, rows in frame[frame.p1.eq(frame.p2)].groupby('p1', sort=False)}
    empty = frame.iloc[:0]
    estimates = {}
    for trait in traits:
        rows = self_rows.get(trait, empty)
        scale = 'h2_liab' if rows.h2_liab.notna().any() else 'h2_obs'
        def finite_mean(column):
            values = rows[column]
            return float(values.mean()) if len(values) and values.map(math.isfinite).all() else float('nan')
        h2, se = finite_mean(scale), finite_mean(scale + '_se')
        estimates[trait] = (h2, se, h2 / se if se > 0 else float('nan'))
    return estimates


def select_complete_traits(frame, pairs, trait_order):
    """Return retained names and exclusions, including an empty/singleton outcome.

One existing valid orientation can cover a missing opposite orientation. Any
invalid observed row makes its unordered pair unusable; contradictory observed
orientations remain fatal. Callers save audits before enforcing the minimum size.
"""
    traits = list(trait_order)
    order = {name: i for i, name in enumerate(traits)}
    validate_duplicate_estimates(frame, order)
    present, failed, reasons = set(), set(), {}
    for row in pairs.itertuples(index=False):
        pair = _pair_key(row.p1, row.p2, order)
        if row.Status != 'missing_result':
            present.add(pair)
            if row.Status != 'valid':
                failed.add(pair)
        if row.Status != 'valid':
            reasons.setdefault(pair, []).append(row.Reason)
    failed |= set(combinations_with_replacement(traits, 2)) - present
    estimates = _self_estimates(frame, traits)
    active, exclusions = set(traits), []

    def exclude(trait, reason, degree):
        h2, se, z = estimates[trait]
        exclusions.append({'Exclusion_Step': len(exclusions) + 1,
            'Manifest_Order': order[trait] + 1, 'Trait': trait, 'Reason': reason,
            'Failed_Pairs_At_Exclusion': degree, 'Self_H2': h2, 'Self_H2_SE': se, 'H2_Z': z})
        active.remove(trait)

    for trait in traits:
        if (trait, trait) in failed:
            reason = '; '.join(dict.fromkeys(reasons.get((trait, trait), ['Missing self-pair'])))
            exclude(trait, reason, sum(trait in pair for pair in failed))
    remaining = {pair for pair in failed if pair[0] != pair[1] and set(pair) <= active}
    while remaining:
        degree = {trait: sum(trait in pair for pair in remaining) for trait in active}
        def priority(trait):
            z = estimates[trait][2]
            return -degree[trait], z if math.isfinite(z) else -math.inf, -order[trait]
        selected = min(active, key=priority)
        exclude(selected, 'Removal to resolve failed/missing pair estimates '
                '(failed-pair degree, self-h2 Z, manifest-order tie-break)', degree[selected])
        remaining = {pair for pair in remaining if selected not in pair}
    columns = ('Exclusion_Step', 'Manifest_Order', 'Trait', 'Reason',
               'Failed_Pairs_At_Exclusion', 'Self_H2', 'Self_H2_SE', 'H2_Z')
    return [trait for trait in traits if trait in active], pd.DataFrame(exclusions, columns=columns)
