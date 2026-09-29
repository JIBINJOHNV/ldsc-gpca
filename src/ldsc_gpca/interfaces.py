"""Validation for the single managed CLI and manifest interface."""
import argparse
import csv
import re
import sys

# These names are rejected, never translated or accepted as alternatives.
REMOVED_OPTIONS = frozenset({
    '--input_file', '--output_folder', '--cores', '--ld', '--wld',
    '--ld_ref_snp_file', '--hapmap_file', '--ldsc_input_folder', '--munge_output',
    '--python_ldsc', '--ldsc_path', '--tolerance', '--no_mhc_exclude',
})
REMOVED_MANIFEST_COLUMNS = frozenset({
    'gwas_name', 'sampleprevalence', 'pop_prevalence', 'populationprevalence',
    'munge_inputs', 'traits',
})


def validate_option_spelling(parser, arguments):
    """Reject removed names before a managed wrapper starts any preparation."""
    for token in (sys.argv[1:] if arguments is None else arguments):
        if token == '--':
            break
        flag = token.partition('=')[0]
        if (flag in REMOVED_OPTIONS or
                (flag.startswith('--') and '-' in flag[2:]) or
                (re.fullmatch(r'-[A-Za-z]\w*', flag) and flag != '-h')):
            parser.error(f'Unrecognized option {flag}. Use the exact names in --help.')


def reject_conflicting_options(parser, arguments):
    """Check repeated scalar options before argparse can silently take the last one."""
    seen = {}
    tokens = list(sys.argv[1:] if arguments is None else arguments)
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == '--':
            break
        flag, equal, inline = token.partition('=')
        action = parser._option_string_actions.get(flag)
        index += 1
        if action is None or isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            continue
        if action.nargs == 0:
            value = action.const
        elif action.nargs is None:
            if equal:
                value = inline
            elif index < len(tokens):
                value = tokens[index]
                index += 1
            else:
                continue  # argparse reports the missing value.
            if action.type:
                try:
                    value = action.type(value)
                except (TypeError, ValueError, argparse.ArgumentTypeError):
                    continue  # argparse reports the invalid value.
        else:
            continue  # No managed option currently accepts multiple values.
        if action.dest in seen and seen[action.dest][1] != value:
            parser.error(f'Conflicting values for {seen[action.dest][0]} and {flag}')
        seen[action.dest] = (flag, value)


def read_manifest(path):
    """Read the managed CSV manifest, retaining exact names and row/column order."""
    with open(path, newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        if not columns or any(not column for column in columns) or len(columns) != len(set(columns)):
            raise ValueError('Manifest must have non-empty, unique column headers')
        rows = list(reader)
    for index, row in enumerate(rows, 2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f'Manifest CSV row {index} has a different number of fields than its header')
    removed = sorted(set(columns).intersection(REMOVED_MANIFEST_COLUMNS))
    if removed:
        raise ValueError('Unsupported manifest headers: ' + ', '.join(removed) +
                         '. Use the column names listed in --help.')
    return columns, rows
