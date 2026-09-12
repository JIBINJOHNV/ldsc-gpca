"""LDSC result validation, compilation and reuse checks."""
import os
import math
import pandas as pd

def compile_results(output_folder, log_files, trait_metadata):
    all_data = []
    headers = []

    print(f"  -> Parsing {len(log_files)} log files...")

    for log in log_files:
        rows_before = len(all_data)
        try:
            with open(log, 'r') as f:
                lines = f.readlines()

            for i, line in enumerate(lines):
                # Look for the header line
                if "gcov_int_se" in line:
                    # Capture headers only once
                    headers = line.split()

                    # Capture data lines until we hit a blank line or "Summary"
                    for table_line in lines[i+1:]:
                        parts = table_line.split()
                        # Ensure it's a data line (starts with a path/trait name)
                        if not parts or "Summary" in table_line:
                            break
                        # Only append if it matches the number of headers found
                        if len(parts) == len(headers):
                            all_data.append(dict(zip(headers, parts)))
                        else:
                            raise RuntimeError(f'Malformed LDSC result row in {log}: {table_line.strip()}')
                    break
        except OSError as e:
            raise RuntimeError(f'Cannot read LDSC log {log}: {e}') from e
        if len(all_data) == rows_before:
            raise RuntimeError(f'No correlation table rows found in LDSC log: {log}')

    if all_data and headers:
        df = pd.DataFrame(all_data)

        # Clean up file paths to show just trait names
        for col in ['p1', 'p2']:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: os.path.basename(x).replace('.sumstats.gz', ''))

        # Remove duplicate rows (header repeats from multiple batches)
        df = df[df['p1'] != 'p1'].drop_duplicates()
        expected = {(ref, target) for ref in trait_metadata.loc[trait_metadata['ref'] == 'yes', 'gwas_name']
                    for target in trait_metadata['gwas_name']}
        actual = set(zip(df['p1'], df['p2']))
        if expected != actual:
            raise RuntimeError(f'LDSC comparison mismatch: missing={sorted(expected - actual)}; unexpected={sorted(actual - expected)}')
        for column in ['rg']:
            values = pd.to_numeric(df[column], errors='coerce')
            bad = ~values.map(math.isfinite)
            if bad.any():
                raise RuntimeError(f'Non-finite LDSC {column} for pairs: {list(zip(df.loc[bad, "p1"], df.loc[bad, "p2"]))}')
        prevalence = trait_metadata.set_index('gwas_name')['pop_prevalence']
        # LDSC labels the entire batch h2_liab even for targets with nan prevalence.
        for name in ['h2_obs', 'h2_obs_se', 'h2_liab', 'h2_liab_se']:
            if name not in df:
                df[name] = float('nan')
        no_conversion = df['p2'].map(prevalence).isna()
        for observed, liability in [('h2_obs', 'h2_liab'), ('h2_obs_se', 'h2_liab_se')]:
            df.loc[no_conversion, observed] = df.loc[no_conversion, observed].fillna(df.loc[no_conversion, liability])
            df.loc[no_conversion, liability] = float('nan')
        df['h2_scale'] = df['p2'].map(lambda name: 'NEF_unconverted' if pd.isna(prevalence[name]) else 'liability')

        output_path = os.path.join(output_folder, 'ldsc_results.csv')
        df.to_csv(output_path, index=False)
        print(f"  -> Successfully compiled {len(df)} correlations.")
    else:
        raise RuntimeError('No correlation results found in current LDSC logs')


def check_saved_filters(metadata, filters, trait):
    """Reject changed filters; legacy files cannot establish filter provenance."""
    saved = metadata.get('filters')
    if saved is None:
        print(f'WARNING: {trait}: previous filter settings are unknown; rerun without --ldsc_only for verified filter provenance.')
    elif saved != filters:
        raise ValueError(f'{trait}: filter settings differ from the saved run; rerun without --ldsc_only')
