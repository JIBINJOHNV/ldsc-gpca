"""Optional real LDSC checks; set LDSC_GPCA_TEST_LDSC_PREFIX to a usable runtime."""
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from ldsc_gpca import pairwise, results


PREFIX = os.environ.get("LDSC_GPCA_TEST_LDSC_PREFIX")
CONDA = os.environ.get("LDSC_GPCA_TEST_CONDA", "conda")


def make_fixture(root):
    """Simulated quantitative traits and LD scores, not real biological results."""
    source, reference = root / "sumstats", root / "reference"
    source.mkdir()
    reference.mkdir()
    rng = np.random.RandomState(294)
    per_chromosome = 1000
    size = 22 * per_chromosome
    ids = [f"rs{i + 1}" for i in range(size)]
    scores = rng.uniform(1, 50, size)
    total_m = 22 * 22727
    sample_size = 50000
    shared_genetic = rng.normal(size=size)
    shared_error = rng.normal(size=size)
    for name, heritability in zip(("A", "B", "C"), (0.2, 0.25, 0.15)):
        genetic = np.sqrt(0.6) * shared_genetic + np.sqrt(0.4) * rng.normal(size=size)
        error = np.sqrt(0.15) * shared_error + np.sqrt(0.85) * rng.normal(size=size)
        z = np.sqrt(sample_size * heritability * scores / total_m) * genetic + error
        pd.DataFrame({"SNP": ids, "A1": "A", "A2": "G", "N": sample_size, "Z": z}).to_csv(
            source / f"{name}.sumstats.gz", sep="\t", index=False, compression="gzip"
        )
    for chromosome in range(1, 23):
        start = (chromosome - 1) * per_chromosome
        stop = start + per_chromosome
        pd.DataFrame({
            "CHR": chromosome, "SNP": ids[start:stop],
            "BP": np.arange(1, per_chromosome + 1) * 1000, "L2": scores[start:stop],
        }).to_csv(reference / f"{chromosome}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        for suffix in (".l2.M", ".l2.M_5_50"):
            (reference / f"{chromosome}{suffix}").write_text("22727\n")
    metadata = pd.DataFrame({
        "gwas_name": ["A", "B", "C"], "ref": "yes",
        "sample_prevalence": float("nan"), "pop_prevalence": float("nan"),
    })
    return source, reference, metadata


@unittest.skipUnless(PREFIX and shutil.which(CONDA), "requires an explicit usable LDSC runtime")
class RealLdscTests(unittest.TestCase):
    def test_native_table_exports_exact_values_and_unchanged_rendering(self):
        # Exercise the installed native table builder, not a mocked implementation.
        script = """
import sys
import runpy
from types import SimpleNamespace
import pandas as pd
import ldscore.sumstats as native
export_rg_tables = runpy.run_path(sys.argv[2])['export_rg_tables']
pd.set_option('display.float_format', '{:.4f}'.format)
args = SimpleNamespace(out=sys.argv[1], samp_prev=None, pop_prev=None)
hsq = SimpleNamespace(tot=0.06123456789123456, tot_se=0.01234567891234567,
                      intercept=1.0123456789123456, intercept_se=1.234567891234567e-08)
cov = SimpleNamespace(intercept=0.21234567891234567, intercept_se=2.345678912345678e-09)
estimate = SimpleNamespace(rg_ratio=0.9999999999999876, rg_se=1.9647123456789123e-06,
                           z=508990.39251234567, p=3.212345678912345e-100, hsq2=hsq, gencov=cov)
expected = native._get_rg_table(['A', 'A'], [estimate], args)
with export_rg_tables(native, pd) as paths:
    observed = native._get_rg_table(['A', 'A'], [estimate], args)
assert observed == expected
assert '0.0000' in observed
data = pd.read_csv(paths[0], float_precision='round_trip')
for column, value in {'rg': estimate.rg_ratio, 'se': estimate.rg_se, 'z': estimate.z,
                      'p': estimate.p, 'h2_obs': hsq.tot, 'h2_obs_se': hsq.tot_se,
                      'h2_int': hsq.intercept, 'h2_int_se': hsq.intercept_se,
                      'gcov_int': cov.intercept, 'gcov_int_se': cov.intercept_se}.items():
    assert data[column].iloc[0] == value, column
print('Exact native estimates preserved; readable rendering unchanged.')
"""
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [CONDA, "run", "--no-capture-output", "--prefix", PREFIX, "python", "-B", "-c",
                 script, str(Path(directory) / "batch"),
                 str(Path(__file__).resolve().parents[1] / "src" / "ldsc_gpca" / "ldsc_export.py")],
                capture_output=True, text=True,
            )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_three_trait_regressions_compile_without_reading_logs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, reference, metadata = make_fixture(root)
            exported = pairwise.parallel_ldsc_analysis(
                3, 100, str(root / "batches"), str(reference) + os.sep, metadata, str(source),
                retries=0, runtime={"conda": CONDA, "prefix": PREFIX},
            )
            self.assertEqual(len(exported), 3)
            self.assertTrue(all(path.endswith(".results.csv") for path in exported))
            self.assertEqual(len(list((root / "batches").glob("*.log"))), 3)
            with patch.object(results, "_read_legacy_log", side_effect=AssertionError("log read")), \
                    patch.object(results.pd, "read_csv", wraps=pd.read_csv) as read:
                results.compile_results(str(root), exported, metadata)
                self.assertEqual(read.call_count, 3)
            compiled = pd.read_csv(root / "ldsc_results.csv", float_precision="round_trip")
            self.assertEqual(len(compiled), 9)
            self.assertEqual(set(zip(compiled.p1, compiled.p2)),
                             {(a, b) for a in ("A", "B", "C") for b in ("A", "B", "C")})
            self_pairs = compiled[compiled.p1 == compiled.p2]
            self.assertTrue(self_pairs.se.between(0, 5e-5, inclusive="neither").all())
            self.assertTrue((self_pairs.rg - 1).abs().lt(0.01).all())
            self.assertTrue(compiled.h2_obs.gt(0).all())
            self.assertTrue(compiled[["se", "h2_obs_se", "h2_int_se", "gcov_int_se"]].gt(0).all().all())
            print("Real LDSC: three batches, nine pairs; positive tiny self-pair SEs preserved.")


if __name__ == "__main__":
    unittest.main()
