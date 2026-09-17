"""Machine-readable result compilation must not read or parse human logs."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import results
from ldsc_gpca.ldsc_export import write_results_csv


def record(p1="A", p2="A"):
    return {
        "p1": f"/inputs/{p1}.sumstats.gz", "p2": f"/inputs/{p2}.sumstats.gz",
        "rg": 0.9999999999999876 if p1 == p2 else 0.6123456789123456,
        "se": 1.9647123456789123e-06, "z": 508990.39251234567,
        "p": 3.212345678912345e-100,
        "h2_obs": 0.06351234567891234, "h2_obs_se": 0.015312345678912345,
        "h2_int": 1.0110123456789123, "h2_int_se": 6.212345678912345e-08,
        "gcov_int": 0.21234567891234567, "gcov_int_se": 6.212345678912345e-08,
    }


def manifest(names=("A",), references=None):
    return pd.DataFrame({
        "gwas_name": names, "ref": references or ["yes"] * len(names),
        "pop_prevalence": [float("nan")] * len(names),
    })


class NumericalCsvTests(unittest.TestCase):
    def compile(self, frame, metadata=None):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "batch.results.csv"
            write_results_csv(frame, source)
            with patch.object(results, "_read_legacy_log", side_effect=AssertionError("log read")), \
                    patch.object(results, "_detailed_correlations", side_effect=AssertionError("log parse")), \
                    patch.object(results.pd, "read_csv", wraps=pd.read_csv) as read:
                results.compile_results(directory, [source], manifest() if metadata is None else metadata)
                read.assert_called_once()
                self.assertEqual(read.call_args.args[0], source)
            return pd.read_csv(Path(directory) / "ldsc_results.csv", float_precision="round_trip")

    def test_round_trips_every_numerical_field_with_one_csv_read_and_no_logs(self):
        original = record()
        compiled = self.compile(pd.DataFrame([original]))
        for column, value in original.items():
            if column not in ("p1", "p2"):
                self.assertEqual(compiled[column].iloc[0], value, column)
        self.assertEqual(compiled.p1.iloc[0], "A")
        self.assertGreater(compiled.se.iloc[0], 0)

    def test_complete_shuffled_three_trait_matrix_and_matching_names(self):
        records = [record(a, b) for a in ("A", "B", "C") for b in ("A", "B", "C")]
        compiled = self.compile(
            pd.DataFrame(records).sample(frac=1, random_state=17), manifest(("C", "A", "B"))
        )
        self.assertEqual(len(compiled), 9)
        self.assertEqual(set(zip(compiled.p1, compiled.p2)),
                         {(a, b) for a in ("A", "B", "C") for b in ("A", "B", "C")})

    def test_preserves_existing_observed_and_liability_mapping(self):
        observed = pd.DataFrame([record("A", "A")])
        liability = pd.DataFrame([record("A", "B")]).rename(
            columns={"h2_obs": "h2_liab", "h2_obs_se": "h2_liab_se"}
        )
        metadata = manifest(("A", "B"), ["yes", "no"])
        metadata.loc[1, "pop_prevalence"] = 0.01
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / "observed.results.csv", Path(directory) / "liability.results.csv"]
            write_results_csv(observed, paths[0])
            write_results_csv(liability, paths[1])
            with patch.object(results.pd, "read_csv", wraps=pd.read_csv) as read:
                results.compile_results(directory, paths, metadata)
                self.assertEqual(read.call_count, 2)
            compiled = pd.read_csv(Path(directory) / "ldsc_results.csv", float_precision="round_trip")
        self.assertEqual(compiled.h2_scale.tolist(), ["NEF_unconverted", "liability"])
        self.assertEqual(compiled.h2_obs.iloc[0], observed.h2_obs.iloc[0])
        self.assertEqual(compiled.h2_liab.iloc[1], liability.h2_liab.iloc[0])
        self.assertTrue(pd.isna(compiled.h2_liab.iloc[0]))

    def test_missing_csv_never_falls_back_to_existing_log(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "batch.log").write_text("readable log exists")
            with patch.object(results, "_read_legacy_log") as legacy:
                with self.assertRaisesRegex(RuntimeError, "No readable-log fallback"):
                    results.compile_results(directory, [Path(directory, "batch.results.csv")], manifest())
                legacy.assert_not_called()
            self.assertFalse(Path(directory, "ldsc_results.csv").exists())

    def test_empty_or_malformed_csv_fails_without_log_fallback(self):
        for content in ("", 'p1,p2,rg\n"unterminated,A,1\n'):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory, "batch.results.csv")
                path.write_text(content)
                with patch.object(results, "_read_legacy_log") as legacy:
                    with self.assertRaisesRegex(RuntimeError, "No readable-log fallback"):
                        results.compile_results(directory, [path], manifest())
                    legacy.assert_not_called()

    def test_p_boundaries_and_smallest_positive_se_are_not_clipped(self):
        for p in (0.0, 1.0):
            row = record()
            row.update(p=p, se=float.fromhex("0x0.0000000000001p-1022"))
            with self.subTest(p=p):
                observed = self.compile(pd.DataFrame([row]))
                self.assertEqual(observed.p.iloc[0], p)
                self.assertEqual(observed.se.iloc[0], row["se"])

    def test_invalid_numeric_and_missing_values_fail_without_imputation(self):
        for column, value, message in (
            ("se", 0, "Non-positive"), ("se", -1e-6, "Non-positive"),
            ("se", float("nan"), "Non-finite"), ("se", float("inf"), "Non-finite"),
            ("h2_obs_se", 0, "Non-positive"), ("h2_int_se", 0, "Non-positive"),
            ("gcov_int_se", 0, "Non-positive"), ("z", float("inf"), "Non-finite"),
            ("rg", "bad", "Non-finite"), ("p", -0.1, "outside"),
            ("p", 1.1, "outside"), ("p1", None, "identifier"),
        ):
            row = record()
            row[column] = value
            with self.subTest(column=column, value=value), self.assertRaisesRegex(RuntimeError, message):
                self.compile(pd.DataFrame([row]))

    def test_missing_schema_empty_file_and_incomplete_pairs_fail(self):
        for frame, metadata, message in (
            (pd.DataFrame([record()]).drop(columns="se"), manifest(), "Missing LDSC result columns"),
            (pd.DataFrame([record()]).drop(columns="h2_obs_se"), manifest(), "complete h2_obs"),
            (pd.DataFrame(columns=list(record())), manifest(), "No correlation results"),
            (pd.DataFrame([record()]), manifest(("A", "B")), "missing="),
            (pd.DataFrame([record("A", "C")]), manifest(), "unexpected="),
        ):
            with self.subTest(message=message), self.assertRaisesRegex(RuntimeError, message):
                self.compile(frame, metadata)

    def test_consistent_duplicates_collapse_conflicting_records_are_not_silently_discarded(self):
        self.assertEqual(len(self.compile(pd.DataFrame([record(), record()]))), 1)
        different = record()
        different["rg"] = 0.9
        # The downstream R duplicate-conflict checks remain authoritative.
        compiled = self.compile(pd.DataFrame([record(), different]))
        self.assertEqual(len(compiled), 2)
        self.assertEqual(compiled.rg.iloc[1], 0.9)


if __name__ == "__main__":
    unittest.main()
