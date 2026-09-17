"""The export hook observes the native table once; it never recomputes estimates."""
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
from ldsc_gpca import ldsc_export, ldsc_runtime


def frame():
    return pd.DataFrame([{
        "p1": "/inputs/A.sumstats.gz", "p2": "/inputs/A.sumstats.gz",
        "rg": 0.9999999999999876, "se": 1.9647123456789123e-06,
        "z": 508990.39251234567, "p": 3.212345678912345e-100,
        "h2_obs": 0.06351234567891234, "h2_obs_se": 0.015312345678912345,
        "h2_int": 1.0110123456789123, "h2_int_se": 6.212345678912345e-08,
        "gcov_int": 1.0110123456789123, "gcov_int_se": 6.212345678912345e-08,
    }])


class ExportTests(unittest.TestCase):
    def test_exports_same_frame_once_and_preserves_native_rendering(self):
        data = frame()
        calls = []
        native_to_string = pd.DataFrame.to_string

        def native_table(paths, estimates, args):
            calls.append((paths, estimates))
            return data.to_string(header=True, index=False) + "\n"

        module = SimpleNamespace(_get_rg_table=native_table)
        with tempfile.TemporaryDirectory() as directory, pd.option_context(
            "display.float_format", "{:.4f}".format
        ):
            args = SimpleNamespace(out=str(Path(directory) / "batch"))
            expected_log = native_to_string(data, header=True, index=False) + "\n"
            with ldsc_export.export_rg_tables(module, pd) as exported:
                rendered = module._get_rg_table(["A", "A"], ["already computed"], args)
                self.assertEqual(rendered, expected_log)
                self.assertEqual(exported, [args.out + ".results.csv"])
            self.assertEqual(len(calls), 1)
            self.assertIs(module._get_rg_table, native_table)
            self.assertIs(pd.DataFrame.to_string, native_to_string)
            observed = pd.read_csv(exported[0], float_precision="round_trip")
            pd.testing.assert_frame_equal(observed, data)
            self.assertIn("0.0000", rendered)
            self.assertGreater(observed.se.iloc[0], 0)

    def test_liability_scale_and_failed_estimates_are_exported_without_imputation(self):
        data = frame().rename(columns={"h2_obs": "h2_liab", "h2_obs_se": "h2_liab_se"})
        failed = data.copy()
        failed.loc[0, ["rg", "se", "z", "p"]] = float("nan")
        data = pd.concat([data, failed], ignore_index=True)
        native_table = lambda paths, estimates, args: data.to_string(index=False)
        module = SimpleNamespace(_get_rg_table=native_table)
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(out=str(Path(directory) / "batch"))
            with ldsc_export.export_rg_tables(module, pd) as exported:
                module._get_rg_table([], [], args)
            observed = pd.read_csv(exported[0], float_precision="round_trip")
            pd.testing.assert_frame_equal(observed, data)

    def test_export_failure_restores_hooks_and_preserves_previous_csv(self):
        module = SimpleNamespace(_get_rg_table=lambda paths, estimates, args: frame().to_string())
        native_table, native_to_string = module._get_rg_table, pd.DataFrame.to_string
        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(out=str(Path(directory) / "batch"))
            destination = Path(args.out + ".results.csv")
            destination.write_text("previous results")
            with patch.object(pd.DataFrame, "to_csv", side_effect=OSError("disk full")):
                with self.assertRaisesRegex(OSError, "disk full"):
                    with ldsc_export.export_rg_tables(module, pd):
                        module._get_rg_table([], [], args)
            self.assertEqual(destination.read_text(), "previous results")
            self.assertEqual(list(Path(directory).glob(".ldsc-results-*")), [])
        self.assertIs(module._get_rg_table, native_table)
        self.assertIs(pd.DataFrame.to_string, native_to_string)

    def test_unsupported_table_and_no_capture_fail(self):
        for getter, message in (
            (lambda *args: pd.DataFrame({"wrong": [1]}).to_string(), "schema"),
            (lambda *args: "already formatted", "exactly one"),
        ):
            module = SimpleNamespace(_get_rg_table=getter)
            with self.subTest(message=message), self.assertRaisesRegex(RuntimeError, message):
                with ldsc_export.export_rg_tables(module, pd):
                    module._get_rg_table([], [], SimpleNamespace(out="unused"))
            self.assertIs(module._get_rg_table, getter)
        with self.assertRaisesRegex(RuntimeError, "Unsupported LDSC runtime"):
            with ldsc_export.export_rg_tables(SimpleNamespace(), pd):
                pass

    def test_main_launches_one_native_job_and_restores_argv(self):
        module = types.ModuleType("ldscore.sumstats")
        module._get_rg_table = lambda paths, estimates, args: frame().to_string(index=False)
        parent = types.ModuleType("ldscore")
        parent.sumstats = module
        original_argv = sys.argv
        with tempfile.TemporaryDirectory() as directory:
            prefix = str(Path(directory) / "batch")

            def native_job(path, run_name):
                self.assertEqual(sys.argv, ["/runtime/ldsc.py", "--rg", "A,A", "--out", prefix])
                module._get_rg_table(["A", "A"], [], SimpleNamespace(out=prefix))

            with patch.dict(sys.modules, {"ldscore": parent, "ldscore.sumstats": module}), \
                    patch.object(ldsc_export.shutil, "which", return_value="/runtime/ldsc.py"), \
                    patch.object(ldsc_export.runpy, "run_path", side_effect=native_job) as run:
                self.assertEqual(ldsc_export.main(["--rg", "A,A", "--out", prefix]), 0)
                run.assert_called_once_with("/runtime/ldsc.py", run_name="__main__")
            self.assertTrue(Path(prefix + ".results.csv").is_file())
        self.assertIs(sys.argv, original_argv)

    def test_main_rejects_missing_fresh_export_even_when_stale_csv_exists(self):
        module = types.ModuleType("ldscore.sumstats")
        module._get_rg_table = lambda *args: ""
        parent = types.ModuleType("ldscore")
        parent.sumstats = module
        with tempfile.TemporaryDirectory() as directory:
            prefix = str(Path(directory) / "batch")
            Path(prefix + ".results.csv").write_text("stale")
            with patch.dict(sys.modules, {"ldscore": parent, "ldscore.sumstats": module}), \
                    patch.object(ldsc_export.shutil, "which", return_value="/runtime/ldsc.py"), \
                    patch.object(ldsc_export.runpy, "run_path", return_value={}):
                with self.assertRaisesRegex(RuntimeError, "fresh numerical CSV"):
                    ldsc_export.main(["--rg", "A,A", "--out", prefix])

    def test_missing_native_script_fails(self):
        with patch.object(ldsc_export.shutil, "which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "not installed"):
                ldsc_export.main(["--help"])

    def test_runtime_uses_existing_child_process_not_a_second_job(self):
        command = ldsc_runtime.ldsc_regression_command(conda="/tools/conda", prefix="/child env")
        self.assertEqual(command[:6], [
            "/tools/conda", "run", "--no-capture-output", "--prefix", "/child env", "python",
        ])
        self.assertEqual(Path(command[6]).name, "ldsc_export.py")
        self.assertTrue(os.path.isfile(command[6]))


if __name__ == "__main__":
    unittest.main()
