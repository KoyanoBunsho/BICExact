"""Regression tests; set BICEXACT_ESTIMATOR to a freshly built executable."""
import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
ESTIMATOR = Path(os.environ.get("BICEXACT_ESTIMATOR", ROOT / ".cache/estimate_hinge_numbers")).resolve()


class EstimatorSecurityTests(unittest.TestCase):
    def run_estimator(self, first, second, chain="A", method="exact"):
        return subprocess.run(
            [str(ESTIMATOR), str(first), str(second), chain, "A", "bic", method],
            capture_output=True, text=True, timeout=10,
        )

    def test_reference_structures(self):
        for method in ("exact", "lh"):
            with self.subTest(method=method):
                result = self.run_estimator(ROOT / "src/pdb3hvp.pdb", ROOT / "src/pdb4hvp.pdb", method=method)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("Protein length: 97", result.stdout)
                self.assertIn("Hinge index is: 45 : 57", result.stdout)

    def test_malformed_structures(self):
        for content in (
            "", "not a structure\n", "ATOM 1 CA GLY A\n",
            "ATOM 1 CA GLY A 1\n", "ATOM 1 CA GLY A 1 nan 2 3\n",
            "ATOM 1 CA GLY A 1 1e999 2 3\n", "ATOM 1 CA GLY A 1 1000001 2 3\n",
            "ATOM 1 CA GLY A 2147483648 1 2 3\n",
        ):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "input.pdb"
                path.write_text(content)
                result = self.run_estimator(path, path)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn(str(path), result.stderr)

    def test_chain_and_format_validation(self):
        sample = ROOT / "src/pdb3hvp.pdb"
        for chain in ("", "ABCDEFGHI", "A;id", "A\nB", "é", "Z"):
            with self.subTest(chain=chain):
                self.assertEqual(self.run_estimator(sample, sample, chain).returncode, 2)
        self.assertEqual(self.run_estimator(ROOT / "README.md", sample).returncode, 2)
        self.assertEqual(self.run_estimator(ROOT / "missing.pdb", sample).returncode, 2)

    def test_residue_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.pdb"
            path.write_text("".join(f"ATOM {i} CA GLY A {i} 1 2 3\n" for i in range(1, 10002)))
            self.assertEqual(self.run_estimator(path, path).returncode, 2)

    def test_missing_ca_coordinates_are_excluded(self):
        # An atom without a CA must not cause unequal matrices or an out-of-bounds read.
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "first.pdb", Path(tmp) / "second.pdb"
            common = "ATOM 1 CA GLY A 1 1 2 3\nATOM 2 CA GLY A 2 2 3 4\n"
            first.write_text(common + "ATOM 3 CA GLY A 3 3 4 5\n")
            second.write_text(common + "ATOM 3 N GLY A 3 3 4 5\n")
            result = self.run_estimator(first, second)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Protein length: 2", result.stdout)

    def test_cif_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.cif"
            path.write_text("".join(f"ATOM {i} C CA . GLY A 1 {i} ? {i}.0 2.0 3.0\n" for i in range(1, 4)))
            result = self.run_estimator(path, path)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Protein length: 3", result.stdout)


class ExperimentCommandSecurityTests(unittest.TestCase):
    def test_filenames_are_not_shell_commands(self):
        # Load only each runner function so the security regression needs no research dependencies.
        for filename, name, args in (
            ("fatcat_simulation.py", "run_fatcat", ("$(touch injected).pdb", "with spaces.pdb", "results")),
            ("fatcat_experiment_shibuya.py", "run_fatcat", ({"p_pdb": "$(touch injected)", "q_pdb": "with spaces"},)),
            ("dyndom_simulation.py", "run_dyndom", ("$(touch injected).command",)),
            ("dyndom_experiment_shibuya.py", "run_dyndom", ({"p_pdb": "$(touch injected)", "q_pdb": "with spaces"},)),
        ):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                tree = ast.parse((ROOT / "src" / filename).read_text())
                function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name)
                namespace = {"subprocess": subprocess, "time": time, "os": os}
                exec(compile(ast.Module(body=[function], type_ignores=[]), filename, "exec"), namespace)
                for executable in ("FATCAT", "DynDom"):
                    stub = Path(tmp) / executable
                    stub.write_text(f"#!{sys.executable}\nimport json,sys\nopen('args.json','w').write(json.dumps(sys.argv[1:]))\n")
                    stub.chmod(0o700)
                previous = Path.cwd()
                try:
                    os.chdir(tmp)
                    namespace[name](*args)
                    self.assertFalse(Path("injected").exists())
                    argv = json.loads(Path("args.json").read_text())
                    self.assertTrue(any("$(touch injected)" in arg for arg in argv))
                    if filename != "dyndom_simulation.py":
                        self.assertTrue(any("with spaces" in arg for arg in argv))
                finally:
                    os.chdir(previous)


if __name__ == "__main__":
    unittest.main()
