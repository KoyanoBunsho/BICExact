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

    def test_fixed_width_residue_matching(self):
        # Adjacent negative fields fill all eight columns. Residues 1000 and
        # 1000A are distinct, and the two files list them in opposite orders.
        ids = [(999, " "), (1000, " "), (1000, "A"), (1001, " ")]
        xyz = [(-110.125, -210.250, -310.5), (10, 20, 30), (0, 5, -8), (20, 1, 7)]
        records = [
            f"ATOM  {i:5d}  CA  GLY A{n:4d}{code}   {x:8.3f}{y:8.3f}{z:8.3f}\n"
            for i, ((n, code), (x, y, z)) in enumerate(zip(ids, xyz), 1)
        ]
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "a.PDB", Path(tmp) / "b.ENT"
            first.write_text("".join(records))
            second.write_text("".join(reversed(records)))
            for method in ("exact", "lh"):
                result = self.run_estimator(first, second, method=method)
                self.assertEqual(result.returncode, 0, result.stderr)
                fields = dict(line.split(": ", 1) for line in result.stdout.splitlines())
                self.assertEqual(fields["Protein length"], "4")
                self.assertEqual(fields["The number of hinges is"], "0")
                self.assertLess(float(fields["RMSD"]), 1e-5)
                self.assertLess(float(fields["RMSDhk is"]), 1e-5)

    def test_cif_declared_columns_and_wrapped_rows(self):
        header = """data_example
_note.text
;
ATOM is metadata here
;
loop_
_atom_site.Cartn_z
_atom_site.label_seq_id
_atom_site.group_PDB
_atom_site.Cartn_y
_atom_site.label_asym_id
_atom_site.label_atom_id
_atom_site.Cartn_x
"""
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "a.CIF", Path(tmp) / "b.pdb"
            second.write_text("ATOM 1 CA GLY A 3 1 2 3\nATOM 2 CA GLY A 5 2 4 9\n")
            rows = """9 5 ATOM 4 'AB' "CA" 2
3 3 ATOM
2 AB CA 1
"""
            first.write_text(header + rows)
            result = self.run_estimator(first, second, chain="AB")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Protein length: 2", result.stdout)
            self.assertLess(float(result.stdout.split("RMSD: ")[1].splitlines()[0]), 1e-6)
            for bad_rows in (rows + "9 6 ATOM", rows.replace("CA 1", "CA 1xyz")):
                first.write_text(header + bad_rows)
                self.assertEqual(self.run_estimator(first, second, chain="AB").returncode, 2)

    def test_mixed_format_author_numbering_and_first_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb, cif = Path(tmp) / "a.pdb", Path(tmp) / "b.cif"
            pdb.write_text("MODEL        1\n"
                           "ATOM      1  CA  GLY A  10A     10.000  20.000  30.000\n"
                           "ENDMDL\nMODEL        2\n"
                           "ATOM      2  CA  GLY A  11      10.000  20.000  30.000\nENDMDL\n")
            tags = ("group_PDB", "label_atom_id", "label_asym_id", "label_seq_id",
                    "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
                    "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num")
            cif.write_text("data_test\nloop_\n" + "".join(f"_atom_site.{tag}\n" for tag in tags)
                           + "ATOM CA Z 1 A 10 A 10 20 30 1\nATOM CA Z 2 A 11 ? 10 20 30 2\n")
            result = self.run_estimator(pdb, cif)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Protein length: 1", result.stdout)
            self.assertIn("The number of hinges is: 0", result.stdout)

    def test_reflection_is_not_a_rotation(self):
        # This centered tetrahedron has covariance 4I. Mirroring one axis
        # gives a minimum squared deviation of 16 under proper rotations.
        tetrahedron = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp) / name for name in ("first.pdb", "second.pdb")]
            for sign, path in zip((1, -1), paths):
                # Small coordinates keep the penalized optimum at zero hinges.
                path.write_text("".join(f"ATOM {i} CA GLY A {i} {sign*x/10} {y/10} {z/10}\n"
                                        for i, (x, y, z) in enumerate(tetrahedron, 1)))
            for method in ("exact", "lh"):
                result = self.run_estimator(*paths, method=method)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("RMSDhk is: 0.2\n", result.stdout)


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
