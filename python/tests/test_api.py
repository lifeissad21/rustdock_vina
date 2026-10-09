import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from rustdock_vina import Vina

ROOT = Path(__file__).resolve().parents[2]
BASIC = ROOT / "reference/example/basic_docking/solution"


class NativeApiTest(unittest.TestCase):
    def setUp(self):
        self.v = Vina(cpu=2, seed=42, verbosity=0)
        self.v.set_receptor(str(BASIC / "1iep_receptor.pdbqt"))
        self.v.set_ligand_from_file(str(BASIC / "1iep_ligand.pdbqt"))
        self.v.compute_vina_maps([15.19, 53.903, 16.917], [20, 20, 20], force_even_voxels=True)

    def tearDown(self):
        self.v.close()

    def test_score_optimize_and_input_validation(self):
        np.testing.assert_allclose(self.v.score(), [-12.513, -17.634, 0, 0, 0, -0.485, 5.121, -0.485], atol=0.001)
        self.assertAlmostEqual(self.v.optimize()[0], -13.170, places=3)
        with self.assertRaises(RuntimeError):
            self.v.set_ligand_from_string("malformed")
        self.assertEqual(self.v.info()["seed"], 42)
        with self.assertRaises(ValueError):
            self.v.compute_vina_maps([0, 0, 0], [1, 1, 1], spacing=0)

    def test_dock_pose_arrays_and_output(self):
        self.v.dock(exhaustiveness=2, n_poses=3, max_evals=2000)
        energies = self.v.energies(n_poses=3, energy_range=100)
        coordinates = self.v.poses(n_poses=3, energy_range=100, coordinates_only=True)
        self.assertEqual(energies.shape[1], 5)
        self.assertEqual(coordinates.shape, (energies.shape[0], 40, 3))
        self.assertTrue(np.isfinite(energies).all())
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / "poses.pdbqt")
            self.v.write_poses(path, n_poses=3, energy_range=100)
            self.assertIn("REMARK VINA RESULT:", Path(path).read_text())
            with self.assertRaises(RuntimeError):
                self.v.write_poses(path)
            self.v.write_poses(path, overwrite=True)
            self.v.write_pose(str(Path(tmp) / "best.pdbqt"))

    def test_maps_roundtrip_and_stateful_bridge(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = str(Path(tmp) / "maps")
            self.v.write_maps(prefix)
            with Vina(cpu=1, seed=42, no_refine=True, verbosity=0) as loaded:
                loaded.set_ligand_from_string((BASIC / "1iep_ligand.pdbqt").read_text())
                loaded.load_maps(prefix)
                self.assertAlmostEqual(loaded.score()[0], -11.557, delta=0.005)
                loaded.randomize(max_steps=5)
                loaded.write_pose(str(Path(tmp) / "random.pdbqt"))


class CliTest(unittest.TestCase):
    def test_batch_accepts_different_atom_types(self):
        binary = ROOT / "target/release" / ("vina.exe" if os.name == "nt" else "vina")
        if not binary.is_file():
            self.skipTest("release vina binary has not been built")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            atom = next(line for line in (BASIC / "1iep_ligand.pdbqt").read_text().splitlines() if line.startswith("ATOM"))
            paths = []
            for kind in ("C", "Cl"):
                path = tmp / (kind + ".pdbqt")
                path.write_text("ROOT\n" + atom[:77] + kind + "\nENDROOT\nTORSDOF 0\n")
                paths.append(str(path))
            result = subprocess.run([str(binary), "--receptor", str(BASIC / "1iep_receptor.pdbqt"), "--batch", *paths, "--dir", str(tmp), "--local_only", "--autobox", "--seed", "42", "--cpu", "1"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((tmp / "C_out.pdbqt").is_file())
            self.assertTrue((tmp / "Cl_out.pdbqt").is_file())


if __name__ == "__main__":
    unittest.main()
