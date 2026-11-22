#!/usr/bin/env python3
"""Unit tests for inspect_complex.py classification and detection."""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(KIT, "scripts"))
import inspect_complex as ic  # noqa: E402

DATA = os.path.join(HERE, "data")


def opts(**overrides):
    base = {
        "model": 1, "altloc": "most-populous", "type": "auto",
        "ligand": None, "receptor_chains": None, "partner_chains": None,
        "drop_glycans": False, "drop_metals": False, "drop_modified": False,
        "drop_unknown_hetatm": False, "drop_other_ligands": False,
        "strict_gaps": False, "peptide_max": 50, "min_ligand_heavy": 3,
        "max_ligand_heavy": 200, "contact_cutoff": 4.0,
        "covalent_cutoff": 1.9, "ss_cutoff": 2.5, "cn_break_cutoff": 2.0,
    }
    base.update(overrides)
    return base


def run(name, **overrides):
    return ic.inspect(os.path.join(DATA, name), opts(**overrides))


class DetectionTests(unittest.TestCase):
    def test_protein_ligand(self):
        data = run("auto_pl.pdb")
        cls = data["classification"]
        self.assertEqual(cls["type"], "protein-ligand")
        self.assertEqual(cls["errors"], [])
        self.assertIsNotNone(cls["ligand"])
        self.assertEqual(cls["ligand"]["resname"], "LIG")

    def test_protein_peptide(self):
        data = run("auto_pep.pdb")
        cls = data["classification"]
        self.assertEqual(cls["type"], "protein-peptide")
        self.assertEqual(cls["partner"]["chains"], ["B"])
        self.assertEqual(cls["receptor"]["chains"], ["A"])

    def test_antibody(self):
        data = run("auto_antibody.pdb")
        cls = data["classification"]
        self.assertEqual(cls["type"], "protein-antibody")
        self.assertEqual(sorted(cls["partner"]["chains"]), ["A", "B"])
        self.assertEqual(cls["receptor"]["chains"], ["C"])

    def test_glycan_error_and_override(self):
        data = run("auto_glycan.pdb")
        codes = [e["code"] for e in data["classification"]["errors"]]
        self.assertIn("glycan_unsupported", codes)
        data = run("auto_glycan.pdb", drop_glycans=True)
        self.assertEqual(data["classification"]["errors"], [])
        self.assertTrue(data["glycans"]["detected"])

    def test_metal_error_and_override(self):
        data = run("auto_zn.pdb")
        codes = [e["code"] for e in data["classification"]["errors"]]
        self.assertIn("metal_unsupported", codes)
        data = run("auto_zn.pdb", drop_metals=True)
        self.assertEqual(data["classification"]["errors"], [])

    def test_mse_convertible(self):
        data = run("auto_mse.pdb")
        self.assertTrue(data["residues"]["modified_convertible"])
        self.assertEqual(data["residues"]["modified_convertible"][0]["resname"], "MSE")

    def test_disulfide(self):
        data = run("auto_ss.pdb")
        self.assertEqual(len(data["disulfides"]), 1)
        self.assertLess(data["disulfides"][0]["sg_distance_A"], 2.5)

    def test_chain_break_same_chain_is_apo(self):
        data = run("auto_gap.pdb")
        self.assertEqual(data["classification"]["type"], "apo-protein")
        self.assertEqual(len(data["chain_breaks"]), 1)
        self.assertGreater(data["chain_breaks"][0]["cn_distance_A"], 2.0)

    def test_ligand_ambiguity_requires_selection(self):
        # Two copies of the same ligand should be rejected without a selector.
        with open(os.path.join(DATA, "auto_pl.pdb")) as fh:
            text = fh.read()
        lig_lines = [l for l in text.splitlines() if l.startswith("HETATM")]
        second_copy = []
        for line in lig_lines:
            # same pose, different residue number: both copies contact protein
            second_copy.append(line[:22] + "   2" + " " + line[27:])
        second = "\n".join(second_copy)
        tmp = os.path.join("/tmp", "auto_two_ligands.pdb")
        with open(tmp, "w") as fh:
            fh.write(text.rstrip("\n") + "\n" + second + "\nTER\nEND\n")
        data = ic.inspect(tmp, opts())
        codes = [e["code"] for e in data["classification"]["errors"]]
        self.assertTrue("multiple_ligand_copies" in codes
                        or "ligand_ambiguous" in codes)


if __name__ == "__main__":
    unittest.main(verbosity=2)
