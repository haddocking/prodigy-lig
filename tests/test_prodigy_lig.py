#!/usr/bin/env python3

"""
Run the prodigy_lig tests.
"""

import json
import unittest
from os.path import join, basename, splitext, dirname

from Bio.PDB import PDBParser, FastMMCIFParser
from prodigy_lig import prodigy_lig

DATA_FOLDER_PATH = join(
    dirname(__file__),
    "..",
    "data"
)

with open(join(DATA_FOLDER_PATH, "test_data.json")) as ifile:
    TEST_DATA = json.load(ifile)


class ProdigyLigOutputTest(unittest.TestCase):
    """Test the output of prodigy_lig."""

    def test_dg_elec_file(self):
        """Test dg-elec & dg-score when reading HADDOCK file."""
        parser = PDBParser(QUIET=1)
        for pdb_file, values in TEST_DATA["elec"].items():
            struc_file_path = join(
                DATA_FOLDER_PATH,
                "elec",
                pdb_file + ".pdb"
            )

            with open(struc_file_path) as ifile:
                electrostatics = prodigy_lig.extract_electrostatics(ifile)

                pl = prodigy_lig.ProdigyLig(
                    parser.get_structure(pdb_file, ifile),
                    values["chains"],
                    electrostatics,
                    10.5
                )

                pl.predict()

                self.assertEqual(values["dg"], "{:.2f}".format(pl.dg_elec))
                self.assertEqual(values["score"], "{:.2f}".format(pl.dg_score))

    def test_dg_elec_cmdline(self):
        """Test dg-elec & dg-score with runtime-defined electrostatics."""
        parser = PDBParser(QUIET=1)
        for pdb_file, values in TEST_DATA["elec"].items():
            struc_file_path = join(
                DATA_FOLDER_PATH,
                "elec",
                pdb_file + ".pdb"
            )

            with open(struc_file_path) as ifile:
                pl = prodigy_lig.ProdigyLig(
                    parser.get_structure(pdb_file, ifile),
                    values["chains"],
                    values["electrostatics"],
                    10.5
                )

                pl.predict()

                self.assertEqual(values["dg"], "{:.2f}".format(pl.dg_elec))
                self.assertEqual(values["score"], "{:.2f}".format(pl.dg_score))

    def test_dg_file(self):
        """Test dg when reading file."""
        parser = PDBParser(QUIET=1)
        for pdb_file, values in TEST_DATA["no-elec"].items():
            struc_file_path = join(
                DATA_FOLDER_PATH,
                "no-elec",
                pdb_file + ".pdb"
            )

            with open(struc_file_path) as ifile:
                electrostatics = prodigy_lig.extract_electrostatics(ifile)

                pl = prodigy_lig.ProdigyLig(
                    parser.get_structure(pdb_file, ifile),
                    values["chains"],
                    electrostatics,
                    10.5
                )

                pl.predict()

                self.assertEqual(values["dg"], "{:.2f}".format(pl.dg))


def main():
    unittest.main()
