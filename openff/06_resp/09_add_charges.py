#!/usr/bin/env python

from typing import List
import argparse
import pathlib
import itertools
import json
import tqdm
import os
import logging
import time
import glob
import sys

import numpy as np
import pandas as pd
from numpy.testing import assert_equal
from qcfractal import interface as ptl
from rdkit import Chem
from openff.toolkit.topology.molecule import Molecule, unit
from openff.toolkit.utils import OpenEyeToolkitWrapper, AmberToolsToolkitWrapper

import qcelemental as qcel

from psiresp.grid import GridOptions
from psiresp.orientation import OrientationEsp
from psiresp.constraint import ConstraintMatrix
from psiresp.resp import RespCharges, RespOptions
from psiresp.charge import MoleculeChargeConstraints, ChargeConstraintOptions
from psiresp import psi4utils

from qm import run_geometry_optimization, compute_energy
from waiter import wait_for_procedures, wait_for_results

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("smiles_index", type=int)
parser.add_argument("--server", default="localhost")
parser.add_argument("--port", default="7777")

METHODS = {
    "hf": {"method": "hf", "basis": "6-31g*"},
    "default": {"method": "b3lyp-d3bj", "basis": "dzvp"},
    "resp2": {"method": "PW6B95", "basis": "cc-pV(D+d)Z"},
}

METHOD_FULL = {
    k: METHODS[k]["method"] + "/" + METHODS[k]["basis"]
    for k in METHODS
}

GRIDS = {
    "hf": GridOptions(vdw_point_density=1, use_radii="msk"),
    "resp2": GridOptions(vdw_point_density=2.5, use_radii="bondi")
}

WEIGHTS = [1, None]

STATES = ["gas", "water"]

class Smiles:

    def __init__(self, smiles_index=0, n_orientations=12):
        self.smiles_index = int(smiles_index)
        self.n_combinations = n_orientations

        with open("01_mnsol_smiles.dat", "r") as f:
            contents = [x.strip() for x in f.readlines()]
        self.smiles = contents[smiles_index]
        logger.info(f"SMILES: {self.smiles}")
        print(f"SMILES: {self.smiles}")

        smiles_parser = Chem.rdmolfiles.SmilesParserParams()
        smiles_parser.removeHs = False
        self.rdmol = Chem.AddHs(Chem.MolFromSmiles(self.smiles, smiles_parser))
        self._offmol = Molecule.from_rdkit(self.rdmol, allow_undefined_stereo=True)

        self.outdir = pathlib.Path(f"08_results/{self.smiles}")
        self.load_orientations()

    def _load_orientations(self, key="hf"):
        for i in range(1, 11):
            orientations = []
            subfile = f'03*_c{i:02d}_o*{key}*.json'
            print(subfile)
            print(str(self.outdir / subfile))
            for file in sorted(glob.glob(str(self.outdir / subfile))):

                print(file)
                orientations.append(qcel.models.Molecule.from_file(file))
            self.orientations.append(orientations)

    def load_orientations(self):
        self.orientations = []
        for key in METHODS:
            self._load_orientations(key=key)
        

    def add_more_charges(self):
        original_file = self.outdir / "05_calculated_charges.csv"
        if not original_file.exists():
            sys.exit()

        outfile = self.outdir / "06_more_charges.csv"
        if outfile.exists():
            sys.exit()

        original_df = pd.read_csv(original_file, index_col=0)
        print(original_df)

        offmol = self.to_offmol()
        smirks = offmol.to_smiles(mapped=True)
        atom_number = list(range(1, len(offmol.atoms) + 1))
        symbols = [at.element.symbol for at in offmol.atoms]

        WRAPPERS = {
            "OpenEye": OpenEyeToolkitWrapper(),
            "AmberTools": AmberToolsToolkitWrapper(),
        }

        chgfile = self.outdir / f"{self.smiles}.json"
        with open(chgfile, "r") as f:
            charges = json.load(f)

        for i, conf_mols in tqdm.tqdm(enumerate(self.orientations, 1)):
            for j, orient in enumerate(conf_mols, 1):
                for tk, wrapper in WRAPPERS.items():
                    chg = self._get_am1bcc_charges(conf=i-1, orient=j-1, wrapper=wrapper)
                    dct = {"Charge model": "AM1", "Conformer": i, "Orientation": j, "Method": tk,
                           "Atom number": atom_number, "Element": symbols,
                           "charge": list(chg), "smirks": smirks}
                    charges.append(dct)
        
        
        with open(chgfile, "w") as f:
            json.dump(charges, f)
        print(f"Wrote to {chgfile}")

        dfs = [original_df]
        for dct in charges:
            df_ = pd.DataFrame({"Charge": dct.pop("charge")})
            for k, v in dct.items():
                df_[k] = v
            dfs.append(df_)
        df = pd.concat(dfs)
   
        df.to_csv(outfile)
        print(f"Wrote to {outfile}")

    def _get_am1bcc_charges(self, conf=0, orient=0, wrapper=None, **kwargs):
        offmol = self.to_offmol(conf=conf, orient=orient)
        offmol.assign_partial_charges("am1-mulliken", use_conformers=offmol.conformers, toolkit_registry=wrapper)
        charges = offmol.partial_charges
        charges = [x / unit.elementary_charge for x in charges]
        return charges


    def to_offmol(self, conf=0, orient=0):
        qcmol = self.orientations[conf][orient]
        
        rdel = [at.GetSymbol() for at in self.rdmol.GetAtoms()]
        assert_equal(rdel, qcmol.symbols) 

        offmol = Molecule.from_rdkit(self.rdmol, allow_undefined_stereo=True)
        offmol._conformers = []
        positions = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
        offmol._add_conformer(unit.Quantity(positions, unit.angstrom))

        return offmol

        
    def run(self):
        self.add_more_charges()



if __name__ == "__main__":
    args = parser.parse_args()

    smiler = Smiles(smiles_index=args.smiles_index)
    smiler.run()