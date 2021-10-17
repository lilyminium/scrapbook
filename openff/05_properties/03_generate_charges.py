#!/usr/bin/env python

import pathlib
import json
import glob
import os
import argparse

import pandas as pd
import numpy as np
import qcelemental as qcel
from rdkit import Chem
from numpy.testing import assert_equal

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.evaluator.backends import ComputeResources
from openff.evaluator.backends.dask import DaskLocalCluster
from openff.toolkit.topology.molecule import Molecule, unit
from openff.evaluator.workflow import Workflow
from openff.toolkit.utils import OpenEyeToolkitWrapper, AmberToolsToolkitWrapper


from utils import property_from_csv

parser = argparse.ArgumentParser("Evaluate properties")
parser.add_argument("indir", type=str)

class TestMolecule:

    def __init__(self, molname="0010", component="solvent", directory="../06_resp/07_results",
                 csv="01_mnsol_data.csv"):
        self.molname = molname
        self.mol_index = int(molname)
        self.component = component

        self.indir = pathlib.Path(directory) / f"mol-{molname}_{component}"
        self.smiles = self.get_smiles()
        self.load_charges()
        self.load_energies()
    
    def load_charges(self, pattern=r"04_w-{weight}_stage-{stage}_{state}_charges.dat"):
        qcmol = self.to_qcmol()
        n_atom = len(qcmol.geometry)
        
        df_dct = {"State": [], "Weight": [], "Stage": [], "Atom": [], "Charge": []}

        for state in ["gas", "water"]:
            for weight in [1, None]:
                for stage in 1, 2:
                    filename = self.indir / pattern.format(stage=stage, weight=weight, state=state)
                    charges = np.loadtxt(str(filename))[:n_atom]
                    setattr(self, f"stage_{stage}_weight_{weight}_{state}_charges", charges)

                    df_dct["Charge"].extend(charges)
                    df_dct["Atom"].extend(np.arange(n_atom))
                    df_dct["Stage"].extend([stage] * n_atom)
                    df_dct["Weight"].extend([weight] * n_atom)
                    df_dct["State"].extend([state] * n_atom)
        
        df = pd.DataFrame(df_dct)
        df.to_csv(f"03_results/resp_charges/{self.smiles}.csv")
        self.charge_df = df

    
    def load_energies(self, pattern="03_c*energy.json"):
        from qcfractal.interface.models.records import ResultRecord
        
        files = sorted(glob.glob(str(self.indir / pattern)))

        self.records = [ResultRecord.parse_file(file) for file in files]
        self.energies = [rec.properties.return_energy for rec in self.records]

    def get_charges(self, delta=0.0, weight=1):
        # delta = % aqueous
        gas = getattr(self, f"stage_2_weight_{weight}_gas_charges")
        water = getattr(self, f"stage_2_weight_{weight}_water_charges")

        charges = (delta * water) + ((1 - delta) * gas)
        return charges

    def get_smiles(self):
        df = pd.read_csv("01_mnsol_data.csv")
        row = df.iloc[self.mol_index]
        if row["Role 1"].lower() == self.component:
            smiles = row["Component 1"]
        else:
            smiles = row["Component 2"]
        return smiles

    def to_qcmol(self, pattern=r"03_c{conf:02d}_o{orient:02d}_mol.json", conf=1, orient=1):
        file = self.indir / pattern.format(conf=conf, orient=orient)
        qcmol = qcel.models.Molecule.from_file(file)
        return qcmol

    def to_offmol(self, conf=1, orient=1):
        qcmol = self.to_qcmol(conf=conf, orient=orient)

        smiles_parser = Chem.rdmolfiles.SmilesParserParams()
        smiles_parser.removeHs = False
        rdmol = Chem.AddHs(Chem.MolFromSmiles(self.smiles, smiles_parser))
        rdel = [at.GetSymbol() for at in rdmol.GetAtoms()]
        assert_equal(rdel, qcmol.symbols)

        offmol = Molecule.from_rdkit(rdmol)
        positions = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
        offmol._add_conformer(unit.Quantity(positions, unit.angstrom))

        return offmol

    def generate_resp2_charges(self, conf=1, orient=1, delta=0.0, weight=1):
        offmol = self.to_offmol(conf=conf, orient=orient)
        charges = self.get_charges(delta=delta, weight=weight)

        offcharges = unit.Quantity(
            np.zeros(shape=len(charges), dtype=np.float64), unit.elementary_charge
        )
        for i, chg in enumerate(charges):
            offcharges[i] = chg * unit.elementary_charge
        offmol.partial_charges = offcharges
        offmol._normalize_partial_charges()
        charges = [x / unit.elementary_charge for x in offmol.partial_charges]
        smiles = offmol.to_smiles(mapped=True)

        return {"charge": list(charges), "smirks": smiles}

    def generate_am1bcc_charges(self, conf=1, orient=1, wrapper=None, **kwargs):
        offmol = self.to_offmol(conf=conf, orient=orient)
        offmol.assign_partial_charges("am1bcc", use_conformers=offmol.conformers, toolkit_registry=wrapper)
        charges = offmol.partial_charges
        charges = [x / unit.elementary_charge for x in charges]
        smiles = offmol.to_smiles(mapped=True)

        return {"charge": list(charges), "smirks": smiles}

    def generate_charges(self):
        offmol = self.to_offmol(conf=1, orient=1)
        smiles = self.get_smiles()

        charges = []
        for weight in [1, None]:
            for delta in np.linspace(0, 1, 21):
                
                dct = {"Charge model": "RESP2", "Weight": weight, "Delta": delta}
                dct.update(self.generate_resp2_charges(delta=delta, weight=weight))
                charges.append(dct)

        WRAPPERS = {
            "OpenEye": OpenEyeToolkitWrapper(),
            "AmberTools": AmberToolsToolkitWrapper(),
        }

        for conf in range(1, 11):
            for orient in range(1, 13):
                for tk, wrapper in WRAPPERS.items():
                    dct = {"Charge model": "AM1BCC", "Conformer": conf, "Orientation": orient, "Toolkit": tk}
                    dct.update(self.generate_am1bcc_charges(conf=conf, orient=orient, wrapper=wrapper))
                    charges.append(dct)
        
        with open(f"03_results/{smiles}.json", "w") as f:
            json.dump(charges, f)


if __name__ == "__main__":
    args = parser.parse_args()

    # set up molecule
    # molname = f"{args.mol_index:04d}"
    molname = args.indir.split("mol-")[-1].split("_")[0]
    component = args.indir.split("mol-")[-1].split("_")[1].split("/")[0]
    mol = TestMolecule(molname=molname, component=component)
    mol.generate_charges()