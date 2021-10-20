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
        self.outdir.mkdir(exist_ok=True, parents=True)

    def compute_conformers(self):
        from rdutils import generate_diverse_conformer_ids, qcmol_from_rdkit

        ids = generate_diverse_conformer_ids(self.rdmol)
        self.conformers = [qcmol_from_rdkit(self.rdmol, confId=i) for i in ids]
        logger.info(f"Generated {len(self.conformers)} conformers")

        for i, qcmol in enumerate(self.conformers, 1):
            with open(self.outdir / f"01_conf-{i:02d}.json", "w") as f:
                f.write(qcmol.json())
    
    def opt_geometry(self, client):
        self.geometry_ids = {}

        for k, v in METHODS.items():
            self.geometry_ids[k] = run_geometry_optimization(client, self.conformers, **v)
            logger.info(f"Generated {self.geometry_ids[k]} geometry ids for {k}")

        with open(self.outdir / f"02_geometry_response_ids.json", "w") as f:
            json.dump(self.geometry_ids, f)

        ids = [x for v in self.geometry_ids.values() for x in v]
        print(f"Waiting for {len(ids)} geometry results")
        wait_for_procedures(client, response_ids=ids)

    def _generate_orientations(self, client, combinations, key="hf"):
        from orientations import generate_orientations
        

        results = client.query_procedures(id=self.geometry_ids[key])
        results = sorted(results, key=lambda x: x.id)
        optmols = [x.get_final_molecule() for x in results]

        response_ids = {"gas": [], "water": []}

        for i, optmol in enumerate(optmols, 1):
            optfile = self.outdir / f"02_c{i:02d}_{key}_opt.json"
            with open(optfile, "w") as f:
                f.write(optmol.json())

            orientation_qcmols = generate_orientations(optmol, combinations)
            self.orientations.append(orientation_qcmols)
            for j, qcmol in enumerate(orientation_qcmols, 1):
                molfile = self.outdir / f"03_c{i:02d}_o{j:02d}_{key}_mol.json"
                with open(molfile, "w") as f:
                    f.write(qcmol.json())

            for state in STATES:
                method_kws = METHODS[key]
                kw_ids = compute_energy(client, orientation_qcmols, **method_kws, state=state)
                response_ids[state].extend(kw_ids)
        
        return response_ids

    def generate_orientations(self, client):
        from orientations import generate_atom_combinations
        combinations = generate_atom_combinations(self.conformers[0], n_combinations=self.n_combinations)

        self.energy_ids = {}
        self.orientations = []

        for key in METHODS:
            self.energy_ids[key] = self._generate_orientations(client, combinations, key=key)
            logger.info(f"Generated {self.energy_ids[key]} orientation ids for {key}")

        with open(self.outdir / f"03_energy_response_ids.json", "w") as f:
            json.dump(self.energy_ids, f)

        ids = [y for v in self.energy_ids.values() for x in v.values() for y in x]
        print(f"Waiting for {len(ids)} orientation results")
        wait_for_results(client, response_ids=ids)
    

    def compute_charges(self, client):
        from esp import compute_charges, compute_orientation_esps
        n_atoms = len(self.conformers[0].geometry)
        atom_number = list(range(1, n_atoms + 1))
        qcmol = self.orientations[0][0]

        dfs = []

        groups = itertools.product(METHODS, STATES, GRIDS)
        for key, state, grid_name in tqdm.tqdm(list(groups)):
            response_ids = self.energy_ids[key]
            ids = response_ids[state]
            grid_options = GRIDS[grid_name]

            esps = compute_orientation_esps(client, ids, grid_options=grid_options)

            for weight in WEIGHTS:
                esps_ = [x.copy(deep=True) for x in esps]
                for esp in esps_:
                    esp.weight = weight

                stage_1_charges, stage_2_charges = compute_charges(qcmol, esps_)
                file1 = self.outdir / f"04_{key}_w-{weight}_g-{grid_name}_{state}_stage_1_respcharges.json"
                file2 = self.outdir / f"04_{key}_w-{weight}_g-{grid_name}_{state}_stage_2_respcharges.json"

                with open(str(file1), "w") as f:
                    f.write(stage_1_charges.json())

                with open(str(file2), "w") as f:
                    f.write(stage_2_charges.json())
                
                df_ = pd.DataFrame({"stage_1_charge": stage_1_charges._charges[:n_atoms]})
                df_["stage_2_charge"] = stage_2_charges._charges[:n_atoms]
                df_["grid"] = grid_name
                df_["weight"] = weight
                df_["state"] = state
                df_["method"] = METHOD_FULL[key]
                df_["atom_number"] = atom_number
                df_["smiles"] = self.smiles
                df_["mapped_smiles"] = self._offmol.to_smiles(mapped=True)
                dfs.append(df_)

        df = pd.concat(dfs)
        dfname = self.outdir / f"04_charges.csv"
        df.to_csv(dfname)
        logger.info(f"Wrote to {dfname}")
        print(f"Wrote to {dfname}")
        self.charge_df = df

    def to_offmol(self, conf=0, orient=0):
        qcmol = self.orientations[conf][orient]
        
        rdel = [at.GetSymbol() for at in self.rdmol.GetAtoms()]
        assert_equal(rdel, qcmol.symbols) 

        offmol = Molecule.from_rdkit(self.rdmol, allow_undefined_stereo=True)
        offmol._conformers = []
        positions = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
        offmol._add_conformer(unit.Quantity(positions, unit.angstrom))

        return offmol

    def _get_resp2_charges(self, delta=0.0, weight=1, key="hf", grid_name="hf"):
        # delta = % aqueous
        df = self.charge_df[self.charge_df.method == METHOD_FULL[key]]

        if weight is None:
            df = df[df.weight.isna()]
        else:
            df = df[df.weight == weight]
        df = df[df.grid == grid_name]
        gas = df[df.state == "gas"].stage_2_charge.values
        water = df[df.state == "water"].stage_2_charge.values

        charges = (delta * water) + ((1 - delta) * gas)

        offmol = self.to_offmol()
        
        offcharges = unit.Quantity(
            np.zeros(shape=len(charges), dtype=np.float64), unit.elementary_charge
        )
        for i, chg in enumerate(charges):
            offcharges[i] = chg * unit.elementary_charge
        offmol.partial_charges = offcharges
        offmol._normalize_partial_charges()
        charges = [x / unit.elementary_charge for x in offmol.partial_charges]

        return charges

    def _get_am1bcc_charges(self, conf=0, orient=0, wrapper=None, **kwargs):
        offmol = self.to_offmol(conf=conf, orient=orient)
        offmol.assign_partial_charges("am1bcc", use_conformers=offmol.conformers, toolkit_registry=wrapper)
        charges = offmol.partial_charges
        charges = [x / unit.elementary_charge for x in charges]
        return charges

    def generate_charges(self):
        offmol = self.to_offmol()
        smirks = offmol.to_smiles(mapped=True)
        atom_number = list(range(1, len(offmol.atoms) + 1))
        symbols = [at.element.symbol for at in offmol.atoms]

        charges = []
        for key in self.energy_ids:
            for grid_name in GRIDS:
                for weight in WEIGHTS:
                    for delta in np.linspace(0, 1, 21):
                        chg = self._get_resp2_charges(delta=delta, weight=weight, key=key, grid_name=grid_name)
                        dct = {"Charge model": "RESP2", "Method": METHOD_FULL[key], "Weight": weight, "Delta": delta,
                               "Element": symbols,
                               "Grid": grid_name, "Atom number": atom_number, "charge": list(chg), "smirks": smirks}
                        charges.append(dct)

        WRAPPERS = {
            "OpenEye": OpenEyeToolkitWrapper(),
            "AmberTools": AmberToolsToolkitWrapper(),
        }

        for i, conf_mols in enumerate(self.orientations, 1):
            for j, orient in enumerate(conf_mols, 1):
                for tk, wrapper in WRAPPERS.items():
                    chg = self._get_am1bcc_charges(conf=i-1, orient=j-1, wrapper=wrapper)
                    dct = {"Charge model": "AM1BCC", "Conformer": i, "Orientation": j, "Method": tk,
                           "Atom number": atom_number, "Element": symbols,
                           "charge": list(chg), "smirks": smirks}
                    charges.append(dct)
        
        chgfile = self.outdir / f"{self.smiles}.json"
        with open(chgfile, "w") as f:
            json.dump(charges, f)
        print(f"Wrote to {chgfile}")

        dfs = []
        for dct in charges:
            df_ = pd.DataFrame({"Charge": dct.pop("charge")})
            for k, v in dct.items():
                df_[k] = v
            dfs.append(df_)
        df = pd.concat(dfs)

        csvfile = self.outdir / "05_calculated_charges.csv"
        df.to_csv(csvfile)
        print(f"Wrote to {csvfile}")

        
    def run(self, client):
        self.compute_conformers()
        self.opt_geometry(client)
        self.generate_orientations(client)
        self.compute_charges(client)
        self.generate_charges()



if __name__ == "__main__":
    args = parser.parse_args()

    smiler = Smiles(smiles_index=args.smiles_index)
    csvfile = smiler.outdir / "04_charges.csv"
    if csvfile.exists():
        import sys
        sys.exit()

    client = ptl.FractalClient(address=f"{args.server}:{args.port}", verify=False)
    logger.info(f"client: {client}")
    smiler.run(client)