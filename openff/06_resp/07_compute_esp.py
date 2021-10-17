#!/usr/bin/env python

from typing import List
import argparse
import pathlib
import itertools
import json
import time
import glob

import numpy as np
from qcfractal import interface as ptl

from psiresp.grid import GridOptions
from psiresp.orientation import OrientationEsp
from psiresp.constraint import ConstraintMatrix
from psiresp.resp import RespCharges, RespOptions
from psiresp.charge import MoleculeChargeConstraints, ChargeConstraintOptions
from psiresp.molecule import Molecule
from psiresp import psi4utils

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("--mol_index", type=int, default=10)
parser.add_argument("--component", default="solvent")
parser.add_argument("--server", default="localhost")
parser.add_argument("--port", default="7777")


GRID_OPTIONS = GridOptions(vdw_point_density=2.5, use_radii="bondi")
RESP_OPTIONS = RespOptions()
CHARGE_CONSTRAINT_OPTIONS = ChargeConstraintOptions()

def compute_esps(records, weight=None):
    esps = []
    for rec in records:
        qcmol = rec.get_molecule()
        grid = GRID_OPTIONS.generate_vdw_grid(qcmol)
        esp = OrientationEsp(
            qcmol=qcmol,
            grid=grid,
            esp=psi4utils.compute_esp(rec, grid),
            energy=rec.properties.return_energy,
            weight=weight
        )
        esps.append(esp)
    return esps


def create_constraint_matrix(mol, esps):
    matrix = ConstraintMatrix.with_n_dim(mol.qcmol.geometry.shape[0])
    for esp in esps:
        mat = esp.get_weighted_matrix()
        matrix += mat
    return matrix


def compute_charges(client, ids, weight=None):
    records = client.query_results(id=ids)
    esps = compute_esps(records, weight=weight)
    qcmol = records[0].get_molecule()
    mol = Molecule(qcmol=qcmol)
    surface_constraints = create_constraint_matrix(mol, esps)
    print("surface constraints", surface_constraints)

    stage_1_constraints = MoleculeChargeConstraints.from_charge_constraints(CHARGE_CONSTRAINT_OPTIONS, molecules=[mol])
    if RESP_OPTIONS.stage_2:
        stage_2_constraints = stage_1_constraints.copy(deep=True)
        stage_1_constraints.charge_equivalence_constraints = []

    stage_1_charges = RespCharges(charge_constraints=stage_1_constraints,
                                  surface_constraints=surface_constraints,
                                  resp_a=RESP_OPTIONS.resp_a1,
                                  **RESP_OPTIONS._base_kwargs)
    stage_1_charges.solve()

    if RESP_OPTIONS.stage_2:
        stage_2_constraints = MoleculeChargeConstraints(molecules=[mol])
        stage_2_constraints.add_constraints_from_charges(stage_1_charges._charges)
        stage_2_charges = RespCharges(charge_constraints=stage_2_constraints,
                                      surface_constraints=surface_constraints,
                                      resp_a=RESP_OPTIONS.resp_a2,
                                      **RESP_OPTIONS._base_kwargs)
        stage_2_charges.solve()

        return stage_1_charges, stage_2_charges
    return stage_1_charges


if __name__ == "__main__":
    args = parser.parse_args()
    client = ptl.FractalClient(address=f"{args.server}:{args.port}", verify=False)
    print(f"client: {client}")

    molname = f"mol-{args.mol_index:04d}_{args.component}"
    response_files = glob.glob(f"{molname}/03_c*_response_ids.json")

    conformer_response_ids = {"gas": [], "water": []}
    for file in response_files:
        print(file)
        with open(file, "r") as f:
            response_ids = json.load(f)
            for k, v in response_ids.items():
                conformer_response_ids[k].extend(v)
    
    print("-- conformer ids --")
    print(conformer_response_ids)


    for state, ids in conformer_response_ids.items():
        for weight in [1, None]:
            stage_1_charges, stage_2_charges = compute_charges(client, ids, weight=weight)
            file1 =  f"{molname}/04_w-{weight}_stage-1_{state}_charges.dat"
            np.savetxt(str(file1), stage_1_charges._charges)
            print(f"Wrote to {file1}")
            file2 = f"{molname}/04_w-{weight}_stage-2_{state}_charges.dat"
            np.savetxt(str(file2), stage_2_charges._charges)
            print(f"Wrote to {file2}")

            file1 = f"{molname}/04_w-{weight}_stage-1_{state}_respcharges.json"
            file2 = f"{molname}/04_w-{weight}_stage-2_{state}_respcharges.json"

            with open(str(file1), "w") as f:
                f.write(stage_1_charges.json())

            print(f"Wrote to {file1}")

            print(stage_2_charges)
            print("------")
            # print(stage_2_charges.keys())
            print("------")
            print(stage_2_charges.charge_constraints)

            print(".>>")
            print(list(stage_2_charges.charge_constraints.charge_sum_constraints[0].atoms)[0].dict())



            first = stage_2_charges.charge_constraints.charge_sum_constraints[0]

            for k, v in first.__dict__.items():
                print(k)
                print(v)
            print(first.dict())


            print("------")
            print(stage_2_charges.surface_constraints.dict())

            with open(str(file2), "w") as f:
                f.write(stage_2_charges.json())

            print(f"Wrote to {file2}")

