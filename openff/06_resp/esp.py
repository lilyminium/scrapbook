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

RESP_OPTIONS = RespOptions()
CHARGE_CONSTRAINT_OPTIONS = ChargeConstraintOptions()

def compute_esps(records, grid_options=None, weight=None):
    esps = []
    for rec in records:
        qcmol = rec.get_molecule()
        grid = grid_options.generate_vdw_grid(qcmol)
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


def compute_charges(client, ids, weight=None, grid_options=None):
    records = client.query_results(id=ids)
    esps = compute_esps(records, grid_options=grid_options, weight=weight)
    qcmol = records[0].get_molecule()
    mol = Molecule(qcmol=qcmol)
    surface_constraints = create_constraint_matrix(mol, esps)

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
