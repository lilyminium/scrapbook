#!/usr/bin/env python

import itertools
import pandas as pd
import numpy as np
import tqdm
import argparse
import os
import sys
import subprocess
import json
import logging
from openeye import oeomega

from rdkit.Chem import AllChem


from openff.toolkit.topology import Molecule
from openff.toolkit.topology.molecule import unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import AmberToolsToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper

parser = argparse.ArgumentParser(description='Check charge range of molecule',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("nmol", type=int, help="molecule number")

PARSLEY_BENCHMARK_SMILES = "/DFS-L/SCRATCH/mobley/lilyw7/charges/parsley_benchmark/04_canonical_conformers/openff_parsley_benchmark_mapped.smi"
MAX_N_CONFS = 500

logger = logging.getLogger(__name__)
logger.setLevel(logging.CRITICAL)  # not super interested in OpenFF's logging

    
    

def get_conformers(i):
    import json

    filename = f"canonical_conformers/mol{i:04d}_confs.npy"

    if not os.path.isfile(filename):

        rdwrapper = RDKitToolkitWrapper()
        oewrapper = OpenEyeToolkitWrapper()
        atwrapper = AmberToolsToolkitWrapper()

        with open(PARSLEY_BENCHMARK_SMILES, "r") as f:
            contents = f.readlines()
        smiles = contents[i - 1].strip()

        rdoffmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        rdoffmol.generate_conformers(toolkit_registry=rdwrapper,
                                    n_conformers=MAX_N_CONFS)
        rdconfs = np.array([c.value_in_unit(unit.angstrom) for c in rdoffmol.conformers])

        oeoffmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        oeoffmol.generate_conformers(toolkit_registry=oewrapper,
                                    n_conformers=MAX_N_CONFS)
        oeconfs = np.array([c.value_in_unit(unit.angstrom) for c in oeoffmol.conformers])

        confs = np.concatenate([rdconfs, oeconfs])
        names = np.array(["RDKit"] * len(rdconfs) + ["OpenEye"] * len(oeconfs))

        np.save(f"canonical_toolkits/mol{i:04d}_tks.npy", names)
        np.save(filename, confs)

if __name__ == "__main__":
    args = parser.parse_args()
    i = args.nmol

    get_conformers(i)
