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
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import AmberToolsToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper

parser = argparse.ArgumentParser(description='Check charge range of molecule',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("nmol", type=int, help="molecule number")
parser.add_argument("norder", type=int, help="order number")

PARSLEY_BENCHMARK_SMILES = "/DFS-L/SCRATCH/mobley/lilyw7/charges/parsley_benchmark/openff_parsley_benchmark.smi"

logger = logging.getLogger(__name__)
logger.setLevel(logging.CRITICAL)  # not super interested in OpenFF's logging

def get_n_samples(n_heavy):
    n_samples = -1
    if n_heavy > 7:
        n_samples = 1_000_000
    # if n_heavy > 10:
    #     n_samples = 200
    # if n_heavy > 20:
    #     n_samples = 50
    # if n_heavy > 40:
    #     n_samples = 20
    # if n_heavy > 60:
    #     n_samples = 5
    return n_samples

def get_permutations(n_mol, mol):
    mapping = {i: i for i in range(mol.n_atoms)}
    heavy_indices = [i for i, atom in enumerate(mol.atoms)
                     if atom.atomic_number != 1]
    
    n_samples = get_n_samples(len(heavy_indices))
    print("n samples", n_samples)
    if n_samples < 0:
        permutations = itertools.permutations(heavy_indices)
    else:
        permutations = set()
        while len(permutations) < n_samples:
            permutations.add(tuple(np.random.permutation(heavy_indices)))

    remappings = []
    unmappings = []
    for permut in permutations:
        remapping = dict(mapping)
        # the below doesn't work becuase they don't like int64s
        # remapping.update(dict(zip(heavy_indices, permut)))
        remapping.update({int(k): int(v) for k, v in zip(heavy_indices, permut)})
        remappings.append(remapping)
        # unmappings.append({v: k for k, v in remapping.items()})
    

    return remappings
    
    

def get_conformers(i, j):

    # with open(PARSLEY_BENCHMARK_SMILES, "r") as f:
    #     contents = f.readlines()
    # smiles = contents[i - 1].strip()

    # # I believe this is canonical
    # mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    # remappings = get_permutations(mol)

    MAX_N_CONFS = 10000
    MOLNAME = f"conformers/mol{i:04d}"
    ALL_RDCONFS = []
    ALL_OECONFS = []

    mapname = f"remappings/mol{i:04d}.smi"
    with open(mapname, "r") as f:
        mapped_smiles = [x.strip() for x in f.readlines()]
    
    offmol = Molecule.from_mapped_smiles(mapped_smiles[j - 1])

    # for j, remapping in tqdm.tqdm(enumerate(remappings, 1)):
    ORDERNAME = f"{MOLNAME}_m{j:05d}"

    print("order name", ORDERNAME)

    if not os.path.isfile(f"{ORDERNAME}_rd.npy"):

    # rdkit
        rdmol = offmol.to_rdkit()
        confs = AllChem.EmbedMultipleConfs(rdmol, numConfs=MAX_N_CONFS,
                                            pruneRmsThresh=1.0,
                                            randomSeed=1)
        rdconfs = np.array([conf.GetPositions() for conf in rdmol.GetConformers()])
        np.save(f"{ORDERNAME}_rd.npy", rdconfs)

    if not os.path.isfile(f"{ORDERNAME}_oe.npy"):

        # openeye
        oemol = offmol.to_openeye()
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(MAX_N_CONFS)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)
        omega.SetEnergyWindow(15.0)  # unit?
        omega.SetRMSThreshold(1.0)
        omega.SetStrictStereo(True)
        status = omega(oemol)

        if status is False:
            omega.SetStrictStereo(False)
            new_status = omega(oemol)
        
        oeconfs = []
        for conf in oemol.GetConfs():
            coord_dict = conf.GetCoords() # wtf this is a dictionary
            position_array = []
            for i in sorted(coord_dict):
                position_array.append(coord_dict[i])
            position_array = np.array(position_array)
            if not (position_array == 0).all():
                oeconfs.append(position_array)
        oeconfs = np.array(oeconfs)
        np.save(f"{ORDERNAME}_oe.npy", oeconfs)
        #     ALL_OECONFS.append(oeconfs)
        #     ALL_RDCONFS.append(rdconfs)
        # return ALL_RDCONFS, ALL_OECONFS


if __name__ == "__main__":
    args = parser.parse_args()
    i = args.nmol
    j = args.norder
    get_conformers(i, j)
