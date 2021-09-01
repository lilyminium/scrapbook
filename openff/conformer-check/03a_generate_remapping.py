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

PARSLEY_BENCHMARK_SMILES = "/DFS-L/SCRATCH/mobley/lilyw7/charges/parsley_benchmark/openff_parsley_benchmark.smi"
MAX_REORDERINGS = 10000

logger = logging.getLogger(__name__)
logger.setLevel(logging.CRITICAL)  # not super interested in OpenFF's logging

def get_n_samples(n_heavy):
    n_samples = -1
    if n_heavy > 7:
        n_samples = MAX_REORDERINGS
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

    print(len(heavy_indices))
    
    n_samples = get_n_samples(len(heavy_indices))
    print("n samples", n_samples)
    if n_samples < 0:
        permutations = itertools.permutations(heavy_indices)
    else:
        permutations = set()
        n_iter = MAX_REORDERINGS * 2
        while len(permutations) < n_samples and n_iter:
            permutations.add(tuple(np.random.permutation(heavy_indices)))
            n_iter -= 1
        all_permut = itertools.permutations(heavy_indices)
        while len(permutations) < n_samples:
            permutations.add(next(all_permut))

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
    
    

def get_conformers(i):
    import os
    outfile = f"remappings/mol{i:04d}.smi"
    if os.path.isfile(outfile):
        import sys
        sys.exit()

    with open(PARSLEY_BENCHMARK_SMILES, "r") as f:
        contents = f.readlines()
    smiles = contents[i - 1].strip()

    # I believe this is canonical
    mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    remappings = get_permutations(i, mol)

    all_smiles = []

    for remapping in tqdm.tqdm(remappings):
        unmapping = {v: k for k, v in remapping.items()}
        offmol = mol.remap(remapping)
        msmi = offmol.to_smiles(mapped=True)
        all_smiles.append(msmi)

    with open(outfile, "w") as f:
        f.write("\n".join(all_smiles))




if __name__ == "__main__":
    args = parser.parse_args()
    i = args.nmol

    get_conformers(i)
