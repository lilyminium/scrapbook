import glob
import numpy as np
from MDAnalysis.analysis import rms

PARSLEY_BENCHMARK_SMILES = "/DFS-L/SCRATCH/mobley/lilyw7/charges/parsley_benchmark/04_canonical_conformers/openff_parsley_benchmark_mapped.smi"


def get_adjacency_matrix(molecule):
    matrix = np.zeros((molecule.n_atoms, molecule.n_atoms), dtype=int)
    for bond in molecule.bonds:
        i, j = bond.atom1_index, bond.atom2_index
        matrix[i, j] = matrix[j, i] = 1
    return matrix

def read_canonical_mol(i=1):
    from openff.toolkit.topology import Molecule
    with open(PARSLEY_BENCHMARK_SMILES, "r") as f:
        smiles = [x.strip() for x in f.readlines()]
    
    smiles = smiles[i - 1]
    offmol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    return offmol

