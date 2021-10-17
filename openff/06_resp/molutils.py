import itertools

import numpy as np
import qcelemental as qcel

from orutils import rigid_orient

def _generate_atom_combinations(symbols):
    symbols = np.asarray(symbols)
    is_H = symbols == "H"
    h_atoms = list(np.flatnonzero(is_H))
    heavy_atoms = list(np.flatnonzero(~is_H) + 1)
    seen = set()

    for comb in itertools.combinations(heavy_atoms, 3):
        seen.add(comb)
        yield comb
        yield comb[::-1]

    for comb in itertools.combinations(heavy_atoms + h_atoms, 3):
        if comb in seen:
            continue
        seen.add(comb)
        yield comb
        yield comb[::-1]

def generate_atom_combinations(qcmol, n_combinations=None):
        atoms = _generate_atom_combinations(qcmol.symbols)
        if n_combinations is None or n_combinations < 0:
            return atoms

        return [next(atoms) for i in range(n_combinations)]


def generate_orientations(qcmol, n_reorientations=4):
    combinations = generate_atom_combinations(qcmol, n_reorientations)
    
    qcmols = []
    for ix in combinations:
        bohr_coordinates = rigid_orient(*ix, qcmol.geometry)
        mol_kwargs = qcmol.dict()
        mol_kwargs["geometry"] = bohr_coordinates
        qcmols.append(qcel.models.Molecule(**mol_kwargs))
    return qcmols