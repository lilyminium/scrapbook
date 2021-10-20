
from typing import List
import itertools

import numpy as np
import qcelemental as qcel

import orutils

def _generate_atom_combinations(symbols: List[str]):
    """Yield combinations of atom indices for transformations

    The method first yields combinations of 3 heavy atom indices.
    Each combination is followed by its reverse. Once the heavy atoms
    are exhausted, the heavy atoms then get combined with the hydrogens.

    Parameters
    ----------
    symbols: list of str
        List of atom elements

    Examples
    --------

    ::

        >>> symbols = ["H", "C", "C", "O", "N"]
        >>> comb = generate_atom_combinations(symbols)
        >>> next(comb)
        (1, 2, 3)
        >>> next(comb)
        (3, 2, 1)
        >>> next(comb)
        (1, 2, 4)
        >>> next(comb)
        (4, 2, 1)

    """
    symbols = np.asarray(symbols)
    is_H = symbols == "H"
    h_atoms = list(np.flatnonzero(is_H))
    heavy_atoms = list(np.flatnonzero(~is_H))
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

    combinations = []
    for i in range(n_combinations):
        try:
            combinations.append(next(atoms))
        except:
            pass
    return combinations

def generate_orientations(qcmol, combinations):
    coordinates = qcmol.geometry 

    all_coordinates = []
    for indices in combinations:
        xyz = coordinates.copy()
        all_coordinates.append(orutils.rigid_orient(*indices, coordinates))
    
    orientations = []
    for coords in all_coordinates:
        # print('COORDS', coords)
        dct = qcmol.dict()
        dct["geometry"] = coords #* qcel.constants.conversion_factor("angstrom", "bohr")
        orientations.append(qcel.models.Molecule(**dct))
    return orientations