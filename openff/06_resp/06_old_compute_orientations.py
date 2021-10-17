#!/usr/bin/env python

from typing import List
import argparse
import pathlib
import itertools

import numpy as np
import qcelemental as qcel
import qcengine as qcng
from qcfractal import FractalSnowflakeHandler, FractalSnowflake, FractalServer
import qcelemental as qcel
import qcengine as qcng
from qcfractal import interface as ptl
from qcfractal.interface.models.records import RecordStatusEnum

from patch_qcfractal import TemporaryPostgres
import orutils
from info import PCM_KEYWORDS

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("infile", type=str)
parser.add_argument("--n_orientations", type=int, default=12)

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

def generate_orientations(qcmol, coordinates=None, n_orientations=4):
    combinations = generate_atom_combinations(qcmol, n_combinations=n_orientations)

    if coordinates is None:
        coordinates = qcmol.geometry

    all_coordinates = []
    for indices in combinations:
        all_coordinates.append(orutils.rigid_orient(*indices, coordinates))
    
    orientations = []
    for coords in all_coordinates:
        dct = qcmol.dict()
        dct["geometry"] = coords
        orientations.append(qcel.models.Molecule(**dct))
    return orientations, combinations



if __name__ == "__main__":
    args = parser.parse_args()

    path = pathlib.Path(args.infile)
    num, confid, state, jobname = path.stem.split("_")
    assert num == "02"
    assert jobname == "opt"

    result = qcel.models.AtomicResult.parse_file(args.infile)
    orientation_qcmols, combinations = generate_orientations(result.molecule,
                                               result.return_result,
                                               args.n_orientations)

    combinations = [list(x) for x in combinations]


    storage = TemporaryPostgres(database_name="test_psiresp")
    with FractalSnowflake(
        storage_project_name="test_psiresp",
        storage_uri=storage.psql.database_uri(),
        reset_database=True,
        start_server=True,
        logging=True,
    ) as server:
        print(server)
        client = ptl.FractalClient(server)
        print(client)

        keywords = {"maxiter": 300}
        if args.state != "gas":
            keywords.update(PCM_KEYWORDS)
        
        keyword_id = client.add_keywords([ptl.models.KeywordSet(values=keywords)])[0]
        print(keyword_id)

        wfn_protocols = {"wavefunction": "orbitals_and_eigenvalues"}

        computation = dict(
            program="psi4",
            basis=args.basis,
            method=args.method,
            driver="energy",
            molecule=orientation_qcmols,
            keywords=keyword_id,
            protocols=wfn_protocols,
        )
        print(computation)
        response = client.add_compute(**computation)

        print(response)
        ret = client.query_results(response.ids)
        print(ret)
        server.await_results()
        print("awaited")
        ret = client.query_results(response.ids)

    for i, record in enumerate(ret, 1):
        subfilename = f"03_c{confid}_o{i:02d}_{state}_energy.json"
        filename = path.parent / subfilename
        with open(filename, "w") as f:
            f.write(record.json())
        print(f"Wrote to {filename}")