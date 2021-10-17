#!/usr/bin/env python

from typing import List
import argparse
import pathlib
import itertools
import json
import time

import numpy as np
import qcelemental as qcel
import qcengine as qcng
from qcfractal import FractalSnowflakeHandler, FractalSnowflake, FractalServer
import qcelemental as qcel
import qcengine as qcng
from qcfractal import interface as ptl
from qcfractal.interface.models.records import RecordStatusEnum
from qcfractal.interface.models.records import ResultRecord, OptimizationRecord

from patch_qcfractal import TemporaryPostgres
import orutils
from info import PCM_KEYWORDS, PROTOCOLS

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("infile", type=str)
parser.add_argument("--server", default="localhost")
parser.add_argument("--port", default="7777")
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
        xyz = coordinates.copy()
        all_coordinates.append(orutils.rigid_orient(*indices, coordinates))
    
    orientations = []
    for coords in all_coordinates:
        # print('COORDS', coords)
        dct = qcmol.dict()
        dct["geometry"] = coords #* qcel.constants.conversion_factor("angstrom", "bohr")
        orientations.append(qcel.models.Molecule(**dct))
    return orientations, combinations


def wait_for_procedures(client, response_ids=[], query_interval=20):
    n_incomplete = len(response_ids)
    while(n_incomplete):
        time.sleep(query_interval)
        results = client.query_procedures(id=response_ids)
        status = [r.status for r in results]
        status = np.array([s.value
                            if isinstance(s, RecordStatusEnum)
                            else s
                            for s in status])
        n_incomplete = (status == "INCOMPLETE").sum()
    results = client.query_procedures(id=response_ids)

    for record in results:
        if record.status == "ERROR":
            print("=== ERROR ===")
            print(record.get_error().error_message)
            print("")
    return results

def wait_for_results(client, response_ids=[], query_interval=20):
    n_incomplete = len(response_ids)
    while(n_incomplete):
        time.sleep(query_interval)
        results = client.query_results(id=response_ids)
        status = [r.status for r in results]
        status = np.array([s.value
                            if isinstance(s, RecordStatusEnum)
                            else s
                            for s in status])
        n_incomplete = (status == "INCOMPLETE").sum()
    results = client.query_results(id=response_ids)

    for record in results:
        if record.status == "ERROR":
            print("=== ERROR ===")
            print(record.get_error().error_message)
            print("")
    return results

def run_geometry_optimization(client, qcmols):
    spec = {
        "keywords": None,
        "qc_spec": {
            "driver": "gradient",
            "method": "PW6B95",
            "basis": "cc-pV(D+d)Z",
            "program": "psi4",
            "protocols": PROTOCOLS
        },
    }

    # Ask the server to compute a new computation
    response = client.add_procedure("optimization", "geometric", spec, qcmols)
    print(response)
    results = wait_for_procedures(client, response_ids=response.ids)
    return results
    


if __name__ == "__main__":
    args = parser.parse_args()

    client = ptl.FractalClient(address=f"{args.server}:{args.port}",
                               verify=False)
    print(f"client: {client}")




    path = pathlib.Path(args.infile)
    confid = path.stem.split("-")[1].split(".")[0]

    optfile = path.parent / f"02_c{confid}_opt.json"

    # try:
    #     result = OptimizationRecord.parse_file(optfile)
    # except:
    mol = qcel.models.Molecule.from_file(args.infile)
    result = run_geometry_optimization(client, [mol])[0]
    with open(str(optfile), "w") as f:
        f.write(result.get_final_molecule().json())
    print(f"Wrote to {optfile}")

    optmol = result.get_final_molecule()
    

    orientation_qcmols, combinations = generate_orientations(optmol,
                                                            optmol.geometry,
                                                            args.n_orientations)

    combinations = [[int(y) for y in x] for x in combinations]
    with open(f"{path.parent}/03_c{confid}_combinations.json", "w") as f:
        json.dump(combinations, f)
    print(f"Wrote to {path.parent}/03_c{confid}_combinations.json")

    for i, qcmol in enumerate(orientation_qcmols, 1):
        fname = path.parent / f"03_c{confid}_o{i:02d}_mol.json"
        with open(fname, "w") as f:
            f.write(qcmol.json())


    response_file = path.parent / f"03_c{confid}_response_ids.json"
    response_ids = {}
    for state in ["gas", "water"]:
        keywords = {"maxiter": 300}
        if state != "gas":
            keywords.update(PCM_KEYWORDS)
        print(keywords)
        keyword_id = client.add_keywords([ptl.models.KeywordSet(values=keywords)])[0]
        print(f"keyword_id: {keyword_id}")

        computation = dict(
            program="psi4",
            basis=result.qc_spec.basis,
            method=result.qc_spec.method,
            driver="energy",
            molecule=orientation_qcmols,
            keywords=keyword_id,
            protocols=PROTOCOLS,
        )

        response = client.add_compute(**computation)
        response_ids[state] = response.ids

    with open(response_file, "w") as f:
        json.dump(response_ids, f)
    print(f"Wrote to {response_file}")

    for state in ["gas", "water"]:
        results = wait_for_results(client, response.ids)

        for i, record in enumerate(results, 1):
            subfilename = f"03_c{confid}_o{i:02d}_{state}_energy.json"
            filename = path.parent / subfilename
            with open(filename, "w") as f:
                f.write(record.json())
            print(f"Wrote to {filename}")

    