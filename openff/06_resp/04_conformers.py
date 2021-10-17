#!/usr/bin/env python

import argparse
import os
import time
import json

from rdkit import Chem
import pandas as pd
import numpy as np
import qcelemental as qcel
import qcengine as qcng

from qcfractal import FractalSnowflakeHandler, FractalSnowflake, FractalServer
import qcelemental as qcel
import qcengine as qcng
from qcfractal import interface as ptl
from qcfractal.interface.models.records import RecordStatusEnum

from info import PCM_KEYWORDS
from rdutils import generate_diverse_conformer_ids
from molutils import generate_orientations
from waiter import wait_for_procedures

ANGSTROM_TO_BOHR = qcel.constants.conversion_factor("angstrom", "bohr")
PROTOCOLS = {"wavefunction": "orbitals_and_eigenvalues"}

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("entry_index", type=int)
parser.add_argument("--component", type=str, default="solvent")
parser.add_argument("--server", default="localhost")
parser.add_argument("--port", default="7777")

def qcmol_from_rdkit(rdmol, confId=-1):
    from rdkit import Chem

    symbols = [a.GetSymbol() for a in rdmol.GetAtoms()]
    if rdmol.GetNumConformers() == 0:
        validate = False
        geometry = np.zeros((len(symbols), 3))
    else:
        validate = None
        geometry = np.array(rdmol.GetConformer(confId).GetPositions())

    geometry *= ANGSTROM_TO_BOHR

    connectivity = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
        for b in rdmol.GetBonds()
    ]

    qcmol = qcel.models.Molecule(symbols=symbols,
                                 geometry=geometry,
                                 validate=validate,
                                 connectivity=connectivity,
                                 molecular_charge=Chem.GetFormalCharge(rdmol),
                                 )
    return qcmol





def get_smiles_and_molname(entry_index=0, component="solvent"):
    df = pd.read_csv("../05_properties/01_mnsol_data.csv")

    row = df.iloc[entry_index]
    if row["Role 1"].lower() == component:
        smiles = row["Component 1"]
    else:
        smiles = row["Component 2"]
    return smiles, f"mol-{entry_index:04d}_{component}"

def create_conformers(smiles, molname):
    smiles_parser = Chem.rdmolfiles.SmilesParserParams()
    smiles_parser.removeHs = False
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ids = generate_diverse_conformer_ids(rdmol)
    conformers = [qcmol_from_rdkit(rdmol, confId=i) for i in ids]
    print(f"Generated {len(conformers)} conformers")
    for i, qcmol in enumerate(conformers, 1):
        with open(f"{molname}/01_conf-{i:02d}.json", "w") as f:
            f.write(qcmol.json())
    return conformers

def run_geometry_optimization(client, qcmols, state="gas"):

    keywords = {"maxiter": 300,
                "geom_maxiter": 3000,}

    if state == "water":
        keywords.update(PCM_KEYWORDS)

    keyword_id = client.add_keywords([ptl.models.KeywordSet(values=keywords)])[0]

    spec = {
        "keywords": None,
        "qc_spec": {
            "driver": "gradient",
            "method": "PW6B95",
            "basis": "cc-pV(D+d)Z",
            "program": "psi4",
            "protocols": PROTOCOLS,
            "keywords": keyword_id,
        },
    }

    # Ask the server to compute a new computation
    response = client.add_procedure("optimization", "geometric", spec, qcmols)
    return response.ids

def get_final_molecule(client, response_ids):
    results = wait_for_procedures(client, response_ids=response_ids)

    results = sorted(results, key=lambda x: int(x.id))
    
    return [r.get_final_molecule() for r in results]



if __name__ == "__main__":
    args = parser.parse_args()

    smiles, molname = get_smiles_and_molname(entry_index=args.entry_index,
                                             component=args.component.lower())

    os.makedirs(molname, exist_ok=True)

    client = ptl.FractalClient(address=f"{args.server}:{args.port}",
                               verify=False)
    conformer_qcmols = create_conformers(smiles, molname)

    # response_ids = {}
    # # for state in ["gas", "water"]:
    # for state in ["water"]:
    #     ids = run_geometry_optimization(client, conformer_qcmols, state=state)
    #     response_ids[state] = ids

    # with open(f"{molname}/02_optimization_ids.json", "w") as f:
    #     json.dump(response_ids, f)


    # for state, ids in response_ids.items():
    #     optmols = get_final_molecule(client, ids)
    #     for i, qcmol in enumerate(optmols, 1):
    #         filename = f"02_c{i:02d}_{state}_opt.json"
    #         path = f"{molname}/{filename}"
    #         with open(path, "w") as f:
    #             f.write(qcmol.json())
    #         print(f"Wrote to {path}")

    

