import argparse
import json
import time
import os

from rdkit import Chem
import pandas as pd
import numpy as np
import qcelemental as qcel
from qcfractal import interface as ptl

from qcfractal.interface.models.records import RecordStatusEnum
from qcfractal import FractalSnowflakeHandler

from rdutils import generate_diverse_conformer_ids
from molutils import generate_orientations

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("entry_index", type=int)
parser.add_argument("--component", type=str, default="solvent")
parser.add_argument("--server", default="localhost")
parser.add_argument("--port", default="7777")

ANGSTROM_TO_BOHR = qcel.constants.conversion_factor("angstrom", "bohr")
PROTOCOLS = {"wavefunction": "orbitals_and_eigenvalues"}

PCM_KEYWORDS = {
    "pcm": "true",
    "pcm_scf_type": "total",
    "pcm__input": r"""
        Units = Angstrom
        Medium {
            SolverType = CPCM
            Solvent = water
        }

        Cavity {
            RadiiSet = Bondi # Bondi | UFF | Allinger
            Type = GePol
            Scaling = True # radii for spheres scaled by 1.2
            Area = 0.3
            Mode = Implicit
        }
        """
}

server = FractalSnowflakeHandler()
print(server)


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


def wait_for_results(client, response_ids=[], query_interval=20):
    # server.await_results()
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
    print(results[0].get_error())
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
    results = wait_for_results(client, response_ids=response.ids)
    
    return [r.get_final_molecule() for r in results]

def optimize_geometry(client, qcmols, molname):
    optimized = run_geometry_optimization(client, qcmols)
    for i, qcmol in enumerate(optimized, 1):
        with open(f"{molname}/02_optconf-{i:02d}.json", "w") as f:
            f.write(qcmol.json())
    return optimized

def get_smiles_and_molname(entry_index=0, component="solvent"):
    df = pd.read_csv("../05_properties/01_mnsol_data.csv")

    row = df.iloc[entry_index]
    if row["Role 1"].lower() == component:
        smiles = row["Component 1"]
    else:
        smiles = row["Component 2"]
    return smiles, f"mol-{entry_index:04d}_{component}"

def submit_energy(client, qcmols, keyword_id=None):
    response = client.add_compute(program="psi4",
                                  driver="energy",
                                  method="PW6B95",
                                  basis="cc-pV(D+d)Z",
                                  protocols=PROTOCOLS,
                                  keywords=keyword_id)
    return response.ids

def create_and_write_orientations(qcmol, conf_id=1):
    orientations = generate_orientations(qcmol)
    for i, orient in enumerate(orientations, 1):
        with open(f"{molname}/03_conf-{conf_id:02d}_orient-{i}.json", "w") as f:
            f.write(orient.json())
    return orientations


def run_resp_calculations(client, smiles, molname, keyword_id):
    conformer_qcmols = create_conformers(smiles, molname)
    optimized_qcmols = optimize_geometry(client, conformer_qcmols, molname)

    energy_calcs = []

    for i, qcmol in enumerate(optimized_qcmols, 1):
        orientation_qcmols = create_and_write_orientations(qcmol)
        gas_ids = submit_energy(client, orientation_qcmols)
        water_ids = submit_energy(client, orientation_qcmols,
                                  keyword_id=keyword_id)
        for j, orient in enumerate(orientation_qcmols, 1):
            entry = {
                "conformer": i,
                "orientation": j,
                "gas_id": gas_ids[j - 1],
                "water_id": water_ids[j - 1],
                "hash": orient.get_hash(),
                "qcmol_str": orient.json(),
            }
            energy_calcs.append(entry)
    
    with open(f"{molname}/energy_calculations.json", "w") as f:
        json.dump(energy_calcs, f)




if __name__ == "__main__":
    args = parser.parse_args()
    

    # client = ptl.FractalClient(address=f"{args.server}:{args.port}",
    #                        verify=False)

    client = server.client()
            
    keyword_set = ptl.models.KeywordSet(values=PCM_KEYWORDS)
    keyword_id = client.add_keywords([keyword_set])[0]

    smiles, molname = get_smiles_and_molname(entry_index=args.entry_index,
                                            component=args.component.lower())

    os.makedirs(molname, exist_ok=True)
    run_resp_calculations(client, smiles, molname, keyword_id)


    
