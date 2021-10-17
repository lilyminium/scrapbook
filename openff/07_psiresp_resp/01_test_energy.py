#!/usr/bin/env python

import argparse
import pathlib

import qcelemental as qcel
import qcengine as qcng

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("infile", type=str)
parser.add_argument("--basis", type=str, default="6-31g*")
parser.add_argument("--method", type=str, default="hf")
parser.add_argument("--state", type=str, default="gas")

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


# Update database
# FractalServer.storage = SQLAlchemySocket
# FractalServer.storage.add_results([ResultRecords])

if __name__ == "__main__":
    args = parser.parse_args()
    mol = qcel.models.Molecule.from_file(args.infile)

    computation = {
        "molecule": mol,
        "driver": "energy",
        "model": {"method": args.method, "basis": args.basis},
        "protocols": {"wavefunction": "orbitals_and_eigenvalues"},
    }

    if args.state != "gas":
        computation["keywords"] = PCM_KEYWORDS

    ret = qcng.compute(computation, "psi4")

    molname = pathlib.Path(args.infile).stem

    if args.method.upper() == "PW6B95":
        suffix = "_resp2"
    else:
        suffix = ""

    jsonfile = f"output/{molname}{suffix}_{args.state}_energy.json"
    with open(jsonfile, "w") as f:
        f.write(ret.json())