#!/usr/bin/env python

import argparse

import qcelemental as qcel
import qcengine as qcng

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("mol_index", type=int, default=10)
parser.add_argument("conf_id", type=int, default=1)
parser.add_argument("--component", type=str, default="solvent")
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

if __name__ == "__main__":
    args = parser.parse_args()

    molname = f"mol-{args.mol_index:04d}_{args.component}"
    confname = f"01_conf-{args.conf_id:02d}"
    molfile = f"{molname}/{confname}.json"
    mol = qcel.models.Molecule.from_file(molfile)

    print(mol)

    computation = {
        "molecule": mol,
        "driver": "gradient",
        "model": {"method": "PW6B95", "basis": "cc-pV(D+d)Z"},
        "protocols": {"wavefunction": "orbitals_and_eigenvalues"},
    }

    if args.state != "gas":
        computation["keywords"] = PCM_KEYWORDS

    print(computation)

    ret = qcng.compute(computation, "psi4")

    jsonfile = f"{molname}/02_{args.conf_id:02d}_{args.state}_opt.json"
    with open(jsonfile, "w") as f:
        f.write(ret.json())

