#!/usr/bin/env python

"""
This script converts raw MNSol data from the MNSolDatabase_v2012
to a usable CSV storing a PhysicalPropertyDataSet
"""


# inspired by https://github.com/openforcefield/openff-sage/blob/main/data-set-curation/physical-property/benchmarks/mnsol-to-csv.py

import argparse

parser = argparse.ArgumentParser("Convert MNSol data into dataset CSV")
parser.add_argument("--xls", type=str,
                    default="MNSolDatabase_v2012/MNSol_alldata.xls",
                    help="Input MNSol Excel file")
parser.add_argument("--json", type=str,
                    default="01_name_to_smiles.json",
                    help="JSON file for converting molecule names to SMILES")
parser.add_argument("--csv", type=str,
                    default="01_mnsol_data.csv",
                    help="Output CSV")
parser.add_argument("--verbose", action="store_true",
                    default=False, help="Toggle verbosity")

import json

from tqdm.auto import tqdm
import pandas as pd

from openff.evaluator import unit
from openff.evaluator.datasets import (MeasurementSource,
                                       PhysicalPropertyDataSet,
                                       PropertyPhase)
from openff.evaluator.properties import SolvationFreeEnergy
from openff.evaluator.substances import Component, ExactAmount, MoleFraction, Substance
from openff.evaluator.thermodynamics import ThermodynamicState


def read_excel(xls="MNSolDatabase_v2012/MNSol_alldata.xls",
               name_json="01_name_to_smiles.json",
               verbose=False):
    df = pd.read_excel(xls)
    if verbose:
        print(f"Found {len(df)} entries from {xls}.")

    # Load in a map between most of the common names found in the MNSol set and the
    # corresponding SMILES patterns generated using pubchem and NIH cactus.
    with open(name_json, "r") as f:
        name_to_smiles = json.load(f)
    
    # add new columns
    df["SoluteSMILES"] = df.SoluteName.apply(name_to_smiles.get)
    df["SolventSMILES"] = df.Solvent.apply(name_to_smiles.get)

    df2 = df[df.SoluteSMILES.notna() & df.SolventSMILES.notna()]
    if verbose:
        print(f"Found {len(df2)} valid entries.")
    return df2


def row_to_property(row):
    STATE = ThermodynamicState(temperature=298.0 * unit.kelvin,
                               pressure=1.0 * unit.atmosphere)
    KCAL_MOL = unit.kilocalorie / unit.mole
    SOURCE = MeasurementSource(doi="10.13020/3eks-j059")

    solvent = Component(smiles=row.SolventSMILES,
                            role=Component.Role.Solvent)
    solute = Component(smiles=row.SoluteSMILES,
                        role=Component.Role.Solute)

    substance = Substance()
    substance.add_component(solvent, MoleFraction(1.0))
    substance.add_component(solute, ExactAmount(1))

    sfe = SolvationFreeEnergy(thermodynamic_state=STATE,
                                phase=PropertyPhase.Liquid,
                                substance=substance,
                                value=row.DeltaGsolv * KCAL_MOL,
                                uncertainty=0.2 * KCAL_MOL,
                                source=SOURCE)
    return sfe


def create_mnsol_dataset(df, output="01_mnsol_data.csv", verbose=False):

    properties = []
    for row in tqdm(df.itertuples(), disable=not verbose):
        sfe = row_to_property(row)
        properties.append(sfe)

    openff_data_set = PhysicalPropertyDataSet()
    openff_data_set.add_properties(*properties, validate=False)
    openff_data_set.to_pandas().to_csv(output, index=False)
    if verbose:
        print(f"Wrote to output to {output}")


if __name__ == "__main__":
    args = parser.parse_args()
    df = read_excel(xls=args.xls, name_json=args.json, verbose=args.verbose)
    create_mnsol_dataset(df, output=args.csv, verbose=args.verbose)