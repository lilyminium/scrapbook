#!/usr/bin/env python


import pathlib
import json
import glob
import os
import argparse
import tqdm

import pandas as pd
import numpy as np
import qcelemental as qcel
from rdkit import Chem
from numpy.testing import assert_equal

# from openff.toolkit.typing.engines.smirnoff import ForceField
# from openff.evaluator.backends import ComputeResources
# from openff.evaluator.backends.dask import DaskLocalCluster
# from openff.toolkit.topology.molecule import Molecule, unit
# from openff.evaluator.workflow import Workflow


# from utils import dataset_from_csv


parser = argparse.ArgumentParser("Evaluate properties")
parser.add_argument("elements", nargs='+', type=str)

def get_entries_with_elements(elements):
    df = pd.read_csv("01_mnsol_data.csv")
    entries = []
    for i, row in df.iterrows():
        rd1 = Chem.MolFromSmiles(row["Component 1"])
        rd2 = Chem.MolFromSmiles(row["Component 2"])
        el1 = {at.GetSymbol() for at in rd1.GetAtoms()}
        el2 = {at.GetSymbol() for at in rd2.GetAtoms()}
        if el1.issubset(elements) and el2.issubset(elements):
            entries.append(i)
    return entries



if __name__ == "__main__":
    args = parser.parse_args()
    elements = set(args.elements)
    entries = get_entries_with_elements(elements)
    print(len(entries))
    print(entries)