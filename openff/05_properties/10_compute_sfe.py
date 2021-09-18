"""
Ok, working backwards to run jobs.
"""

import argparse

parser = argparse.ArgumentParser("Compute SFE")
parser.add_argument("--index", type=int, default=0,
                    help="Index of entry")
parser.add_argument("--csv", type=str,
                    default="01_mnsol_data.csv",
                    help="Input CSV")
parser.add_argument("--verbose", type=bool, action="store_true",
                    default=False, help="Toggle verbosity")
parser.add_argument("--dirname", type=str,
                    default="data", help="directory prefix")



import pandas as pd

from openff.evaluator.datasets import (
    PhysicalPropertyDataSet,
)

from openff.evaluator.client import RequestOptions
from openff.evaluator import unit
from openff.evaluator.storage import LocalFileStorage
from openff.evaluator.workflow import Workflow


def entry_from_csv(csv, entry=0):
    df = pd.read_csv(csv)
    dataset = PhysicalPropertyDataSet.from_pandas(df)
    return dataset.properties[entry]

def estimate_entry(entry, working_directory=""):
    metadata = {}

    property_type = type(entry)
    property_schema = property_type.default_simulation_schema()
    workflow = Workflow.from_schema(property_schema, metadata=metadata)
    result = workflow.execute(root_directory=working_directory)



if __name__ == "__main__":
    args = parser.parse_args()
    entry = entry_from_csv(args.csv, entry=args.index)