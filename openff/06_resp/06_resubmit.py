#!/usr/bin/env python

from typing import List
import argparse
import pathlib
import itertools
import json
import time
import glob

import numpy as np
from qcfractal import interface as ptl
from qcfractal import FractalSnowflakeHandler, FractalSnowflake, FractalServer
import qcelemental as qcel
import qcengine as qcng
from qcfractal.interface.models.records import RecordStatusEnum


from psiresp.grid import GridOptions
from psiresp.orientation import OrientationEsp
from psiresp.constraint import ConstraintMatrix
from psiresp.resp import RespCharges, RespOptions
from psiresp.charge import MoleculeChargeConstraints, ChargeConstraintOptions
from psiresp.molecule import Molecule
from psiresp import psi4utils

parser = argparse.ArgumentParser("Run resp")
parser.add_argument("--resubmit", action="store_true", default=False)
parser.add_argument("--pattern", type=str, default="mol*/02_optimization_ids.json")


if __name__ == "__main__":
    args = parser.parse_args()
    client = ptl.FractalClient(address="hpc3-l18-05:7777", verify=False)

    response_files = sorted(glob.glob(args.pattern))
    conformer_response_ids = []
    for file in response_files:
        with open(file, "r") as f:
            response_ids = json.load(f)
            for k, v in response_ids.items():
                conformer_response_ids.extend(v)
    
    response_ids = sorted(conformer_response_ids)
    print(response_ids)

    results = client.query_results(id=response_ids)

    if not results:
        results = client.query_procedures(id=response_ids)
    status = [s.status for s in results]
    status = np.array([s.value
                            if isinstance(s, RecordStatusEnum)
                            else s
                            for s in status])

    print(status)
    ids = np.array([r.id for r in results])

    n_status = len(results)
    print(f"Complete: {sum(status == 'COMPLETE')} / {n_status}")
    print(f"Errored: {sum(status == 'ERROR')} / {n_status}")
    print(f"Incomplete: {sum(status == 'INCOMPLETE')} / {n_status}")

    ix = np.where(status == "ERROR")[0]
    # print(ix)

    if len(ix):
        record = results[ix[0]]
        # print(record)
        print(record.get_error().error_message)

    if args.resubmit:
        n = client.modify_tasks("restart", list(ids[status == "ERROR"]))
        print(f"Updated: {n.n_updated}")
