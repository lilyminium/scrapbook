import time

import numpy as np

from qcfractal import FractalSnowflakeHandler, FractalSnowflake, FractalServer
import qcelemental as qcel
import qcengine as qcng
from qcfractal import interface as ptl
from qcfractal.interface.models.records import RecordStatusEnum

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

    response_ids = [x for x in response_ids if x is not None]
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