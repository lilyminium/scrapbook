#!/usr/bin/env python


import pandas as pd

from openff.evaluator.datasets import (
    PhysicalPropertyDataSet,
)

from openff.evaluator.client import RequestOptions
from openff.evaluator import unit
from openff.evaluator.storage import LocalFileStorage
from openff.evaluator.workflow import Workflow

def dataset_from_csv(csv="01_mnsol_data.csv"):
    df = pd.read_csv(csv)
    df = df.iloc[:1]
    dataset = PhysicalPropertyDataSet.from_pandas(df)
    
    return dataset

def entry_from_csv(csv="01_mnsol_data.csv", entry=0):
    df = pd.read_csv(csv)
    dataset = PhysicalPropertyDataSet.from_pandas(df)
    return dataset.properties[entry]

def estimate_entry(entry, working_directory=""):
    metadata = {"substance": entry.substance,
                "force_field_path": "openff-1.0.0.offxml",
                "thermodynamic_state": entry.thermodynamic_state,
                "parameter_gradient_keys": [],
               }

    property_type = type(entry)
    property_schema = property_type.default_simulation_schema()
    workflow = Workflow.from_schema(property_schema.workflow_schema, metadata=metadata)
    return workflow
    result = workflow.execute(root_directory=working_directory)

entry = entry_from_csv()


from openff.evaluator.properties import Density, EnthalpyOfVaporization, SolvationFreeEnergy
from openff.evaluator.client import RequestOptions

# Create an options object which defines how the data set should be estimated.
estimation_options = RequestOptions()
# Specify that we only wish to use molecular simulation to estimate the data set.
estimation_options.calculation_layers = ["SimulationLayer"]

estimation_options.add_schema("SimulationLayer", "SolvationFreeEnergy", SolvationFreeEnergy.default_simulation_schema())

from openff.evaluator.backends import ComputeResources
from openff.evaluator.backends.dask import DaskLocalCluster

calculation_backend = DaskLocalCluster(
    number_of_workers=1,
    resources_per_worker=ComputeResources(
        number_of_threads=1,
#         number_of_gpus=1,
#         preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA
    ),
)
calculation_backend.start()


from openff.evaluator.server import EvaluatorServer

evaluator_server = EvaluatorServer(calculation_backend=calculation_backend)
evaluator_server.start(asynchronous=True)

data_set = dataset_from_csv()

from openff.evaluator.forcefield import SmirnoffForceFieldSource

force_field_path = "openff-1.3.0.offxml"
force_field_source = SmirnoffForceFieldSource.from_path(force_field_path)

from openff.evaluator.client import EvaluatorClient
evaluator_client = EvaluatorClient()

request, exception = evaluator_client.request_estimate(
    property_set=data_set,
    force_field_source=force_field_source,
    options=estimation_options,
)

assert exception is None

results, exception = request.results(synchronous=True, polling_interval=30)
assert exception is None

print(results)

print(dir(results))

for x in results.__dict__:
    if "uncertainty" in x:
        print(x, getattr(results, x))