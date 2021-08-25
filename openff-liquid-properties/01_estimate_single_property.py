#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import os
from openff.evaluator.datasets.thermoml import ThermoMLDataSet

from openff.evaluator.backends import ComputeResources
from openff.evaluator.backends.dask import DaskLocalCluster
from openff.evaluator.server import EvaluatorServer
from openff.evaluator.client import EvaluatorClient
from openff.evaluator.forcefield import SmirnoffForceFieldSource

from openff.evaluator.client import RequestOptions
from openff.evaluator import unit
from openff.evaluator.storage import LocalFileStorage
from openff.evaluator.properties import (Density,
                                         EnthalpyOfVaporization,
                                         ExcessMolarVolume,
                                         DielectricConstant,
                                         EnthalpyOfMixing)

force_field_path = "openff_unconstrained-1.3.0.offxml"
force_field_source = SmirnoffForceFieldSource.from_path(force_field_path)


if "OE_LICENSE" in os.environ:
    tk = "oe"
else:
    tk = "at"

if os.environ["ELF10"] == "false":
    method = "am1bcc"
else:
    method = "am1bccelf10"

BIBFILE = "data/00_physprop-benchmark-sources.bib"
RAW_JSON = "data/01_raw_thermoml_dataset_from_sources.json"
RAW_CSV = "data/01_raw_thermoml_dataset_from_sources.csv"
ENTRY_DIRNAME = "data/01_raw_data_entries"
ESTIMATED_DIRNAME = f"data/02_raw_data_estimated_{tk}_{method}"
WORKING = f"working/{tk}_{method}_working-data"
STORED = f"working/{tk}_{method}_stored-data"

os.makedirs(ESTIMATED_DIRNAME, exist_ok=True)
os.makedirs(WORKING, exist_ok=True)
os.makedirs(STORED, exist_ok=True)

ENTRY_NAME = os.path.join(ENTRY_DIRNAME, "entry_{:05d}.csv")
JSON_NAME = os.path.join(ESTIMATED_DIRNAME, "entry_{:05d}.json")

ABSOLUTE_TOLERANCES = {
    Density: 0.45 * unit.kilogram * unit.meter ** -3,
    DielectricConstant: 1.5 * unit.dimensionless,
    EnthalpyOfVaporization: 0.65 * unit.kilojoule / unit.mole,
    EnthalpyOfMixing: 0.02 * unit.kilojoule / unit.mole,
    ExcessMolarVolume: 2e-8 * unit.meter ** 3 / unit.mole,
}

PROPERTIES = [Density, EnthalpyOfVaporization, ExcessMolarVolume, DielectricConstant, EnthalpyOfMixing,]


parser = argparse.ArgumentParser("Run evaluator on single entry")
parser.add_argument("entry", type=int)




if __name__ == "__main__":
    args = parser.parse_args()

    # read DF from file
    filename = ENTRY_NAME.format(args.entry)
    df = pd.read_csv(filename, index_col=0)
    df = df.dropna(axis=1)
    print(df)
    print(df.columns)
    dataset = ThermoMLDataSet.from_pandas(df)
    proptypes = {type(x) for x in dataset}

    # set up request options
    estimation_options = RequestOptions()
    for prop in PROPERTIES:
        if prop in proptypes:
            schema = prop.default_simulation_schema(absolute_tolerance=ABSOLUTE_TOLERANCES[prop])
            estimation_options.add_schema("SimulationLayer", prop, schema)

    # set up server shit
    resources = ComputeResources(number_of_gpus=1,
                                 preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA)
    backend = DaskLocalCluster(number_of_workers=1, resources_per_worker=resources)
    backend.start()
    storage = LocalFileStorage(root_directory=STORED)
    server = EvaluatorServer(calculation_backend=backend,
                             storage_backend=storage,
                             working_directory=WORKING)
    server.start(asynchronous=True)
    print(server)

    # set up client and FINALLY submit job?
    client = EvaluatorClient()
    request, exception = client.request_estimate(
        property_set=dataset,
        force_field_source=force_field_source,
        options=estimation_options,
    )
    assert exception is None, str(exception)

    # Wait for the results.
    results, exception = request.results(synchronous=True, polling_interval=30)
    assert exception is None, str(exception)

    print(results)

    json_name = JSON_NAME.format(args.entry)
    results.estimated_properties.json(json_name, format=True)
    print(f"Wrote to {json_name}")