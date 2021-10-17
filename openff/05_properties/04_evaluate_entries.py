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

from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.evaluator.backends import ComputeResources
from openff.evaluator.backends.dask import DaskLocalCluster
from openff.toolkit.topology.molecule import Molecule, unit
from openff.evaluator.workflow import Workflow


from utils import property_from_csv


parser = argparse.ArgumentParser("Evaluate properties")
parser.add_argument("entry_index", type=int)
parser.add_argument("--component", type=str, default="solvent")

METHODS = {
    "hf": {"method": "hf", "basis": "6-31g*"},
    "default": {"method": "b3lyp-d3bj", "basis": "dzvp"},
    "resp2": {"method": "PW6B95", "basis": "cc-pV(D+d)Z"},
}

METHOD_FULL = {
    k: METHODS[k]["method"] + "/" + METHODS[k]["basis"]
    for k in METHODS
}

REVERSE_METHOD_FULL = {
    v: k for k, v in METHOD_FULL.items()
}

class TestEntry:
    def __init__(self, entry_index, csv="01_mnsol_data.csv"):
        self.entry_index = int(entry_index)
        self.prefix = f"entry_{entry_index:04d}"
        self.outdir = pathlib.Path(f"04_results/{self.prefix}")
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.wkdir = self.outdir / "evaluator_working-data"
        self.wkdir.mkdir(parents=True, exist_ok=True)
        self.results_dir = self.outdir / "evaluated_results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.ffdir = self.outdir / "forcefields"
        self.ffdir.mkdir(parents=True, exist_ok=True)
        self.property = property_from_csv(csv, index=self.entry_index)

        comp1, comp2 = self.property.substance

        self.component_1_df = pd.read_csv(f"../06_resp/08_results/{comp1.smiles}/05_calculated_charges.csv")
        self.component_2_df = pd.read_csv(f"../06_resp/08_results/{comp2.smiles}/05_calculated_charges.csv")

    
    def run(self, backend):
        ffs = self.generate_resp2_ffs()
        ffs += self.generate_am1bcc_ffs()

        for ff in tqdm.tqdm(ffs):
            self.execute_single_workflow(backend, **f)



    def execute_single_workflow(self, backend, name=None, ff=None):
        ff_filename = self.ffdir / f"{name}_ff.json"
        workflow = self._generate_workflow(ff=ff, filename=ff_filename)
        results = workflow.execute(root_directory=str(self.wkdir),
                                    calculation_backend=backend,
                                    ).result()
        jsonfile = str(self.results_dir / f"{name}_results.json")
        results.json(jsonfile, format=True)
        

    def _generate_workflow(self, ff=None, filename=None):
        from openff.evaluator.forcefield import SmirnoffForceFieldSource
        property_type = type(self.property)
        filename = os.path.abspath(filename)

        # Set up schema. Here I know the force field is SMIRNOFF
        schema = property_type.default_simulation_schema().workflow_schema
        schema.replace_protocol_types({"BaseBuildSystem": "BuildSmirnoffSystem"})

        # force field
        ff_source = SmirnoffForceFieldSource.from_object(ff)
        ff_source.json(filename)

        metadata = Workflow.generate_default_metadata(self.property, filename, [])
        workflow = Workflow.from_schema(schema, metadata=metadata)
        return workflow



    def generate_resp2_ffs(self):
        COLS = ["Method", "Weight", "Delta", "Grid"]

        df1 = self.component_1_df[self.component_1_df == "RESP2"]
        df2 = self.component_2_df[self.component_2_df == "RESP2"]

        forcefields = []

        for (method, weight, delta, grid), df1_ in df1.groupby(by=COLS):
            df1_ = df1_.sort_values("charge", inplace=False)
            df2_ = df2[(df2.Method == method) & (df2.Weight == weight) & (df2.Delta == delta) & (df2.Grid == grid)]
            df2_ = df2_.sort_values("charge", inplace=False)

            assert len(df1_.charge.values) == max(df1_["Atom number"].values)
            assert len(df2_.charge.values) == max(df2_["Atom number"].values)

            chg1 = [x * unit.elementary_charge for x in df1_.charge.values]
            smirks1 = df1_.smirks.values[0]

            chg2 = [x * unit.elementary_charge for x in df2_.charge.values]
            smirks2 = df2_.smirks.values[0]

            ff = ForceField("1.3.0.offxml")
            handler = ff.get_parameter_handler("LibraryCharges")
            handler.add_parameter({"charge": chg1, "smirks": smirks1})
            handler.add_parameter({"charge": chg2, "smirks": smirks2})


            forcefields.append({"name": f"RESP2_m-{REVERSE_METHOD_FULL[method]}_g-{grid}_w-{weight}_d-{delta}" 
                                "ff": ff})

        return forcefields
        

    def generate_am1bcc_ffs(self):
        COLS = ["Method", "Conformer", "Orientation"]

        df1 = self.component_1_df[self.component_1_df == "AM1BCC"]
        df2 = self.component_2_df[self.component_2_df == "AM1BCC"]

        forcefields = []

        for (tk, conf, orient), df1_ in df1.groupby(by=COLS):
            df1_ = df1_.sort_values("charge", inplace=False)
            df2_ = df2[(df2.Method == tk) & (df2.Conformer == conf) & (df2.Orientation == orient)]
            df2_ = df2_.sort_values("charge", inplace=False)

            assert len(df1_.charge.values) == max(df1_["Atom number"].values)
            assert len(df2_.charge.values) == max(df2_["Atom number"].values)

            chg1 = [x * unit.elementary_charge for x in df1_.charge.values]
            smirks1 = df1_.smirks.values[0]

            chg2 = [x * unit.elementary_charge for x in df2_.charge.values]
            smirks2 = df2_.smirks.values[0]

            ff = ForceField("1.3.0.offxml")
            handler = ff.get_parameter_handler("LibraryCharges")
            handler.add_parameter({"charge": chg1, "smirks": smirks1})
            handler.add_parameter({"charge": chg2, "smirks": smirks2})


            forcefields.append({"name": f"AM1BCC_{tk}_c-{conf}_o-{orient}" 
                                "ff": ff})

        return forcefields


if __name__ == "__main__":
    args = parser.parse_args()
    entry = TestEntry(args.entry_index)
    resources = ComputeResources(number_of_gpus=1, number_of_threads=12,
                                 preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA)

    with DaskLocalCluster(number_of_workers=1, resources_per_worker=resources) as backend:
        entry.run(backend)
