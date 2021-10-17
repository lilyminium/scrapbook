#!/usr/bin/env python

import pathlib
import json
import glob
import os
import argparse

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
parser.add_argument("mol_index", type=int)
parser.add_argument("--entry_index", type=int, default=None)
parser.add_argument("--component", type=str, default="solvent")

if "OE_LICENSE" in os.environ:
    tk = "oe"
else:
    tk = "at"

class TestMolecule:

    def __init__(self, molname="0010", component="solvent", directory="../06_resp/07_results",
                 entry_index=None,
                 csv="01_mnsol_data.csv"):
        self.molname = molname
        self.mol_index = int(molname)
        self.entry_index = self.mol_index if entry_index is None else int(entry_index)
        self.component = component

        self.prefix =  f"mol-{molname}_{component}_entry-{self.entry_index}"
        self.indir = pathlib.Path(directory) / f"mol-{molname}_{component}"
        self.outdir = pathlib.Path(f"03_results/{self.prefix}")
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.ffdir = self.outdir / "forcefields"
        self.ffdir.mkdir(parents=True, exist_ok=True)
        self.chargedir = self.outdir / "charges"
        self.chargedir.mkdir(parents=True, exist_ok=True)
        self.property = property_from_csv(csv, index=self.entry_index)
        self.smiles = self.get_smiles()
        self.load_charges()
    
    def load_charges(self, pattern=r"04_w-{weight}_stage-{stage}_{state}_charges.dat"):
        qcmol = self.to_qcmol()
        n_atom = len(qcmol.geometry)
        
        df_dct = {"State": [], "Weight": [], "Stage": [], "Atom": [], "Charge": []}

        for state in ["gas", "water"]:
            for weight in [1, None]:
                for stage in 1, 2:
                    filename = self.indir / pattern.format(stage=stage, weight=weight, state=state)
                    charges = np.loadtxt(str(filename))[:n_atom]
                    setattr(self, f"stage_{stage}_weight_{weight}_{state}_charges", charges)

                    df_dct["Charge"].extend(charges)
                    df_dct["Atom"].extend(np.arange(n_atom))
                    df_dct["Stage"].extend([stage] * n_atom)
                    df_dct["Weight"].extend([weight] * n_atom)
                    df_dct["State"].extend([state] * n_atom)
        
        df = pd.DataFrame(df_dct)
        df.to_csv(str(self.outdir / f"{self.molname}_{self.component}_charges.csv"))
        self.charge_df = df

    
    def load_energies(self, pattern="03_c*energy.json"):
        from qcfractal.interface.models.records import ResultRecord
        
        files = sorted(glob.glob(str(self.indir / pattern)))

        self.records = [ResultRecord.parse_file(file) for file in files]
        self.energies = np.array([rec.properties.return_energy for rec in self.records])

    def get_charges(self, delta=0.0, weight=1):
        # delta = % aqueous
        gas = getattr(self, f"stage_2_weight_{weight}_gas_charges")
        water = getattr(self, f"stage_2_weight_{weight}_water_charges")

        charges = (delta * water) + ((1 - delta) * gas)
        return charges

    def get_smiles(self):
        df = pd.read_csv("01_mnsol_data.csv")
        row = df.iloc[self.mol_index]
        if row["Role 1"].lower() == self.component:
            smiles = row["Component 1"]
        else:
            smiles = row["Component 2"]
        return smiles


    def to_qcmol(self, pattern=r"03_c{conf:02d}_o{orient:02d}_mol.json", conf=1, orient=1):
        file = self.indir / pattern.format(conf=conf, orient=orient)
        qcmol = qcel.models.Molecule.from_file(file)
        return qcmol

    def to_offmol(self, conf=1, orient=1):
        qcmol = self.to_qcmol(conf=conf, orient=orient)

        smiles_parser = Chem.rdmolfiles.SmilesParserParams()
        smiles_parser.removeHs = False
        rdmol = Chem.AddHs(Chem.MolFromSmiles(self.smiles, smiles_parser))
        rdel = [at.GetSymbol() for at in rdmol.GetAtoms()]
        assert_equal(rdel, qcmol.symbols)

        offmol = Molecule.from_rdkit(rdmol)
        positions = qcmol.geometry * qcel.constants.conversion_factor("bohr", "angstrom")
        offmol._add_conformer(unit.Quantity(positions, unit.angstrom))

        return offmol

    def generate_resp2_forcefield(self, conf=1, orient=1, delta=0.0, weight=1):
        offmol = self.to_offmol(conf=conf, orient=orient)
        charges = self.get_charges(delta=delta, weight=weight)
        chargefile = f"{self._resp2_name(delta=delta, weight=weight)}_charges.txt"
        np.savetxt(str(self.chargedir / chargefile), charges)
        charges = [x * unit.elementary_charge for x in charges]
        smiles = offmol.to_smiles(mapped=True)

        ff = ForceField("1.3.0.offxml")
        handler = ff.get_parameter_handler("LibraryCharges")
        handler.add_parameter({"charge": charges, "smirks": smiles})
        return ff

    def generate_am1bcc_forcefield(self, conf=1, orient=1, **kwargs):
        offmol = self.to_offmol(conf=conf, orient=orient)
        offmol.assign_partial_charges("am1bcc", use_conformers=offmol.conformers)
        charges = offmol.partial_charges

        _charges = [x / unit.elementary_charge for x in charges]
        chargefile = f"{self._am1bcc_name(conf=conf, orient=orient)}_charges.txt"
        np.savetxt(str(self.chargedir / chargefile), _charges)
        smiles = offmol.to_smiles(mapped=True)

        ff = ForceField("1.3.0.offxml")
        handler = ff.get_parameter_handler("LibraryCharges")
        handler.add_parameter({"charge": charges, "smirks": smiles})
        return ff



    def _generate_workflow(self, ff, filename):
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

    def generate_resp2_workflow(self, delta=0.0, weight=1):
        filename = f"{self._resp2_name(delta=delta, weight=weight)}_ff.json"
        filename = str(self.ffdir / filename)
        ff = self.generate_resp2_forcefield(delta=delta, weight=weight)
        return self._generate_workflow(ff, filename)

    def generate_am1bcc_workflow(self, conf=1, orient=1):
        filename = f"{self._am1bcc_name(conf=conf, orient=orient)}_ff.json"
        filename = str(self.ffdir / filename)
        ff = self.generate_am1bcc_forcefield(conf=conf, orient=orient)
        return self._generate_workflow(ff, filename)

    def _resp2_name(self, delta=0.0, weight=1):
        return f"{self.prefix}_weight-{weight}_delta-{delta:.02f}_resp2"

    def _am1bcc_name(self, conf=1, orient=1):
        return f"{self.prefix}_c{conf}_o{orient}_am1bcc_{tk}"


if __name__ == "__main__":
    args = parser.parse_args()

    # set up molecule
    molname = f"{args.mol_index:04d}"
    emol = TestMolecule(molname=molname, component=args.component, entry_index=args.entry_index)


    working_directory = emol.outdir / "evaluator_working-data"
    working_directory.mkdir(exist_ok=True)
    storage_directory = emol.outdir / "evaluator_stored-data"
    storage_directory.mkdir(exist_ok=True)

    resources = ComputeResources(number_of_gpus=1, number_of_threads=12,
                                 preferred_gpu_toolkit=ComputeResources.GPUToolkit.CUDA)


    with DaskLocalCluster(number_of_workers=1, resources_per_worker=resources) as backend:
        for weight in [1, None]:
            for delta in np.linspace(0, 1, 21):
                workflow = emol.generate_resp2_workflow(delta=delta, weight=weight)
                results = workflow.execute(root_directory=str(working_directory),
                                           calculation_backend=backend,
                                           ).result()
                jsonfile = f"{emol._resp2_name(delta=delta, weight=weight)}_results.json"
                file = str(storage_directory / jsonfile)
                results.json(file, format=True)
                print(f"Wrote results to {file}")

        for conf in range(1, 11):
            for orient in range(1, 13):
                workflow = emol.generate_am1bcc_workflow(conf=conf, orient=orient)
                results = workflow.execute(root_directory=str(working_directory),
                                           calculation_backend=backend,
                                           ).result()
                jsonfile = f"{emol._am1bcc_name(conf=conf, orient=orient)}_results.json"
                file = str(storage_directory / jsonfile)
                results.json(file, format=True)
                print(f"Wrote results to {file}")


