#!/usr/bin/env python

from utils import property_from_csv

from openff.evaluator.datasets import PhysicalProperty
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField


def generate_charges(offmol: Molecule):


def generate_charged_forcefield(physprop: PhysicalProperty,
                                charge_function: callable):
    for component in physprop.substances:
        mol =


def workflow_from_entry(entry: PhysicalProperty):
    from openff.evaluator.forcefield import SmirnoffForceFieldSource
    property_type = type(entry)

    # Set up schema. Here I know the force field is SMIRNOFF
    schema = property_type.default_simulation_schema().workflow_schema
    schema.replace_protocol_types({"BaseBuildSystem": "BuildSmirnoffSystem"})

    # force field
    FF_FILE = "force-field.json"
    ff_source = SmirnoffForceFieldSource.from_object(ff)
    ff_source.json(FF_FILE)

    metadata = Workflow.generate_default_metadata(entry, FF_FILE, [])
    workflow = Workflow.from_schema(schema, metadata=metadata)
    return workflow
