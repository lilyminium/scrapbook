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

from psiresp.grid import GridOptions
from psiresp.orientation import OrientationEsp
from psiresp.constraint import ConstraintMatrix
from psiresp.resp import RespCharges, RespOptions
from psiresp.charge import MoleculeChargeConstraints, ChargeConstraintOptions, ChargeSumConstraint
from psiresp.molecule import Molecule, Atom
from psiresp import psi4utils

from psiresp.tests.datafiles import DMSO
import qcelemental as qcel
qcmol = qcel.models.Molecule.from_file(DMSO, dtype="xyz")
mol = Molecule(qcmol=qcmol)
atom = Atom(molecule=mol, index=0)
a = ChargeSumConstraint(charge=0.3, atoms=[atom])
print(a.dict())