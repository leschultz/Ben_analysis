from PyQt5 import QtGui  # Added to be able to import ovito

from trajsteps import trajectorysteps
from infoparser import inputinfo
from outimport import outdata

import pymatgen as mg
import pandas as pd

import os

# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue

    # Location of files
    system = os.path.join(path, 'test.out')  # Output file
    traj = os.path.join(path, 'traj.lammpstrj')  # Trajectories
    dep = os.path.join(path, 'dep.in')  # Input file

    # Important paramters from the input file
    param = inputinfo(dep)

    # Thermodynamic data from test.out file
    dfsys = outdata(system)

    # The steps where trajectories where exported
    trajsteps = trajectorysteps(traj)

    # Match dfsys steps with trajsteps
    df = dfsys[dfsys['Step'].isin(trajsteps)]
