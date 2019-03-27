from trajsteps import trajectorysteps
from infoparser import inputinfo
from outimport import outdata

import pymatgen as mg
import pandas as pd
import numpy as np

import os


def volume(r):
    '''
    Calculate the atomic volume from a radius. This assumes that the
    volume is spherical.
    '''

    return 4./3.*np.pi*r**3.


# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue

    # Determine the composition based on binary name
    dirs = path.split('/')
    print(dirs)
    elements = dirs[1].split('-')
    els = [mg.Element(i) for i in elements]

    atomicradii = [i.atomic_radius*10. for i in els]  # in Am
    atomicvol = [volume(i) for i in atomicradii]  # in Am^3
    print(els)
    print(atomicradii)
    print(atomicvol)
    
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

    print(path)
    print(df)
