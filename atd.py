import pymatgen as mg
import pandas as pd
import numpy as np

import traj
import test
import dep

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
    elements = dirs[1].split('-')
    els = [mg.Element(i) for i in elements]

    atomicradii = [i.atomic_radius*10. for i in els]  # in Am
    atomicvol = [volume(i) for i in atomicradii]  # in Am^3
    
    # Location of files
    file_system = os.path.join(path, 'test.out')  # Output file
    file_trajs = os.path.join(path, 'traj.lammpstrj')  # Trajectories
    file_dep = os.path.join(path, 'dep.in')  # Input file

    # Important paramters from the input file
    depparam = dep.info(file_dep)

    # Information from traj.lammpstrj file
    dftraj, counts = traj.info(file_trajs)

    # Thermodynamic data from test.out file
    dfsys = test.info(file_system)

    # Merge dfsys and dftraj on matching steps
    df = pd.merge(dfsys, dftraj, on=['Step'])

    # Find the box lengths
    df['dx'] = df['xhi']-df['xlo']  # Am
    df['dy'] = df['yhi']-df['ylo']  # Am
    df['dz'] = df['zhi']-df['zlo']  # Am

    # Find the box volume
    df['Volume'] = df['dx']*df['dy']*df['dz']  # Am^3

    print(counts)
