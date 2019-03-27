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


testdotout = 'test.out'
trajdotlammpstrj = 'traj.lammpstrj'
depdotin = 'dep.in'

savedirname = 'analysis_data'

# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue
    if savedirname in path:
        continue

    # Location of files
    file_system = os.path.join(path, testdotout)  # Output file
    file_trajs = os.path.join(path, trajdotlammpstrj)  # Trajectories
    file_dep = os.path.join(path, depdotin)  # Input file

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

    # Element counts from actual element
    elements = {}
    for key, count in counts.items():
        element = mg.Element(depparam['elements'][key])
        atomicradii = element.atomic_radius  # in Am
        atomicvol = volume(atomicradii)  # in Am^3

        elements[depparam['elements'][key]] = {
                                               'counts': count,
                                               'radius': atomicradii,
                                               'volume': atomicvol,
                                               }

    dfel = pd.DataFrame(elements).T

    # Calculate the atomic packing density (APD)
    numerator = np.sum(dfel['counts']*dfel['volume'])

    df['APD'] = numerator/df['Volume']

    # Cooling data starts after hold1
    dfcool = df[df['Step'] >= depparam['hold1']]

    # Create a directory for the analysis files
    savepath = os.path.join(path, savedirname)
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    dfcool.to_csv(os.path.join(savepath, 'apd.txt'), index=False)  # Export APD
