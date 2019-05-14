from scipy.signal import argrelextrema
from matplotlib import pyplot as pl

from job import job

import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # The job directories
datadir_name = sys.argv[2]  # The name of directory containing data
trajdotlammpstrj = sys.argv[3]  # Trajectories
depdotin = sys.argv[4]  # Input file
export_dir = sys.argv[5]  # The export directory
datadirname = sys.argv[6]  # Name of data directory

# Loop for each path
for item in os.walk(jobs_dir):

    path = item[0]

    split = path.split('/')

    # Filter for paths that contain jobs
    if datadir_name not in split[-1]:
        continue

    run = job(path, export_dir, datadirname)

    traj = os.path.join(path, trajdotlammpstrj)
    dep = os.path.join(path, depdotin)

    run.apd_last(traj, dep)

    print('-'*79)
