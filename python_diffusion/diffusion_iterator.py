#!/usr/bin/env python3

from scipy.signal import argrelextrema
from matplotlib import pyplot as pl

from job import job

import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # The job directories
job_name = sys.argv[2]  # The generic name for jobs
export_dir = sys.argv[3]  # The export directory
datadirname = sys.argv[4]  # Name of data directory
plotdirname = sys.argv[5]  # Name of plot directory

trajdotlammpstrj = sys.argv[6]  # Trajectories
testdotout = sys.argv[7]  # LAMMPS print to screen
depdotin = sys.argv[8]  # Input file
alpha = float(sys.argv[9])  # The significance level for t-test

# Loop for each path
for item in os.walk(jobs_dir):

    path = item[0]

    split = path.split('/')

    # Filter for paths that contain jobs
    if job_name not in split[-1]:
        continue

    run = job(path, export_dir, datadirname, plotdirname)

    run.input_file(depdotin)
    run.sys(testdotout)
    run.box(trajdotlammpstrj)

    run.diffusion(alpha=alpha)

    print('-'*79)
