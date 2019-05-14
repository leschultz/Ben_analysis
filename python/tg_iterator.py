from scipy.signal import argrelextrema
from matplotlib import pyplot as pl

from job import job

import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # The job directories
jobs_name = sys.argv[2]  # The generic job name

datadirname = sys.argv[3]  # Name of data directory
plotdirname = sys.argv[4]  # Name of plot directory

# The name of important files for each job
trajdotlammpstrj = sys.argv[5]  # Trajectories
testdotout = sys.argv[6]  # LAMMPS print to screen
depdotin = sys.argv[7]  # Input file

# Loop for each path
for item in os.walk(jobs_dir):

    path = item[0]

    split = path.split('/')

    # Filter for paths that contain jobs
    if 'job' not in split[-1]:
        continue

    run = job(path, datadirname, plotdirname)

    run.input_file(depdotin)
    dfcool = run.sys(testdotout)
    run.box(trajdotlammpstrj)

    temps = dfcool['Temp'].values
    min_t = temps[-6]  # Minimum number of points for spline to work
    max_t = max(temps)

    range_t = np.linspace(min_t, max_t, 100)

    tgs = []
    for t in range_t:
        tg = run.etg(max_temp=t, write=False, plot=False, verbose=False)
        tgs.append(tg)

    tgs = np.array(tgs)

    fig, ax = pl.subplots()

    ax.plot(range_t, tgs, marker='.', linestyle='none', label='Data')

    ax.set_xlabel('Upper Temperature Cutoff [K]')
    ax.set_ylabel('Tg [K]')

    ax.grid()

    fig.tight_layout()

    # Go full screen
    mng = pl.get_current_fig_manager()
    mng.full_screen_toggle()

    plot_path = os.path.join(path, plotdirname)

    click = fig.ginput(n=-1, mouse_add=1, mouse_pop=3)
    tcut = click[-1][0]

    ax.axvline(
               tcut,
               linestyle=':',
               color='k',
               label='Upper Temperature Cutoff: '+str(tcut)+' [K]'
               )

    ax.legend(loc='upper center')

    fig.savefig(os.path.join(plot_path, 'tg_upper_t_cutoff'))

    pl.close('all')

    print('Upper cutoff temperature set to '+str(tcut)+' [K]')

    run.etg(max_temp=tcut)

    print('-'*79)
