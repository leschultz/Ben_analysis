from job import job

from matplotlib import pyplot as pl

from scipy.signal import argrelextrema
import numpy as np

import sys
import os

datadirname = 'analysis_data'
plotdirname = 'analysis_plots'
minfile = os.path.join('100K_Structure_minimization', 'finaltraj.lammpstrj')
mininfile = os.path.join(
                         '100K_Structure_minimization',
                         '100k_minimize_template.in'
                         )

# The name of important files for each job
trajdotlammpstrj = 'traj.lammpstrj'
testdotout = 'test.out'
depdotin = 'dep.in'

min_t = float(sys.argv[2])  # Min temp to start optimization
max_t = float(sys.argv[3])  # Max temp to start optimization

# Loop for each path
for item in os.walk(sys.argv[1]):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue
    if datadirname in path:
        continue
    if plotdirname in path:
        continue
    if 'minimization' in path:
        continue

    run = job(path, datadirname, plotdirname)

    run.input_file(depdotin)
    dfcool = run.sys(testdotout)
    run.box(trajdotlammpstrj)

    temps = dfcool['Temp'].values
    min_t = temps[-6]  # Minimum number of point for spline to work
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
