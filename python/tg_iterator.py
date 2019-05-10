from job import job

from matplotlib import pyplot as pl
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

    #try:
    run = job(path, datadirname, plotdirname)

    run.input_file(depdotin)
    dfcool = run.sys(testdotout)
    run.box(trajdotlammpstrj)

    range_t = np.linspace(min_t, max_t, 100)

    tgs = []
    for t in range_t:
        tg = run.etg(max_temp=t, write=False, plot=False)
        tgs.append(tg)

    tgs = np.array(tgs)
    index = np.argmax(tgs)

    tcut = range_t[index]

    run.etg(max_temp=tcut)

    fig, ax = pl.subplots()

    ax.plot(range_t, tgs, marker='.', linestyle='none', label='Data')
    ax.axvline(tcut, color='k', label='Upper Cutoff Temperature [K]')

    ax.set_xlabel('Upper Temperature Cutoff [K]')
    ax.set_ylabel('Tg [K]')

    ax.grid()
    ax.legend()

    fig.tight_layout()

    plot_path = os.path.join(path, plotdirname)
    fig.savefig(os.path.join(plot_path, 'tg_upper_t_cutoff'))

    pl.close('all')
    #except Exception:
    #    pass

    print('-'*79)
