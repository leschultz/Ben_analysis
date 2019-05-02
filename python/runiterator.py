from job import job

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

    try:
        run = job(path)

        run.input_file(depdotin)
        run.sys(testdotout)
        run.box(trajdotlammpstrj)

        #run.apd()
        run.etg()
        #run.vtg()
        #run.apd_single(
        #               os.path.join(path, minfile),
        #               os.path.join(path, mininfile)
        #               )
        #run.vp()

        run.save_data()

    except Exception:
        pass

    print('-'*79)
