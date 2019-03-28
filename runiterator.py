from job import job
import os

datadirname = 'analysis_data'
plotdirname = 'analysis_plots'


# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue
    if datadirname in path:
        continue
    if plotdirname in path:
        continue

    run = job(path)

    run.apd()
    run.etg()
    run.vtg()
    run.save_data()

    print('-'*79)
