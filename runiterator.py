from job import job
import os

savedirname = 'analysis_data'

# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if 'job' not in path:
        continue
    if savedirname in path:
        continue

    run = job(path)

    run.apd()
    run.etg()
    run.vtg()
    run.save_data()

    print('-'*79)
