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

    print(path)

    run = job(path)
    run.apd()
    run.save_data(['system', 'box', 'apd'])
