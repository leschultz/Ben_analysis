from matplotlib import pyplot as pl

import pandas as pd

import os

savedirname = 'analysis_data'
apddottxt = 'apd.txt'

# Loop for each path
for item in os.walk('../'):

    path = item[0]

    # Filter for paths that contain jobs
    if savedirname not in path:
        continue

    df = pd.read_csv(os.path.join(path, apddottxt))
    print(df)
    print(path)
