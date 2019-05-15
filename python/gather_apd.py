#!/usr/bin/env python3

'''
Construct APD dataframe for analysis data.
'''

import pandas as pd
import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # Directory where run analysis data is stored
apd_file = sys.argv[2]  # Tg file name
export_dir = sys.argv[3]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

rows = []
for item in os.walk(jobs_dir):
    path = item[0]
    files = item[-1]

    # Check if files are present
    if apd_file not in files:
        continue

    split = path.split('/')
    split.pop(-1)
    system = split[-5]
    composition = split[-4]
    steps = int(split[-3])
    job = split[-2]

    apd_path = os.path.join(path, apd_file)

    apd = np.loadtxt(apd_path)

    row = [system, composition, steps, job, apd]

    rows.append(row)

df = pd.DataFrame(rows)
df.columns = ['system', 'composition', 'steps', 'job', 'apd']

df.to_csv(os.path.join(export_dir, 'apd_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'apd_df.html'), index=False)
