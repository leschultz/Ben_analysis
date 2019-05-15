#!/usr/bin/env python3

'''
Construct Tg dataframe for analysis data.
'''

import pandas as pd
import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # Directory where run analysis data is stored
tg_file = sys.argv[2]  # Tg file name
temp_cutoff_file = sys.argv[3]  # Temperature cutoff file name
export_dir = sys.argv[4]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

rows = []
for item in os.walk(jobs_dir):
    path = item[0]
    files = item[-1]

    # Check if files are present
    if tg_file not in files:
        continue
    if temp_cutoff_file not in files:
        continue

    split = path.split('/')
    system = split[-5]
    composition = split[-4]
    steps = int(split[-3])
    job = split[-2]

    tg_path = os.path.join(path, tg_file)
    temp_cutoff_path = os.path.join(path, temp_cutoff_file)

    tg = np.loadtxt(tg_path)
    temp_cutoff = np.loadtxt(temp_cutoff_path)

    row = [system, composition, steps, job, tg, temp_cutoff]

    rows.append(row)

df = pd.DataFrame(rows)
df.columns = ['system', 'composition', 'steps', 'job', 'tg', 'temp_cutoff']

df.to_csv(os.path.join(export_dir, 'tg_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'tg_df.html'), index=False)
