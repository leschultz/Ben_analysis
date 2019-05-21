#!/usr/bin/env python3

'''
Construct ICO at Tg dataframe for analysis data.
'''

import pandas as pd
import numpy as np

import sys
import os

jobs_dir = sys.argv[1]  # Directory where run analysis data is stored
ico_tg_file = sys.argv[2]  # ico file name
export_dir = sys.argv[3]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

rows = []
for item in os.walk(jobs_dir):
    path = item[0]
    files = item[-1]

    # Check if files are present
    if ico_tg_file not in files:
        continue

    split = path.split('/')
    system = split[-5]
    composition = split[-4]
    tg = float(split[-3])
    job = split[-2]

    ico_tg_path = os.path.join(path, ico_tg_file)

    ico_tg = np.loadtxt(ico_tg_path)

    row = [system, composition, tg, job, ico_tg]

    rows.append(row)

df = pd.DataFrame(rows)
df.columns = ['system', 'composition', 'tg', 'job', 'ico_tg']

df.to_csv(os.path.join(export_dir, 'ico_tg_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'ico_tg_df.html'), index=False)
