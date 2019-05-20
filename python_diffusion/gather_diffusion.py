#!/usr/bin/env python3

'''
Construct diffusion dataframe for analysis data.
'''

from functools import reduce

import pandas as pd

import sys
import os

jobs_dir = sys.argv[1]  # Directory where run analysis data is stored
diffusion_file = sys.argv[2]  # Tg file name
export_dir = sys.argv[3]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

rows = []
for item in os.walk(jobs_dir):
    path = item[0]
    files = item[-1]

    # Check if files are present
    if diffusion_file not in files:
        continue

    split = path.split('/')
    system = split[-5]
    composition = split[-4]
    steps = float(split[-3])
    job = split[-2]

    diffusion_path = os.path.join(path, diffusion_file)

    diffusion = pd.read_csv(diffusion_path)
    diffusion['system'] = system
    diffusion['composition'] = composition
    diffusion['job'] = job

    rows.append(diffusion)

df = pd.concat(rows)
df = df.reset_index(drop=True)
print(df)

df.to_csv(os.path.join(export_dir, 'diffusion_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'diffusion_df.html'), index=False)
