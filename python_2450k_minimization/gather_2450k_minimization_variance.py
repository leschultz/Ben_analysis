#!/usr/bin/env python3

'''
Construct diffusion dataframe for analysis data.
'''

from functools import reduce

import pandas as pd

import sys
import os

jobs_dir = sys.argv[1]  # Directory where run analysis data is stored
variance_file = sys.argv[2]  # Tg file name
export_dir = sys.argv[3]  # The export directory

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

rows = []
for item in os.walk(jobs_dir):
    path = item[0]
    files = item[-1]

    # Check if files are present
    if variance_file not in files:
        continue

    split = path.split('/')
    system = split[-7]
    composition = split[-6]
    steps = float(split[-5])
    job = split[-4]

    variance_path = os.path.join(path, variance_file)

    variance = pd.read_csv(variance_path)
    variance['system'] = system
    variance['composition'] = composition
    variance['steps'] = steps
    variance['job'] = job

    rows.append(variance)

df = pd.concat(rows)
df = df.reset_index(drop=True)

df.to_csv(
          os.path.join(export_dir, 'variance_2450k_minimization_df.txt'),
          index=False
          )
df.to_html(
           os.path.join(export_dir, 'variance_2450k_minimization_df.html'),
           index=False
           )
