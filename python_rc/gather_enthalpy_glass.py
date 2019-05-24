#!/usr/bin/env python3

import pandas as pd

import sys
import os

jobs_dir = sys.argv[1]  # The job directories
export_dir = sys.argv[2]  # The export directory
crystal_ref_run = sys.argv[3]  # The name of the input file
logdotlammps = sys.argv[4]  # The name of the log file

# Make export directory
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

# Loop for each path
df = []
for item in os.walk(jobs_dir):

    path = item[0]
    files = item[2]

    if not ((crystal_ref_run in files) and (logdotlammps in files)):
        continue

    if crystal_ref_run in files:
        with open(os.path.join(path, crystal_ref_run)) as f:
            for line in f:
                line = line.strip().split(' ')
                if 'thermo_style' == line[0]:
                    line = [i for i in line if i != '']
                    columns = line[2:]

    data = []
    if logdotlammps in files:

        with open(os.path.join(path, logdotlammps)) as f:
            for line in f:
                line = line.strip().split(' ')
                line = [i for i in line if '' != i]

                if len(line) != len(columns):
                    continue

                try:
                    row = [float(i) for i in line]
                    data.append(row)

                except Exception:
                    pass

    split = path.split('/')

    phase = split[-6]
    system = split[-5]
    composition = split[-4]
    steps = split[-3]
    job = split[-2]

    dfdata = pd.DataFrame(data, columns=columns)
    enthalpy = dfdata['enthalpy'].values[-1]

    row = [system, composition, steps, job, enthalpy]
    df.append(row)

df = pd.DataFrame(df)
df.columns = ['system', 'composition', 'steps', 'job', 'enthalpy']

df.to_csv(os.path.join(export_dir, 'enthalpy_glass_df.txt'), index=False)
df.to_html(os.path.join(export_dir, 'enthalpy_glass_df.html'), index=False)
