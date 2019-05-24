#!/usr/bin/env python3

import pandas as pd

import sys
import os

jobs_dir = sys.argv[1]  # The job directories
export_dir = sys.argv[2]  # The export directory
crystal_ref_run = sys.argv[3]  # The name of the input file
logdotlammps = sys.argv[4]  # The name of the log file

# Loop for each path
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

    df = pd.DataFrame(data, columns=columns)
    enthalpy = df['enthalpy'].values[-1]

    # Save paths
    export_path = os.path.join(export_dir, path.strip('../'))

    if not os.path.exists(export_path):
        os.makedirs(export_path)

    write_name = os.path.join(export_path, 'enthalpy.txt')
    with open(write_name, 'w+') as outfile:
        outfile.write(str(enthalpy))
