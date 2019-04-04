from scipy import stats as st

from os.path import join

import pandas as pd
import numpy as np

import sys
import os

datadirname = 'analysis_data'
externaldirname = 'data_external'

etgfile = 'tg_e.txt'
vtgfile = 'tg_v.txt'

cxpath = [sys.argv[1], externaldirname, 'All_jobs_with_crystal_info.csv']
dfcx = pd.read_csv(join(*cxpath))

# Modify Ben's columns for merging datasets
alloys = dfcx['alloy']
elements = [i.split('-')[1] for i in alloys]
comp = dfcx['Comp']
comp = ['{:.2f}'.format(np.round(i, 2)) for i in comp]
comp = [i+j for i, j in zip(elements, comp)]
dfcx['Comp'] = comp

dfcx = dfcx[['alloy', 'Comp', 'Thold', 'Job', 'Run Crystallized']]
dfcx.columns = [
                'System',
                'Composition [decimal]',
                'Steps [-]',
                'Job',
                'Crystallization'
                ]

# Build a dataframe for all collected Tg values
columns = [
           'System',
           'Composition [decimal]',
           'Steps [-]',
           'Job',
           'Tg from E-3kT Curve [K]',
           'Tg from Specific Volume Curve [K]',
           ]

dftg = pd.DataFrame(columns=columns)

count = 0
for item in os.walk(sys.argv[1]):

    if datadirname not in item[0]:
        continue

    if '-' not in item[0]:
        continue

    names = item[0].split('/')[-5:-1]

    epath = join(item[0], etgfile)
    if os.path.exists(epath):
        etg = np.loadtxt(epath, dtype=float)
    else:
        etg = np.nan

    vpath = join(item[0], vtgfile)
    if os.path.exists(vpath):
        vtg = np.loadtxt(vpath, dtype=float)
    else:
        vtg = np.nan

    row = names+[etg, vtg]

    dftg.loc[count] = row  # Append a row to the dataframe

    count += 1

dftg['Steps [-]'] = dftg['Steps [-]'].apply(pd.to_numeric)

# The columns to merge crystallization and tg dataframes
mergecolumns = [
                'System',
                'Composition [decimal]',
                'Steps [-]',
                'Job',
                ]

df = pd.merge(dftg, dfcx, on=mergecolumns)

# Alphabetically sort the dataframe
df = df.sort_values(
                    by=[
                        'System',
                        'Composition [decimal]',
                        'Steps [-]',
                        'Job'
                        ]
                    )

df = df.reset_index(drop=True)

# Create the path to work in
workdir = join(sys.argv[1], 'analysis_data')
if not os.path.exists(workdir):
    os.makedirs(workdir)

df.to_html(join(sys.argv[1], 'analysis_data/alltg.html'), index=False)
df.to_csv(join(sys.argv[1], 'analysis_data/alltg.txt'), index=False)

# Filter by data that has not crystallized
df = df[df['Crystallization'] == False]
df = df[columns]  # Remove crystallization column
df = df.loc[:, df.columns != 'Job']  # Remove job column

dfmean = df.groupby(columns[:-3]).agg([np.average])
dfsem = df.groupby(columns[:-3]).agg([st.sem])

df = pd.merge(dfmean, dfsem, how='inner', on=mergecolumns[:-1])
df = pd.DataFrame(df.to_records())

df.to_html(join(sys.argv[1], 'analysis_data/meantg.html'), index=False)
df.to_csv(join(sys.argv[1], 'analysis_data/meantg.txt'), index=False)
