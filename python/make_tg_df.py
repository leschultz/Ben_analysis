from scipy import stats as st

from os.path import join

import pandas as pd
import numpy as np

import sys
import os

datadirname = 'analysis_data'
externaldirname = 'cryxdata'

etgfile = 'tg_e.txt'
etgcutfile = 'tg_e_t_cutoff.txt'

cxpath = [externaldirname, 'crystalinfo.txt']
dfcx = pd.read_csv(join(*cxpath))

# Build a dataframe for all collected Tg values
columns = [
           'System',
           'Composition [decimal]',
           'Steps [-]',
           'Job',
           'Tg from E-3kT Curve [K]',
           'Upper Temperature Cutoff [K]'
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

    ecutpath = join(item[0], etgcutfile)
    if os.path.exists(ecutpath):
        etgcut = np.loadtxt(ecutpath, dtype=float)
    else:
        etgcut = np.nan


    row = names+[etg, etgcut]

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
df = df.sort_values(by=mergecolumns)

df = df.reset_index(drop=True)

# Create the path to work in
workdir = join(sys.argv[1], datadirname)
if not os.path.exists(workdir):
    os.makedirs(workdir)

df.to_html(join(*[sys.argv[1], datadirname, 'alltg.html']), index=False)
df.to_csv(join(*[sys.argv[1], datadirname, 'alltg.txt']), index=False)

# Filter by data that has not crystallized
df = df[df['Crystallization'] == False]
df = df[columns]  # Remove crystallization column
df = df.loc[:, df.columns != 'Job']  # Remove job column

group = df.groupby(columns[:-3])

nsig = 2
dfmean = []

for item in group:
    system = item[0][0]
    comp = item[0][1]
    step = item[0][2]

    etg = item[1]['Tg from E-3kT Curve [K]']
    en = len(etg)

    emean = np.mean(etg)
    estd = np.std(etg)
    efiltered = etg[abs(etg-emean) < nsig*estd]
    enremoved = en-len(efiltered)

    ecut = item[1]['Upper Temperature Cutoff [K]']
    ecutfiltered = ecut[abs(etg-emean) < nsig*estd]

    if efiltered.size == 0:
        enewmean = np.nan
        enewsem = np.nan
        ecutnewmean = np.nan
    else:
        enewmean = np.mean(efiltered)
        enewsem = st.sem(efiltered)
        ecutnewmean = np.mean(ecutfiltered)

    row = [
           system,
           comp,
           step,
           enewmean,
           enewsem,
           ecutnewmean,
           en,
           enremoved,
           ]

    dfmean.append(row)

meancolumns = mergecolumns[:-1]
meancolumns += [
                'Mean Tg from E-3kT Curve [K]',
                'Sem Tg from E-3kT Curve [K]',
                'Mean Upper Temperature Cutoff [K]',
                'Jobs from Tg from E-3kT Curve [K]',
                'Jobs outside '+str(nsig)+' sigma from Tg from E-3kT Curve [K]'
                ]

df = pd.DataFrame(dfmean, columns=meancolumns)

df.to_html(join(*[sys.argv[1], datadirname, 'meantg.html']), index=False)
df.to_csv(join(*[sys.argv[1], datadirname, 'meantg.txt']), index=False)
