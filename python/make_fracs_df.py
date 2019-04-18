from scipy import stats as st

from os.path import join

import pandas as pd
import numpy as np

import sys
import os

datadirname = 'analysis_data'
externaldirname = 'cryxdata'

fracsfile = 'fracs.txt'

cxpath = [externaldirname, 'crystalinfo.txt']
dfcx = pd.read_csv(join(*cxpath))

# Build a dataframe for all collected APD values
columns = [
           'System',
           'Composition [decimal]',
           'Steps [-]',
           'Job',
           'max variance',
           'variety'
           ]

dffracs = pd.DataFrame(columns=columns)

count = 0
for item in os.walk(sys.argv[1]):

    if datadirname not in item[0]:
        continue

    if '-' not in item[0]:
        continue

    names = item[0].split('/')[-5:-1]

    fracspath = join(item[0], fracsfile)
    if os.path.exists(fracspath):
        data = pd.read_csv(join(item[0], fracsfile))
        meanvariance = np.mean(data['max variance'])
        meanvariety = np.mean(data['variety'])

    else:
        meanvariance = np.nan
        meanvariety = np.nan

    row = names+[meanvariance, meanvariety]

    dffracs.loc[count] = row  # Append a row to the dataframe

    count += 1

dffracs['Steps [-]'] = dffracs['Steps [-]'].apply(pd.to_numeric)

# The columns to merge crystallization and tg dataframes
mergecolumns = [
                'System',
                'Composition [decimal]',
                'Steps [-]',
                'Job',
                ]

df = pd.merge(dffracs, dfcx, on=mergecolumns)

# Alphabetically sort the dataframe
df = df.sort_values(by=mergecolumns)
df = df.reset_index(drop=True)

# Create the path to work in
workdir = join(sys.argv[1], datadirname)
if not os.path.exists(workdir):
    os.makedirs(workdir)

df.to_html(join(*[sys.argv[1], datadirname, 'allfracs.html']), index=False)
df.to_csv(join(*[sys.argv[1], datadirname, 'allfracs.txt']), index=False)

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
    variance = item[1]['max variance']
    n0 = len(variance)

    mean = np.mean(variance)
    std = np.std(variance)
    filtered = variance[abs(variance-mean) < nsig*std]
    nremoved0 = n0-len(filtered)

    if filtered.size == 0:
        newmean0 = np.nan
        newsem0 = np.nan
    else:
        newmean0 = np.mean(filtered)
        newsem0 = st.sem(filtered)

    variety = item[1]['variety']
    n1 = len(variety)

    mean = np.mean(variety)
    std = np.std(variety)
    filtered = variety[abs(variety-mean) < nsig*std]
    nremoved1 = n1-len(filtered)

    if filtered.size == 0:
        newmean1 = np.nan
        newsem1 = np.nan
    else:
        newmean1 = np.mean(filtered)
        newsem1 = st.sem(filtered)

    row = [
           system,
           comp,
           step,
           newmean0,
           newsem0,
           n0,
           nremoved0,
           newmean1,
           newsem1,
           n1,
           nremoved1
           ]

    dfmean.append(row)

meancolumns = mergecolumns[:-1]
meancolumns += [
                'mean max variance',
                'sem max variance',
                'jobs max variance',
                'jobs outside '+str(nsig)+' sigma for max variance'
                ]

meancolumns += [
                'mean variety',
                'sem variety',
                'jobs variety',
                'jobs outside '+str(nsig)+' sigma for variety'
                ]

df = pd.DataFrame(dfmean, columns=meancolumns)

df.to_html(join(*[sys.argv[1], datadirname, 'meanfracs.html']), index=False)
df.to_csv(join(*[sys.argv[1], datadirname, 'meanfracs.txt']), index=False)
