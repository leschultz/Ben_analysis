from scipy import stats as st

from os.path import join

import pandas as pd
import numpy as np

import sys
import os

df1file = 'All_jobs_with_crystal_info_0_2pct.csv'
df2file = 'All_jobs_with_crystal_info_4_10pct.csv'

def load(name):
    '''
    Load and change the crystallization dataframes to a standard.

    inputs:
        name = The file with the path containing crystallization data

    outputs:
        dfcx = The crystallization dataframe
    '''

    dfcx = pd.read_csv(name)

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

    return dfcx

df1 = load(df1file)
df2 = load(df2file)

df = pd.concat([df1, df2])
df.to_csv('crystalinfo.txt', index=False)
