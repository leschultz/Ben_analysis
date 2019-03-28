import pandas as pd
import numpy as np

datadirname = 'analysis_data'
externaldirname = 'data_external'

dfcx = pd.read_csv('../../'+externaldirname+'/All_jobs_with_crystal_info.csv')
dftg = pd.read_csv('../../'+datadirname+'/Tg.txt')

# Modify Ben's columns for mergin datasets
alloys = dfcx['alloy']
elements = [i.split('-')[1] for i in alloys]
comp = dfcx['Comp']
comp = ['{:.2f}'.format(np.round(i, 2)) for i in comp]
comp = [i+j for i, j in zip(elements, comp)]
dfcx['Comp'] = comp

df = dfcx[['alloy', 'Comp', 'Thold', 'Job', 'Run Crystallized']]
df.columns = [
              'System',
              'Composition [decimal]',
              'Steps [-]',
              'Job',
              'Crystallization'
              ]

mergecolumns = [
                'System',
                'Composition [decimal]',
                'Steps [-]',
                'Job',
                ]

dfmerge = pd.merge(dftg, df, on=mergecolumns)

dfmerge.to_html('../../'+datadirname+'/tg_and_crystallization.html', index=False)
dfmerge.to_csv('../../'+datadirname+'/tg_and_crystallization.txt', index=False)
