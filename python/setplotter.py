from matplotlib import pyplot as pl

from scipy import stats

import pandas as pd
import numpy as np

import sys
import os

# The name of the file containing Tg
tg_file = 'tg_e.txt'

# Simulation paramters
timestep = float(sys.argv[1])  # The time step in [ps]
dt = float(sys.argv[2])  # The change in temperature in [K]

set_path = sys.argv[3]  # Path for set of jobs
data_path = sys.argv[4]  # Path to save analysis data
plot_path = sys.argv[5]  # Path to save analysis plots

df = []
for item in os.walk(set_path):

    # Gather Tg values
    if tg_file not in item[2]:
        continue

    path = item[0]

    splits = path.split('/')
    job = splits[-2]
    steps = int(splits[-3])

    tg = float(np.loadtxt(os.path.join(path, tg_file)))

    df.append([steps, job, tg])

cols = ['steps', 'job', 'tg']
df = pd.DataFrame(df, columns=cols)

df['rate'] = dt/(df['steps']*timestep*10**-12)

dfmean = df.groupby([cols[0]], as_index=False).mean()
dfsem = df.groupby([cols[0]], as_index=False).sem()

# Create a directory for the analysis files
if not os.path.exists(data_path):
    os.makedirs(data_path)

if not os.path.exists(plot_path):
    os.makedirs(plot_path)

dfmean.to_csv(os.path.join(data_path, 'tg_mean.txt'), index=False)

x = dfmean['rate'].values
y = dfmean['tg'].values
yerr = dfsem['tg'].values

fig, ax = pl.subplots()

ax.errorbar(
            x,
            y,
            yerr,
            ecolor='r',
            linestyle='none',
            marker='.',
            label='Mean Between Jobs'
            )

ax.set_xlabel('Cooling Rate [K/s]')
ax.set_ylabel('Tg [K]')

ax.set_xscale('log')

ax.legend()
ax.grid()

fig.tight_layout()
fig.savefig(os.path.join(plot_path, 'tg_mean'))
