from matplotlib import pyplot as pl

import pandas as pd
import numpy as np

import sys
import os
import re

dtemp = float(sys.argv[1])  # The change in temperature [K]
timestep = float(sys.argv[2])  # The timestep [ps]

dataframe = sys.argv[3]  # The name and path of the dataframe
savefolder = sys.argv[4]  # The name of the save folder for the images

if not os.path.exists(savefolder):
    os.makedirs(savefolder)

df = pd.read_csv(dataframe)
group = df.groupby(['system', 'steps'])

for i in group:

    data = i[1]

    system = i[0][0]
    elements = system.split('-')

    rate = dtemp/(i[0][1]*timestep*1e-12)

    y = data['tg_mean'].values
    yerr = data['tg_sem'].values

    x = data['composition']

    x = [float(re.findall("\d+\.\d+", i)[0]) for i in x]

    fig, ax = pl.subplots()

    ax.errorbar(
                x,
                y,
                yerr,
                marker='.',
                ecolor='r',
                linestyle='none',
                label='Cooling Rate: '+"{:.2E}".format(rate)+' [K/s]'
                )

    ax.set_xlabel('Fraction of '+elements[1]+' in '+elements[0])
    ax.set_ylabel('Tg [K]')

    ax.grid()
    ax.legend(loc='best')

    fig.tight_layout()

    name = i[0][0]+'_'+str(i[0][1])

    fig.savefig(savefolder+'composition_vs_tg_'+name)
