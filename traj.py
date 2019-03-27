'''
Import the steps that are printed onto a trajectory file.
'''

import pandas as pd


def info(name):
    '''
    Gather the steps where trajectories were dumped and the volume dimensions.

    inputs:
        name = The location of the trajectory file

    outputs:
        df = The parsed parameters of interest
    '''

    steps = []

    xlo = []
    xhi = []

    ylo = []
    yhi = []

    zlo = []
    zhi = []

    with open(name) as file:
        for line in file:
            if 'TIMESTEP' in line:
                line = next(file)
                line = int(line.strip('\n'))
                steps.append(line)

                # Skip lines to get box dimensions
                line = next(file)
                line = next(file)
                line = next(file)

                # x coordinates
                line = next(file)
                x = line.split(' ')
                xlo.append(float(x[0]))
                xhi.append(float(x[1]))

                # y coordinates
                line = next(file)
                y = line.split(' ')
                ylo.append(float(y[0]))
                yhi.append(float(y[1]))

                # z coordinates
                line = next(file)
                z = line.split(' ')
                zlo.append(float(z[0]))
                zhi.append(float(z[1]))

    param = {
             'Step': steps,
             'xlo': xlo,
             'xhi': xhi,
             'ylo': ylo,
             'yhi': yhi,
             'zlo': zlo,
             'zhi': zhi,
             }

    df = pd.DataFrame(param)

    return df
