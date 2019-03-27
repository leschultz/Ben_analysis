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
        counts = The number of types of atoms
    '''

    steps = []

    xlo = []
    xhi = []

    ylo = []
    yhi = []

    zlo = []
    zhi = []

    # Gather steps and simulation box boundaries
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

    # Gather the number of elements
    counts = {}  # The counts for the types of atoms
    with open(name) as file:
        terminatecount = -1
        startcount = 0
        for line in file:

            # Break the loop after one frame
            if terminatecount == 1:
                break

            if (startcount == 1) & ('TIMESTEP' not in line):
                values = line.split(' ')
                element_type = values[1]

                if counts.get(element_type) is None:
                    counts[element_type] = 1

                else:
                    counts[element_type] += 1

            if 'TIMESTEP' in line:
                terminatecount += 1

            if 'ATOMS id type' in line:
                startcount += 1

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

    return df, counts
