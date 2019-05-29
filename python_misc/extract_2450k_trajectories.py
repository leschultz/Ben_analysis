#!/usr/bin/env python3

from os.path import join

import sys
import os

datadir = sys.argv[1]  # Directory containing runs

jobname = sys.argv[2]  # The generic job name
stepset = sys.argv[3]  # The set of steps to use
inputname = sys.argv[4]  # The name for the LAMMPS input file
trajname = sys.argv[5]  # The name for the trajectories file

exportdir = sys.argv[6]  # Directory to save trajectories

# Create generator object to iterate through jobs
for path, subdirs, files in os.walk(datadir):

    split = path.split('/')
    if jobname not in split[-1]:
        continue

    if stepset != split[-2]:
        continue

    export = join(*[exportdir, path.strip(datadir), '2450k_minimization'])
    print(export)

    if not os.path.exists(export):
        os.makedirs(export)

    # Paths to files
    inputfile = join(path, inputname)
    trajfile = join(path, trajname)

    # Gather the steps for when he 2450 K hold occurs
    holds = []
    with open(inputfile) as f:
        for line in f:
            values = line.split(' ')
            values = [i for i in values if i != '']

            if values[0] == 'run':
                holds.append(int(values[1].strip('\n')))

    hold1 = sum(holds[:3])  # Preparation holds
    hold2 = holds[3]  # Hold at 2450 K

    # Parse trajectory file for 2450 K hold
    with open(trajfile) as f:

        case = 0  # Determine what the next line is
        for line in f:
            if 'ITEM: TIMESTEP' in line:
                case = 1
                continue

            if case == 1:
                step = int(line)
                case = 0

                # Condition to cut trajectories
                if (step >= hold1) & (step <= hold2):
                    condition = True

                    name = join(export, 'frame_'+str(step)+'.lammpstrj')

                    # The file to dump data into
                    w = open(name, 'w')
                    w.write('ITEM: TIMESTEP\n')

                else:
                    condition = False

            if not condition:
                continue

            w.write(line)

        w.close()
