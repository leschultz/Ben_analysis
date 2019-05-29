#!/usr/bin/env python3

from shutil import copyfile
from os.path import join

import sys
import os

datadir = sys.argv[1]  # Directory containing runs
jobname = sys.argv[2]  # Generic job name
template = sys.argv[3]  # Template file path
potentialsdir = sys.argv[4]  # Where potentials are stored

# Gather poteantials
potentials = {}
for path, subdirs, files in os.walk(potentialsdir):
    for f in files:
        split = f.split('.')
        system = split[0]

        if ('eam' and 'fs') in f:
            filetype = 'eam/fs'

        else:
            filetype = 'eam/alloy'

        potentials[system] = {'file': f, 'filetype': filetype}

# Create generator object to iterate through jobs
for path, subdirs, files in os.walk(datadir):

    split = path.split('/')

    if '2450k_minimization' != split[-2]:
        continue

    if jobname not in path:
        continue

    system = split[-6].replace('-', '')

    elements = [i for i in split[-6].split('-')]
    elements = ' '.join(elements)
    timestep = split[-1]
    potfile = potentials[system]['file']
    filetype = potentials[system]['filetype']

    print(path)

    parentpot = len(split)*'../'+'potentials'
    minfile = join(path, 'minimization.in')

    with open(template, 'rt') as i:
        with open(minfile, 'wt') as j:
            for line in i:
                if 'TIMESTEP' in line:
                    j.write(line.replace('TIMESTEP', timestep))
                elif 'POTENTIAL' in line:
                    replacement = join(parentpot, potfile)+' '+elements
                    j.write(line.replace('POTENTIAL', replacement))

                elif 'STYLE' in line:
                    j.write(line.replace('STYLE', filetype))

                else:
                    j.write(line)
