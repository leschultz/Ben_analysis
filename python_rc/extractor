#!/usr/bin/env python3

'''
A script to extract files from Ben's runs
'''

from shutil import copyfile
from os.path import join
import tarfile
import sys
import os


def iterator(path, outdir, compressed, folders, files, depth):
    '''
    Search for all possible files to copy for Tg analysis.

    inputs:
        path = The path with all the runs
        outdir = The export path
        compressed = List of compressed files
        files = List of needed files
        depth = The depth for directory naming convention

    outputs:
        A collection of data files
    '''

    # Look for all directories as generator object
    paths = os.walk(path)
    extensions = [i.split('.')[-1] for i in files]

    # Loop for each path
    for item in paths:
        for tar in compressed:
            if tar not in item[2]:
                continue

            archive = tarfile.open(join(item[0], tar), 'r')
            for member in archive.getmembers():
                if member.name not in files:
                    continue

                name = item[0].split('/')
                name = name[-depth:]
                name = join(outdir, *name)

                # Create the path
                if not os.path.exists(name):
                    os.makedirs(name)

                name = join(name, member.name)

                print('Creating: '+name)
                content = archive.extractfile(member).read()
                open(name, 'wb').write(content)

        for folder in folders:
            if folder not in item[0]:
                continue

            for f in item[2]:
                if f not in files:
                    continue

                name = item[0].split('/')
                name = name[-depth-1:]+item[1]
                name = join(outdir, *name)

                # Create the path
                if not os.path.exists(name):
                    os.makedirs(name)

                name = join(name, f)

                print('Creating: '+name)
                copyfile(join(item[0], f), name)


# The name of important files for each job
compressed = [
              'inputs.tar.gz',
              'outputs.tar.gz'
              ]
folders = [
           '100K_Structure_minimization'
           ]
files = [
         'traj.lammpstrj',
         'test.out',
         './dep.in',
         'finaltraj.lammpstrj',
         '100k_minimize_template.in'
         ]

# Extract data
iterator(sys.argv[1], sys.argv[2], compressed, folders, files, 5)
