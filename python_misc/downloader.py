#!/usr/bin/env python3

'''
A script to copy data from crystallization runs
'''

from shutil import copyfile
from os.path import join
import tarfile
import sys
import os

source = sys.argv[1]
destination = sys.argv[2]

def iterator(source, destination):
    '''
    Search for all possible files to copy for Tg analysis.

    inputs:
        source = path of original files to be copied
        destination = path to save copy of files

    outputs:
        A collection of data files
    '''

    # Loop for each path
    for item in os.walk(source):
        path = item[0]
        files = item[2]

        condition = ('fcc' in path) or ('bcc' in path) or ('hcp' in path)
        if not condition:
            continue

        if not files:
            continue

        folders = path.split('/')[-3:]
        dump = join(*[destination]+folders)

        sor = list(map(lambda i: join(path, i), files))
        des = list(map(lambda i: join(dump, i), files))

        if not os.path.exists(dump):
            os.makedirs(dump)

        [copyfile(i, j) for i, j in zip(sor, des)]

        print('Copying from: '+path)

iterator(source, destination)
