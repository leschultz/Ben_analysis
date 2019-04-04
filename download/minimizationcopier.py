from shutil import copytree

import os


def jobiterator(path, savedir):
    '''
    Search for all possible files to copy for Tg analysis.

    inputs:
        path = The path with all the runs
        name = The name fo the save directory

    outputs:
        A collection of data files
    '''

    # Look for all directories as generator object
    paths = os.walk(path)

    # Loop for each path
    for item in paths:

        if 'minimization' not in item[0]:
            continue

        name = item[0].replace(path, '')
        print(name)
        save = '../../.'+name
        copytree(item[0], save)


# Download data into current directory
path = '/home/nerve/Documents/UW/gdrive/DMREF/md/Rc_database/comp_0pct_2pct'
savedir = '../../.'

if not os.path.exists(savedir):
    os.makedirs(savedir)

jobiterator(path, savedir)
