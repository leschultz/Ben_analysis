from shutil import copy, copytree

import tarfile
import sys
import os


# The name of important files for each job
trajdotlammpstrj = 'traj.lammpstrj'
testdotout = 'test.out'
depdotin = 'dep.in'
compressed = 'outputs.tar.gz'


def filecopy(copyname, filelist, copypath, savepath):
    '''
    Copy the file if it exists in path.

    inputs:
        copyname = The name of the file of interest
        filelist = The list of files available
        copypath = The path from where the file is copied from
        savepath = The path to save the file

    outputs:
        A copied file in a directory
    '''

    # Grab the output if the file exists
    if (copyname in filelist):
        exist = copyname
        filename = os.path.join(copypath, copyname)

        # Create the path to work in
        if not os.path.exists(savepath):
            os.makedirs(savepath)

        copy(filename, savepath)


def copydata(item, outdir):
    '''
    Copy the data.

    inputs:
        item = The output from a os.walk generator
        outdir = The name for the export path

    outputs:
        A collection of data files
    '''

    if ('job' in item[0]) & ('minimization' not in item[0]):

        # Were parsed data will be stored
        data = []

        # Create a name from the path
        name = item[0].split('/')
        name = name[-4:]
        name = os.path.join(*name)

        # Print status
        print('Copying data: '+name)

        savepath = os.path.join(*[outdir, name])

        # Grab the output archive file that contains run system data
        if (compressed in item[-1]) & (depdotin in item[-1]):
            filename = os.path.join(*[item[0], compressed])
            inputfile = os.path.join(*[item[0], depdotin])

            # Open the archive
            archive = tarfile.open(filename, 'r')

            # Iterate for each file in the archive
            for member in archive.getmembers():

                # Open the file containing system data
                if testdotout in str(member):
                    member = archive.extractfile(member)
                    content = member.read()

                # Open the file containing the trajectory data
                if trajdotlammpstrj in str(member):
                    member = archive.extractfile(member)
                    traj = member.read()

            # Create the path to work in
            if not os.path.exists(savepath):
                os.makedirs(savepath)

            # Save the information up one directory
            copy(inputfile, savepath)
            open(os.path.join(savepath, testdotout), 'wb').write(content)
            open(os.path.join(savepath, trajdotlammpstrj), 'wb').write(traj)

        # Grab the files if they exist
        else:

            if testdotout in item[-1]:
                filecopy(testdotout, item[-1], item[0], savepath)
            if depdotin in item[-1]:
                filecopy(depdotin, item[-1], item[0], savepath)
            if trajdotlammpstrj in item[-1]:
                filecopy(trajdotlammpstrj, item[-1], item[0], savepath)

    if 'minimization' in item[0]:

        # Create a name from the path
        name = item[0].split('/')
        name = name[-5:]
        name = os.path.join(*name)

        savepath = os.path.join(*[outdir, name])

        copytree(item[0], savepath)


def jobiterator(path, outdir):
    '''
    Search for all possible files to copy for Tg analysis.

    inputs:
        path = The path with all the runs
        outdir = The export path

    outputs:
        A collection of data files
    '''

    # Look for all directories as generator object
    paths = os.walk(path)

    # Loop for each path
    for item in paths:
        copydata(item, outdir)


# Download data into current directory
jobiterator(sys.argv[1], sys.argv[2])
