import pymatgen as mg
import pandas as pd
import numpy as np

import traj
import test
import dep

import os

def volume_sphere(r):
    '''
    Calculate the atomic volume from a radius. This assumes that the
    volume is spherical.
    '''

    return 4./3.*np.pi*r**3.


class job:
    '''
    Setup all the data per job for analysis.
    '''

    def __init__(self, path):
        '''
        Load all data needed from a job for analysis.
        '''

        # The location of the job
        self.path = path

        print('Analysis for: '+path)

        # The name of important files for each job
        trajdotlammpstrj = 'traj.lammpstrj'
        testdotout = 'test.out'
        depdotin = 'dep.in'

        # The name of saving directories
        self.datadirname = 'analysis_data'
        self.plotdirname = 'analysis_plots'

        # Location of files
        file_system = os.path.join(path, testdotout)  # Output file
        file_trajs = os.path.join(path, trajdotlammpstrj)  # Trajectories
        file_dep = os.path.join(path, depdotin)  # Input file

        # Important paramters from the input file
        depparams = dep.info(file_dep)

        # Information from traj.lammpstrj file
        self.dftraj, counts = traj.info(file_trajs)

        # Thermodynamic data from test.out file
        self.dfsys = test.info(file_system)

        self.timestep = depparams['timestep']
        self.hold1 = depparams['hold1']
        self.hold2 = depparams['hold2']
        self.hold3 = depparams['hold3']
        self.increment = depparams['increment']
        self.iterations = depparams['iterations']
        self.tempstart = depparams['tempstart']
        self.deltatemp = depparams['deltatemp']
        self.elements = depparams['elements']

        # Add the times for each dataframe
        self.dfsys['time'] = self.dfsys['Step']*self.timestep
        self.dftraj['time'] = self.dftraj['Step']*self.timestep

        # Element counts from actual element
        elements = {}
        for key, count in counts.items():
            element = mg.Element(self.elements[key])
            atomicradii = element.atomic_radius  # in Am
            atomicvol = volume_sphere(atomicradii)  # in Am^3

            elements[self.elements[key]] = {
                                            'counts': count,
                                            'radius': atomicradii,
                                            'volume': atomicvol,
                                            }

        dfelprops = pd.DataFrame(elements).T

        self.dfelprops = dfelprops  # The properties of elements

    def volume(self):
        '''
        Calculate the volume from box dimensions.
        '''

        print('Calculating volume')

        df = pd.DataFrame()

        # Find the box lengths
        df['dx'] = self.dftraj['xhi']-self.dftraj['xlo']  # Am
        df['dy'] = self.dftraj['yhi']-self.dftraj['ylo']  # Am
        df['dz'] = self.dftraj['zhi']-self.dftraj['zlo']  # Am

        self.dfvol = pd.DataFrame()

        # Find the box volume
        self.dfvol['time'] = self.dftraj['time']
        self.dfvol['Volume'] = df['dx']*df['dy']*df['dz']  # Am^3

        return self.dfvol

    def apd(self):
        '''
        Calculate the atomic packing density.
        '''

        print('Calculating APD')

        try:
            self.dfvol
        except Exception:
            job.volume(self)

        # Merge dfsys and dftraj on matching steps
        df = pd.merge(self.dfsys, self.dfvol, on=['time'])

        # Calculate the atomic packing density (APD)
        numerator = np.sum(self.dfelprops['counts']*self.dfelprops['volume'])

        # Create a dataframe for APD
        self.dfapd = pd.DataFrame()
        self.dfapd['time'] = self.dfvol['time']
        self.dfapd['APD'] = numerator/self.dfvol['Volume']

        return self.dfapd

    def save_data(self, system=True, box=True, apd=True):
        '''
        Save all data as csv.

        inputs:
            system = True for saving dataframe
            box = True for saving dataframe
            apd = True for saving dataframe
        '''

        print('Saving data')

        # Create a directory for the analysis files
        savepath = os.path.join(self.path, self.datadirname)
        if not os.path.exists(savepath):
            os.makedirs(savepath)

        # Save system data
        if system:
            self.dfsys.to_csv(
                              os.path.join(savepath, 'system.txt'),
                              index=False
                              )

        # Save information of simulation box
        if box:
            self.dftraj.to_csv(
                               os.path.join(savepath, 'boxboundary.txt'),
                               index=False
                               )

        # Save APD data
        if apd:
            self.dfapd.to_csv(
                              os.path.join(savepath, 'apd.txt'),
                              index=False
                              )
