from PyQt5 import QtGui  # Added to be able to import ovito

from matplotlib import pyplot as pl

from scipy.stats import linregress

import pymatgen as mg
import pandas as pd
import scipy as sc
import numpy as np

import ast
import os

from ovito.modifiers import CalculateDisplacementsModifier
from ovito.modifiers import PythonScriptModifier
from ovito.io import import_file

import traj
import test
import dep


def self_diffusion(x, y):
    '''
    Calculate self diffusion from MSD curve.

    inputs:
        x = time
        y = MSD
    outputs:
        d = diffusion coefficient [*10^-4 cm^2 s^-1]
    '''

    m, _, _, _, _ = linregress(x, y)  # Fit linear line
    d = m/6.0  # Divide by degrees of freedom

    return d


def msdmodify(frame, data):
    '''
    Access the per-particle displacement magnitudes computed by an existing
    Displacement Vectors modifier that precedes this custom modifier in the
    data pipeline. This loops over for all and each particle type.

    inputs:
        frame = the frame for trajectories considered
        data = a pipeline variable for data
    '''

    elements = [i.id for i in data.particles['Particle Type'].types]

    # Calculate diplacements
    dispmag = data.particles['Displacement Magnitude'].array

    # Compute MSD for all atoms
    msd = {'all': np.sum(dispmag**2)/len(dispmag)}

    # Compute MSD for a type of atom
    for item in elements:
        index = (data.particles['Particle Type'] == item)
        msd[item] = np.sum(dispmag[index]**2)/len(dispmag[index])

    # Export msd
    data.attributes['msd'] = msd

def gather_msd(file_trajs, start, stop):
    '''
    Calculate MSD for all and each element type.

    inputs:
        file_trajs = name of the trajectory file with the path
        start = the starting frame
        stop = the stopping frame
    outputs:
        dfmsd = dataframe for msd
    '''

    # Load input data and create an ObjectNode with a data pipeline.
    node = import_file(file_trajs, multiple_frames=True)

    # Calculate per-particle displacements with respect to a start
    modifier = CalculateDisplacementsModifier()
    modifier.assume_unwrapped_coordinates = True
    modifier.reference.load(file_trajs)
    modifier.reference_frame = start
    node.modifiers.append(modifier)

    # Insert custom modifier into the data pipeline.
    node.modifiers.append(PythonScriptModifier(function=msdmodify))

    # The variables where data will be held
    msd = []

    # Compute the MSD for each frame of interest
    for frame in range(start, stop+1):

        out = node.compute(frame)
        msd.append(ast.literal_eval(out.attributes['msd']))

    dfmsd = pd.DataFrame(msd)

    return dfmsd


class job:
    '''
    Setup all the data per job for analysis.
    '''

    def __init__(self, path, export, data_path=False, plot_path=False):
        '''
        Create all the paths needed to save analysis data.

        inputs:
            self = The object reference
            path = The path to the data
            export = The path to export analysis
            data_path = The name of the folder to save analysis data
            plot_path = The name off the folder to save analysis plots
        '''

        # The location of the job
        self.path = path

        # Save paths
        export_path = os.path.join(export, path.strip('../'))

        if data_path:
            self.datapath = os.path.join(export_path, data_path)
            if not os.path.exists(self.datapath):
                os.makedirs(self.datapath)

        if plot_path:
            self.plotpath = os.path.join(export_path, plot_path)
            if not os.path.exists(self.plotpath):
                os.makedirs(self.plotpath)

        print('Analysis for: '+path)

    def input_file(self, depdotin):
        '''
        Gather parameters from the input file.

        inputs:
            self = The object reference
            depdotin = The name of the input file
        '''

        file_dep = os.path.join(self.path, depdotin)  # Input file
        self.file_dep = file_dep

        # Important paramters from the input file
        depparams = dep.info(file_dep)

        self.trajdumprate = depparams['trajdumprate']
        self.timestep = depparams['timestep']
        self.runsteps = depparams['runsteps']
        self.hold1 = depparams['hold1']
        self.hold2 = depparams['hold2']
        self.hold3 = depparams['hold3']
        self.increment = depparams['increment']
        self.elements = depparams['elements']

    def sys(self, testdotout):
        '''
        Gather thermodynamic data from the out file.

        inputs:
            self = The object reference
            testdotout = The name of the out file
        '''

        try:
            self.file_dep
        except Exception:

            message = 'Need to specify input file first.'
            raise ValueError(message)

        file_system = os.path.join(self.path, testdotout)  # Output file
        self.file_system = file_system

        # Thermodynamic data from test.out file
        self.dfsys = test.info(file_system)

        self.dfsys['time'] = self.dfsys['Step']*self.timestep

        # Use data at and after start of cooling
        condition = self.dfsys['Step'] >= self.hold1
        dfcool = self.dfsys[condition]

        return dfcool

    def box(self, trajdotlammpstrj):
        '''
        Gather trajectories from the trajectory file.

        inputs:
            self = The object reference
            trajdotlammpstrj = The name of the input file
        '''

        try:
            self.file_dep

        except Exception:
            message = 'Need to specify input file.'
            raise ValueError(message)

        file_trajs = os.path.join(self.path, trajdotlammpstrj)  # Trajectories
        self.file_trajs = file_trajs

        # Information from traj.lammpstrj file
        self.dftraj, counts = traj.info(file_trajs)

        self.dftraj['time'] = self.dftraj['Step']*self.timestep

        # The number of frames where temperatures where recorded
        frames = list(range(self.dftraj.shape[0]))
        self.dftraj['frame'] = np.array(frames)+1

    def msd(self, write=True, plot=True, verbose=True):
        '''
        Calculate MSD.

        inputs:
            self = The object reference
            write = Whether or not to save the fractions and temperatures
            plot = Whether or not to plot the fractions and temperatures
            verbose = Wheter or not to print calculation status

        outputs:
            msd = Dataframe containing msd information
        '''

        if verbose:
            print('Calculating MSD')

        # Find the interval for the isothermal hold
        cutoff = sum(self.runsteps[:5])
        condition = (self.dftraj['Step'] >= cutoff)

        # Grab trajectory information from interval
        df = self.dftraj[condition]
        df = df.reset_index(drop=True)

        # The beggining frame
        frames = df['frame'].values
        start = frames[0]
        stop = frames[-1]

        dfmsd = gather_msd(self.file_trajs, start, stop)
        dfmsd.columns = [list(dfmsd.columns)[0]]+self.elements
        
        dfmsd['time'] = df['time']-df['time'][0]

        if write:
            msd_path = os.path.join(self.datapath, 'msd.txt')
            dfmsd.to_csv(msd_path, index=False)

        if plot:

            fig, ax = pl.subplots()

            plotcols = list(dfmsd.columns.difference(['time']))

            x = dfmsd['time'].values
            for col in plotcols:
                y = dfmsd[col].values

                ax.plot(
                        x,
                        y,
                        linestyle='none',
                        marker='.',
                        label='element: '+col
                        )

            ax.grid()
            ax.legend()

            ax.set_xlabel('Time [ps]')
            ax.set_ylabel(r'MSD $[A^{2}]$')

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'msd.png'))

        return msd

    def diffusion(self, write=True, plot=True, verbose=True):
        '''
        Calculate diffusion from multiple time origins (MTO).

        inputs:
            self = The object reference
            write = Whether or not to save the fractions and temperatures
            plot = Whether or not to plot the fractions and temperatures
            verbose = Wheter or not to print calculation status

        outputs:
            diff = The mean diffusion
            dfdiff = The MTO diffusion dataframe
        '''

        if verbose:
            print('Calculating MTO diffusion')

        # Find the interval for the isothermal hold
        cutoff = sum(self.runsteps[:5])
        condition = (self.dftraj['Step'] >= cutoff)

        # Grab trajectory information from interval
        df = self.dftraj[condition]
        df = df.reset_index(drop=True)

        # Reset time
        df['time'] = df['time']-df['time'][0]

        frames = df['frame'].values

        # Split data in half
        cut = frames.shape[0]//2
        split1 = frames[:cut]
        split2 = frames[cut:cut+split1.shape[0]]

        # Each of the time origins
        time_origins = df['time'].values[:cut]
        time_endings = df['time'].values[cut:cut+time_origins.shape[0]]

        # Collect diffusion coefficients
        data = []
        count = 0
        for start, stop in zip(split1, split2):
            dfmsd = gather_msd(self.file_trajs, start, stop)
            dfmsd.columns = [list(dfmsd.columns)[0]]+self.elements

            # Remove first value which is always zero
            dfmsd = dfmsd.loc[1:, :]

            # Calculate diffusion for each element and all
            d = dfmsd.apply(lambda x: self_diffusion(time_origins, x))
            data.append(d)

            count += 1
            if count > 1:
                break

        dfdif = pd.DataFrame(data)
        dfdif['start'] = time_origins[:2]
        dfdif['stop'] = time_endings[:2]

        if write:
            dfdif.to_csv(
                         os.path.join(self.datapath, 'diffusion_mo.txt'),
                         index=False
                         )

        if plot:

            fig, ax = pl.subplots()

            plotcols = list(dfdif.columns.difference(['start', 'stop']))

            x = dfdif['start'].values
            for col in plotcols:
                y = dfdif[col].values

                ax.plot(
                        x,
                        y,
                        linestyle='none',
                        marker='.',
                        label='element: '+col
                        )

            ax.grid()
            ax.legend()

            ax.set_xlabel('Time Origin from '+str(time_endings[-1])+' [ps]')
            ax.set_ylabel(r'MSD $[10^{-4} cm^2 s^-1]$')

            fig.tight_layout()
            pl.show()
            fig.savefig(os.path.join(self.plotpath, 'diffusion_mo.png'))

        return dfdif
