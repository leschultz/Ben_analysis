from PyQt5 import QtGui  # Added to be able to import ovito

from matplotlib import pyplot as pl

from scipy.signal import argrelextrema
import scipy as sc

import pymatgen as mg
import pandas as pd
import numpy as np

import os

from ovito.modifiers import VoronoiAnalysisModifier
from ovito.io import import_file

import traj
import dep


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

        self.timestep = depparams['timestep']
        self.runsteps = depparams['runsteps']
        self.hold1 = depparams['hold1']
        self.hold2 = depparams['hold2']
        self.hold3 = depparams['hold3']
        self.increment = depparams['increment']
        self.elements = depparams['elements']

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
        self.dftraj['frame'] = frames

    def vp_variance(
                self,
                threshold=0.1,
                write=True,
                plot=True,
                first=10,
                verbose=True
                ):
        '''
        Calculate the maximum variance of clusters.

        inputs:
            self = The object reference
            threshold = The maximum length for a VP edge
            write = Whether or not to save the fractions and temperatures
            plot = Whether or not to plot the fractions and temperatures
            verbose = Wheter or not to print calculation status

        outputs:
            variance = The variance of clusters as a function of number
                       included
        '''

        if verbose:
            print(
                  'Calculating VP variety and ' +
                  'variance at constant temperature hold'
                  )

        try:
            self.file_trajs

        except Exception:
            message = 'Need to specify trajectory file.'
            raise ValueError(message)

        try:
            self.file_dep

        except Exception:
            message = 'Need to specify input file.'
            raise ValueError(message)

        # Find the interval for the isothermal hold
        cutoff1 = sum(self.runsteps[:3])
        cutoff2 = sum(self.runsteps[:4])

        condition = (self.dftraj['Step'] >= cutoff1)
        condition = condition & (self.dftraj['Step'] <= cutoff2)

        # Grab trajectory information from interval
        df = self.dftraj[condition]
        df = df.reset_index(drop=True)

        # Reset time
        df['time'] = df['time']-df['time'][0]
        df['frame'] = df['frame']-df['frame'][0]

        # Load input data and create an ObjectNode with a data pipeline.
        node = import_file(self.file_trajs, multiple_frames=True)

        voro = VoronoiAnalysisModifier(
                                       compute_indices=True,
                                       use_radii=False,
                                       edge_threshold=threshold
                                       )

        node.modifiers.append(voro)

        all_indexes = []
        frames = 0
        for frame in df['frame']:
            out = node.compute(frame)
            indexes = out.particle_properties['Voronoi Index'].array
            all_indexes.append(indexes)
            frames += 1

        # Combine all the frames
        all_indexes = [pd.DataFrame(i) for i in all_indexes]
        df = pd.concat(all_indexes)
        df = df.fillna(0)  # Make cure indexes are zero if not included
        df = df.astype(int)  # Make sure all counts are integers

        # Count the number of unique VP
        coords, counts = np.unique(df.values, axis=0, return_counts=True)

        # Standard notation
        coords = coords[:, 2:]
        coords = np.array([tuple(np.trim_zeros(i)) for i in coords])

        # Sort counts and indexes by descending order
        indexes = counts.argsort()[::-1]
        coords = coords[indexes]
        counts = counts[indexes]

        total = sum(counts)  # The total number of atoms for all frames

        # Gather fraction values
        fractions = counts/total

        # Calculate variance from list including ordered fractions of VP
        variance = []
        for i in range(1, counts.shape[0]):
            variance.append(np.var(fractions[:i]))

        variance = np.array(variance)  # Numpy array
        vptypes = np.arange(1, len(variance)+1)  # Types of vp

        df_variance = {'number_of_vp': vptypes, 'variance': variance}
        df_variance = pd.DataFrame(df_variance)
        df_variance['number_of_vp'] = df_variance['number_of_vp'].astype(int)

        max_index = np.argmax(variance)
        max_variance = pd.DataFrame(df_variance.loc[max_index, :]).T
        max_variance['number_of_vp'] = max_variance['number_of_vp'].astype(int)

        # Create a directory for the analysis files
        if write:

            # Export all variances calculated
            df_variance.to_csv(
                               os.path.join(self.datapath, 'variance.txt'),
                               index=False
                               )

            max_variance.to_csv(
                                os.path.join(self.datapath, 'max_variance.txt'),
                                index=False
                                )


        if plot:

            # Plot variance
            for i in [counts.shape[0], first]:
                fig, ax = pl.subplots()

                ax.plot(
                        vptypes[:i],
                        variance[:i],
                        marker='.',
                        linestyle='none',
                        label='Data for '+str(frames)+' frames'
                        )

                coordinate = (
                              max_variance['number_of_vp'].values[0],
                              max_variance['variance'].values[0]
                              )

                vlinelabel = (
                              'Max Variance: ' +
                              str(coordinate)
                              )
                ax.axvline(
                           max_variance['number_of_vp'].values[0],
                           color='k',
                           linestyle=':',
                           label=vlinelabel
                           )

                ax.set_ylabel('Variance of VP [-]')
                ax.set_xlabel('Number of VP Types [-]')

                ax.grid()
                ax.legend()

                fig.tight_layout()

                plotname = (
                            'variance_top_' +
                            str(i) +
                            '_from_' +
                            str(counts.shape[0]) +
                            '.png'
                            )

                fig.savefig(os.path.join(self.plotpath, plotname))

                pl.close('all')

            # Plot distribution of cluster types seen
            fig, ax = pl.subplots()

            x = [str(i) for i in coords[:first]][::-1]
            y = fractions[:first][::-1]

            ax.barh(x, y)

            ax.set_xlabel('Fraction of VP [-]')
            ax.set_ylabel('The '+str(first)+' Most Frequent VP [-]')

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'variety.png'))

            pl.close('all')

        return variance
