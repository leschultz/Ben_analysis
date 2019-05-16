from PyQt5 import QtGui  # Added to be able to import ovito

from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

from matplotlib import pyplot as pl

from scipy.signal import argrelextrema
import scipy as sc

import pymatgen as mg
import pandas as pd
import numpy as np

import os

from ovito.modifiers import VoronoiAnalysisModifier
from ovito.io import import_file

from line_intersector import opt

import traj
import test
import dep


def volume_sphere(r):
    '''
    Calculate of a sphere given a radius.
    '''

    return 4.0/3.0*sc.pi*r**3.0


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
        self.iterations = depparams['iterations']
        self.deltatemp = depparams['deltatemp']
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
        self.dftraj['frame'] = frames

        # Element counts from actual element
        self.natoms = 0  # Count the total number of atoms
        elements = {}
        for key, count in counts.items():
            self.natoms += count

            element = mg.Element(self.elements[key])

            try:
                atomicradii = element.metallic_radius  # in Am

            except Exception:
                atomicradii = element.atomic_radius  # in Am

            atomicvol = volume_sphere(atomicradii)  # in Am^3

            elements[self.elements[key]] = {
                                            'counts': count,
                                            'radius': atomicradii,
                                            'volume': atomicvol,
                                            }

        dfelprops = pd.DataFrame(elements).T

    def find_tl(
                self,
                edges=5,
                faces=10,
                lim=0.1,
                threshold=0.1,
                min_temp=800,
                write=True,
                plot=True,
                verbose=True
                ):
        '''
        Compute the liquidus temperature based on VP curve.

        inputs:
            self = The object reference
            edges = The number of VP edges
            faces = The number of minimum faces for the specified edges
            lim = The threshold for the upper and lower cut values
            threshold = The maximum length for a VP edge
            write = Whether or not to save the fractions and temperatures
            plot = Whether or not to plot the fractions and temperatures
            verbose = Wheter or not to print calculation status


        outputs:
            tl = The liquidus temperature
        '''

        if verbose:
            print('Calculating Tl from VP curve')

        edges -= 1  # Compensate for indexing

        # Merge dfsys and dftraj on matching steps
        df = pd.merge(
                      self.dfsys.loc[:, self.dfsys.columns != 'Volume'],
                      self.dftraj, on=['Step', 'time']
                      )

        condition = df['Step'] >= self.hold1
        df = df[condition]
        df = df.reset_index(drop=True)

        # Load input data and create an ObjectNode with a data pipeline.
        node = import_file(self.file_trajs, multiple_frames=True)

        voro = VoronoiAnalysisModifier(
                                       compute_indices=True,
                                       use_radii=False,
                                       edge_threshold=threshold
                                       )

        node.modifiers.append(voro)

        fractions = []
        count = 0
        for frame in df['frame']:
            out = node.compute(frame)
            indexes = out.particle_properties['Voronoi Index'].array

            # Count unique VP types
            coords, counts = np.unique(indexes, axis=0, return_counts=True)

            # Sort counts and indexes by descending order
            indexes_sorted = counts.argsort()[::-1]
            coords = coords[indexes_sorted]
            counts = counts[indexes_sorted]

            if count == 0:
                vp_common_first = coords[0]
                n_common_first = len(vp_common_first)
                count += 1

            n = len(coords[0])
            if n_common_first < n:
                vp_common = np.pad(vp_common_first, (0, n-n_common_first), 'constant')
            elif n_common_first > n:
                vp_common = vp_common_first[:-(n_common_first-n)]
            else:
                vp_common = vp_common_first

            first_index = np.where((coords == vp_common).all(axis=1))
            fraction = counts[first_index]/self.natoms  # Calculate fraction

            fractions.append(fraction)

        df['fractions'] = fractions
        df = df.sort_values(by=['Temp'])

        x = df['Temp'].values
        y = df['fractions'].values

        # Cutoff region
        cut = x >= min_temp
        xcut = x[cut]
        ycut = y[cut]

        # Spline fit of cut region
        k, s = (5, 1)
        spl = UnivariateSpline(x=xcut, y=ycut, k=k, s=s)
        xfitcut = np.linspace(xcut[0], xcut[-1], 100)
        yfitcut = spl(xfitcut)

        tl, endpoints, middle_rmse = opt(xfitcut, yfitcut)

        # Standard notation
        vp_tracked = vp_common_first[2:]
        vp_tracked = tuple(np.trim_zeros(vp_tracked))

        if plot:

            fig, ax = pl.subplots()

            ax.plot(x, y, marker='.', linestyle='none', label='data')
            ax.axvline(
                       tl,
                       color='r',
                       linestyle=':',
                       label='Tl='+str(tl)+' [K]'
                       )

            ax.axvline(
                       min_temp,
                       color='k',
                       label='cutoff = '+str(min_temp)+' [K]'
                       )

            ax.plot(
                    xfitcut,
                    yfitcut,
                    linestyle='-',
                    label='Univariate Spline (k='+str(k)+', s='+str(s)+')'
                    )

            xlabel = 'Temperature [K]'
            ylabel = 'Fractions of most common VP [-]: '+str(vp_tracked)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax.grid()
            ax.legend()

            fig.tight_layout()
            pl.show()

        return tl

    def etg(
            self,
            max_temp=1000,
            write=True,
            plot=True,
            verbose=True
            ):
        '''
        Calculate the glass transition temperature based on E-3kt.

        inputs:
            self = The object reference
            max_temp = The maximum temperature for analysis
            write = Whether or not to save Tg
            plot = Whether or not to save plot of data
            verbose = Wheter or not to print calculation status

        outputs:
            tg = The Tg
        '''

        if verbose:
            print('Calculating Tg from E-3kT')

        try:
            self.file_system
        except Exception:
            message = 'Need to specify LAMMPS output file.'
            raise ValueError(message)

        try:
            self.file_trajs
        except Exception:
            message = 'Need to specify trajectory file.'
            raise ValueError(message)

        k = sc.constants.physical_constants['Boltzmann constant in eV/K']

        df = pd.DataFrame()

        e = self.dfsys['TotEng']/self.natoms-3.*k[0]*self.dfsys['Temp']

        df['Temp'] = self.dfsys['Temp']
        df['E-3kT'] = e

        # Use data at and after start of cooling
        condition = self.dfsys['Step'] >= self.hold1
        dfcool = df[condition]
        dfcool = dfcool.sort_values(by=['Temp'])

        x = dfcool['Temp'].values
        y = dfcool['E-3kT'].values

        # Cutoff region
        condition = x <= max_temp
        xcut = x[condition]
        ycut = y[condition]

        # Spline fit of cut region
        k, s = (5, 1)
        spl = UnivariateSpline(x=xcut, y=ycut, k=k, s=s)
        xfitcut = np.linspace(xcut[0], xcut[-1], 100)
        yfitcut = spl(xfitcut)

        tg, endpoints, middle_rmse = opt(xfitcut, yfitcut)

        if write:

            # Export the glass transition temperature
            write_name = os.path.join(self.datapath, 'etg.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(tg))

            # Export the upper temperature cutoff
            write_name = os.path.join(self.datapath, 'etg_temp_cutoff.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(max_temp))

        if plot:

            fig, ax = pl.subplots(2)

            ax[0].plot(
                       x,
                       y,
                       marker='.',
                       linestyle='none',
                       color='b',
                       label='data'
                       )

            ax[0].axvline(
                          max_temp,
                          linestyle='--',
                          color='r',
                          label='cutoff = '+str(max_temp)+' [K]'
                          )

            ax[0].set_ylabel('E-3kT [K/atom]')
            ax[0].grid()
            ax[0].legend()

            ax[1].plot(
                       xcut,
                       ycut,
                       marker='.',
                       linestyle='none',
                       color='b',
                       label='data'
                       )

            ax[1].plot(
                       xfitcut,
                       yfitcut,
                       linestyle='-',
                       label='Univariate Spline (k='+str(k)+', s='+str(s)+')'
                       )

            ax[1].axvline(
                          tg,
                          linestyle='--',
                          color='r',
                          label='Tg = '+str(tg)+' [K]'
                          )

            ax[1].grid()
            ax[1].legend()

            ax[1].set_xlabel('Temperature [K]')
            ax[1].set_ylabel('E-3kT [K/atom]')

            fig.set_size_inches(15, 10, forward=True)
            fig.savefig(os.path.join(self.plotpath, 'etg.png'))
            pl.close('all')

            fig, ax = pl.subplots()

            ax.plot(endpoints, middle_rmse, label='Mean RMSE')

            ax.axvline(
                       tg,
                       linestyle='--',
                       color='r',
                       label='Tg = '+str(tg)+' [K]'
                       )

            ax.set_xlabel('End Temperature [K]')
            ax.set_ylabel('RMSE')
            ax.grid()
            ax.legend()

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'etg_rmse.png'))
            pl.close('all')

        return tg

    def apd_last(self, traj_path, in_path, write=True, verbose=True):
        '''
        Calculate the APD from the last trajectories.

        inputs:
            traj_path = Path with the trajectory snapshots name
            in_path = The path to the input file.
            write = Whether or not to save the APD
            verbose = Wheter or not to print calculation status

        outputs:
            apd_last = The APD for the last trajectory snapsot
        '''

        if verbose:
            print('Calculating APD from: '+traj_path)

        df, counts = traj.info(traj_path)

        # Find the box lengths
        df['dx'] = df['xhi']-df['xlo']  # Am
        df['dy'] = df['yhi']-df['ylo']  # Am
        df['dz'] = df['zhi']-df['zlo']  # Am

        df['volume'] = df['dx']*df['dy']*df['dz']  # Am^3

        # Gather input file parameters
        allelements = {}
        with open(in_path) as f:
            for line in f:
                line = line.strip().split(' ')
                if 'pair_coeff' in line:
                    line = [i for i in line if i != ''][4:]

                    count = 1
                    for item in line:
                        allelements[count] = item
                        count += 1

        # Element counts from actual element
        self.natoms = 0  # Count the total number of atoms
        elements = {}
        for key, count in counts.items():
            self.natoms += count

            element = mg.Element(allelements[key])

            try:
                atomicradii = element.metallic_radius  # in Am

            except Exception:
                atomicradii = element.atomic_radius  # in Am

            atomicvol = volume_sphere(atomicradii)  # in Am^3

            elements[allelements[key]] = {
                                          'counts': count,
                                          'radius': atomicradii,
                                          'volume': atomicvol,
                                          }

        dfelprops = pd.DataFrame(elements).T

        # Calculate the atomic packing density (APD)
        d = np.sum(dfelprops['counts']*dfelprops['volume'])
        d /= df['volume']
        df['apd'] = d

        apd_last = df['apd'].values[-1]

        print(self.datapath)
        if write:

            write_name = os.path.join(self.datapath, 'apd_last.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(apd_last))

        return apd_last

    def vp(
           self,
           threshold=0.1,
           write=True,
           plot=True,
           first=10,
           verbose=True
           ):
        '''
        Do the cluster analysis for n_5 >= 10 Voronoi polyhedra (VP).

        inputs:
            self = The object reference
            threshold = The maximum length for a VP edge
            write = Whether or not to save the fractions and temperatures
            plot = Whether or not to plot the fractions and temperatures
            verbose = Wheter or not to print calculation status

        outputs:
            max_number = The number of VP types considered for maximum variance
            max_variance = The maximum variance
            variety =  The variety of clusters
        '''

        if verbose:
            print(
                  'Calculating VP variety and ' +
                  'variance at highest temperature hold'
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

        try:
            self.file_system

        except Exception:
            message = 'Need to specify LAMMPS output file.'
            raise ValueError(message)

        # Merge dfsys and dftraj on matching steps
        df = pd.merge(
                      self.dfsys.loc[:, self.dfsys.columns != 'Volume'],
                      self.dftraj, on=['Step', 'time']
                      )

        interval = (sum(self.runsteps[:3]), sum(self.runsteps[:4]))
        condition = (df['Step'] >= interval[0]) & (df['Step'] <= interval[1])
        df = df[condition]
        df = df.reset_index(drop=True)

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

        # The average number of types of VP seen over all atoms from all frames
        variety = counts.shape[0]/total

        # Calculate variance from list including ordered fractions of VP
        variance = []
        for i in range(1, counts.shape[0]):
            variance.append(np.var(fractions[:i]))

        variance = np.array(variance)  # Numpy array

        max_index = np.argmax(variance)
        max_variance = variance[max_index]
        max_number = max_index+1

        # Create a directory for the analysis files
        if write:

            # Export all variances calculated
            np.savetxt(
                       os.path.join(self.datapath, 'variance.txt'),
                       variance,
                       )

            # Export the variety calculated
            write_name = os.path.join(self.datapath, 'variety.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(variety))

        if plot:

            # Plot variance
            fig, ax = pl.subplots()

            ax.plot(
                    np.arange(1, len(variance)+1),
                    variance,
                    marker='.',
                    linestyle='none',
                    label='Data for '+str(frames)+' frames'
                    )

            max_label = 'Maximum variance: '+str((max_number, max_variance))
            ax.axvline(
                       max_number,
                       linestyle=':',
                       color='r',
                       label=max_label
                       )

            ax.set_ylabel('Variance of VP [-]')
            ax.set_xlabel('Number of VP Types [-]')

            ax.grid()
            ax.legend()

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'variance.png'))

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

        return max_number, max_variance, variety
