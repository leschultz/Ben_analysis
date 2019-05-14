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

    def __init__(self, path, data_path, plot_path):
        '''
        Create all the paths needed to save analysis data.
        '''

        # The location of the job
        self.path = path

        # Save paths
        self.datapath = os.path.join(self.path, data_path)
        self.plotpath = os.path.join(self.path, plot_path)

        # Create a directory for the analysis files
        if not os.path.exists(self.datapath):
            os.makedirs(self.datapath)

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

    def volume(self):
        '''
        Calculate the volume from box dimensions.

        inputs:
            self = The object reference

        outputs:
            dfvol = The volumes for trajectory snapshots
        '''

        print('Calculating volume')

        try:
            self.file_trajs

        except Exception:
            message = 'Need to specify trajectory file.'
            raise ValueError(message)

        df = pd.DataFrame()

        # Find the box lengths
        df['dx'] = self.dftraj['xhi']-self.dftraj['xlo']  # Am
        df['dy'] = self.dftraj['yhi']-self.dftraj['ylo']  # Am
        df['dz'] = self.dftraj['zhi']-self.dftraj['zlo']  # Am

        dfvol = pd.DataFrame()

        # Find the box volume
        dfvol['time'] = self.dftraj['time']
        dfvol['Volume'] = df['dx']*df['dy']*df['dz']  # Am^3

        return dfvol

    def find_tl(
                self,
                edges=5,
                faces=10,
                lim=0.1,
                threshold=0.1,
                write=True,
                plot=True
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

        outputs:
            tl = The liquidus temperature
        '''

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
        for frame in df['frame']:
            out = node.compute(frame)
            indexes = out.particle_properties['Voronoi Index'].array

            indexes = indexes[:, edges]  # Gather edge bin
            count = sum(indexes >= faces)  # Count condition
            fraction = count/self.natoms  # Calculate fraction

            fractions.append(fraction)

        df['fractions'] = fractions
        df = df.sort_values(by=['Temp'])

        x = df['Temp'].values
        y = df['fractions'].values

        k = 5
        s = 1
        spl = UnivariateSpline(x=x, y=y, k=k, s=s)
        xfit = np.linspace(np.min(x), np.max(x), 100)
        yfit = spl(xfit)

        if plot:

            fig, ax = pl.subplots()

            ax.plot(x, y, marker='.', linestyle='none', label='data')
            ax.plot(
                    xfit,
                    yfit,
                    color='g',
                    label='Univariate Spline (k='+str(k)+', s='+str(s)+')'
                    )

            xlabel = 'Temperature [K]'
            ylabel = r'Fractions of $n_{'+str(edges+1)+'} >= '+str(faces)+'$'
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax.grid()
            ax.legend()

            fig.tight_layout()
            pl.show()

    def etg(self, max_temp=1000, write=True, plot=True, verbose=True):
        '''
        Calculate the glass transition temperature based on E-3kt.

        inputs:
            self = The object reference
            max_temp = The maximum temperature for analysis
            write = Whether or not to save Tg
            plot = Whether or not to save plot of data

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

        # Spline fit of cut region
        k, s = (5, 1)
        condition = x <= max_temp
        xcut = x[condition]
        ycut = y[condition]

        spl = UnivariateSpline(x=xcut, y=ycut, k=k, s=s)
        xfitcut = np.linspace(xcut[0], xcut[-1], 100)
        yfitcut = spl(xfitcut)

        tg, left, right, ldata, rdata, middle_rmse = opt(xfitcut, yfitcut)

        if write:

            # Export the glass transition temperature
            write_name = os.path.join(self.datapath, 'tg_e.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(tg))

            # Export the upper temperature cutoff
            write_name = os.path.join(self.datapath, 'tg_e_t_cutoff.txt')
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
            fig.savefig(os.path.join(self.plotpath, 'etg'))
            pl.close('all')

            fig, ax = pl.subplots()

            ax.plot(ldata[:, 0], ldata[:, 1], label='left fits')
            ax.plot(rdata[:, 0], rdata[:, 1], label='right fits')
            ax.plot(ldata[:, 0], middle_rmse, label='left and right fits')

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
            fig.savefig(os.path.join(self.plotpath, 'etg_rmse'))
            pl.close('all')

        return tg

    def apd_single(self, traj_path, in_path, write=True):
        '''
        Calculate the APD from the last trajectories.

        inputs:
            traj_path = Path with the trajectory snapshots name
            in_path = The path to the input file.
            write = Whether or not to save the APD

        outputs:
            apd_last = The APD for the last trajectory snapsot
        '''

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

        if write:

            # Export the glass transition temperature
            write_name = os.path.join(self.datapath, 'atp_single.txt')
            with open(write_name, 'w+') as outfile:
                outfile.write(str(apd_last))

        return apd_last

    def vp(
           self,
           threshold=0.1,
           write=True,
           plot=True
           ):
        '''
        Do the cluster analysis for n_5 >= 10 Voronoi polyhedra (VP).

        inputs:
            self = The object reference
            threshold = The maximum length for a VP edge
            write = Whether or not to save the fractions and temperatures
            saveplot = Whether or not to plot the fractions and temperatures

        outputs:
            dfv = A dataframe containing VP variance and variety
        '''

        print(
              'Calculating VP variety and variance at highest temperature hold'
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

        variety = []
        variance = []
        for frame in df['frame']:
            out = node.compute(frame)
            indexes = out.particle_properties['Voronoi Index'].array
            coords, counts = np.unique(indexes, axis=0, return_counts=True)

            dfvp = pd.DataFrame()
            dfvp['VP'] = [tuple(i) for i in coords]
            dfvp['Step'] = int(df['Step'][df['frame'] == frame])
            dfvp['frame'] = int(df['frame'][df['frame'] == frame])
            dfvp['time'] = float(df['time'][df['frame'] == frame])
            dfvp['counts'] = counts
            dfvp['fracs'] = counts/len(indexes)
            dfvp['time'] = float(df['time'][df['frame'] == frame])

            dfvp = dfvp.sort_values(by=['fracs'], ascending=False)

            variances = []
            for i in range(1, len(counts)+1):
                top = dfvp['fracs'].values
                top = top[:i]

                variances.append(np.var(top))

            variances = np.array(variances)
            variety.append(len(dfvp)/len(indexes))
            variance.append(variances.max())

        dfv = pd.DataFrame()
        dfv['Step'] = df['Step']
        dfv['time'] = df['time']
        dfv['variety'] = variety
        dfv['max variance'] = variance

        # Create a directory for the analysis files
        if write:
            dfv.to_csv(
                       os.path.join(self.datapath, 'top_vp_fractions.txt'),
                       index=False
                       )

        if plot:

            # Plot variance
            fig, ax = pl.subplots()

            ax.plot(
                    dfv['time'],
                    dfv['max variance'],
                    marker='.',
                    linestyle='none'
                    )

            ax.set_ylabel('Maximum Variance of VP [-]')
            ax.set_xlabel('Time [ps]')
            ax.grid()
            ax.legend()

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'fracs_variance'))

            pl.close('all')

            # Plot variety
            fig, ax = pl.subplots()

            ax.plot(
                    dfv['time'],
                    dfv['variety'],
                    marker='.',
                    linestyle='none'
                    )

            ax.set_ylabel('Variety of VP [-]')
            ax.set_xlabel('Time [ps]')
            ax.grid()
            ax.legend()

            fig.tight_layout()
            fig.savefig(os.path.join(self.plotpath, 'fracs_variety'))

            pl.close('all')

        return dfv
