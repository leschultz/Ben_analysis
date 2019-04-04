from scipy.interpolate import LSQUnivariateSpline as spline
from matplotlib import pyplot as pl

from scipy.signal import argrelextrema

import numpy as np

import math
import os


def knees(x, y):
    '''
    Find the knee based on how a second derivative behaves.
    Data has to be strictly increasing.

    inputs:
        x = The x-axis data
        y = The y-axis data

    outputs:
        yspline = The y-values for fitted data
        ddyspline = The second derivative for fitted data
        splineindex = The index where the knee should occur
    '''

    # Setup the number of knots for the spline fit
    t = [np.mean(i) for i in np.array_split(x, 4)]

    # Fit a spline and find the derivatives
    s = spline(
               x=x,
               y=y,
               k=5,
               t=t,
               )

    dds = s.derivative(2)  # Take a second derivative

    n = 1000
    xnew = np.linspace(x[0], x[-1], n)  # Make smooth
    yspline = s(xnew)
    ddyspline = dds(xnew)

    # The interval for analysis has to be cut because of spline end behavior
    cut = math.floor(n*0.5)

    # Find the local maxima
    localmaxes = argrelextrema(ddyspline[:cut], np.greater)

    # Find the largest of the local maxima
    if len(localmaxes[0]) > 0:
        splineindex = localmaxes[0][np.argmax(ddyspline[localmaxes])]

    # Find the global maxima
    else:
        splineindex = np.min(np.where(ddyspline == np.max(ddyspline))[0])

    return xnew, yspline, ddyspline, splineindex


def plotknee(xdata, ydata, xfit, yfit, ddyfit, kneeindex, path, name):
    '''
    Plot data to visualize where the knee of a curve occurs.

    inputs:
        xdata = The x-axis data
        ydata = The y-axis data
        xfit = The x-axis data for the spline fit
        yfit = The y-axis data for the spline fit
        ddyfit = The y-axis data for the second derivative of yfit
        kneeindex = The index of the data of the knee
        path = The path to save the plot
        name = The beginning save name for the plots

    outputs:
        saves plots in directory where it is run
    '''

    xlabel = 'Temperature in [K]'

    if name == 'etg':
        ylabel = 'E-3kT [eV/atom]'

    if name == 'vtg':
        ylabel = 'Specific Volume [A^3/atom]'

    xlabel = 'Temperature in [K]'
    ddylabel = 'Second Derivative'

    legenddigits = 6  # The display digit length for plots legends

    fig, ax = pl.subplots(2, 1)

    ax[0].set_title('Analysis of the rate of change of the slope')

    splineintersection = (
                          float(str(xfit[kneeindex])[:legenddigits]),
                          float(str(yfit[kneeindex])[:legenddigits])
                          )

    ax[0].plot(xdata, ydata, marker='.', label='Data')
    ax[0].plot(xfit, yfit, label='Spline Fit')
    ax[0].axvline(
                  x=xfit[kneeindex],
                  color='k',
                  label='Knee at '+str(splineintersection)
                  )
    ax[0].axhline(y=yfit[kneeindex], color='k')
    ax[0].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel)
    ax[0].legend()
    ax[0].grid()

    splineintersectionerr = (
                             float(str(xfit[kneeindex])[:legenddigits]),
                             float(str(ddyfit[kneeindex])[:legenddigits])
                             )

    # Odd end behavior with spline so cut ends from display
    cutstart = 1
    cutstop = -10
    ax[1].plot(xfit[cutstart:cutstop], ddyfit[cutstart:cutstop])
    ax[1].axvline(
                  x=xfit[kneeindex],
                  color='k',
                  label='Knee at '+str(splineintersectionerr)
                  )
    ax[1].axhline(y=ddyfit[kneeindex], color='k')
    ax[1].set_ylabel(ddylabel)
    ax[1].set_xlabel(xlabel)
    ax[1].legend()
    ax[1].grid()

    savename = os.path.join(path, name)
    fig.tight_layout()
    fig.savefig(savename)

    pl.close('all')
