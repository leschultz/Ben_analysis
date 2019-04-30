#!/usr/bin/python3

from scipy.interpolate import UnivariateSpline as interpolate
from matplotlib import pyplot as pl

from scipy.stats import linregress
from os.path import join

import numpy as np
import os


def rmse(act, pred):
    act = np.array(act)
    pred = np.array(pred)

    e = (pred-act)**2
    e = e.mean()
    e = np.sqrt(e)

    return e


def opt(x, y, tol=0.1):
    '''
    '''

    x = np.array(x)
    y = np.array(y)

    left_rmse = []
    right_rmse = []

    n = len(x)
    indexes = list(range(n-1))
    ldata = []
    rdata = []
    for i in indexes:
        xl = x[:n-i]
        xr = x[i:]

        yl = y[:n-i]
        yr = y[i:]

        ml, il, _, _, _ = linregress(xl, yl)
        mr, ir, _, _, _ = linregress(xr, yr)

        yfitl = ml*xl+il
        yfitr = mr*xr+ir

        rmsel = rmse(yl, yfitl)
        rmser = rmse(yr, yfitr)

        if i > 0:
            ldata.append((xl[-1], rmsel))
            rdata.append((xr[0], rmser))

        left_rmse.append(rmsel)
        right_rmse.append(rmser)

    left_rmse = np.array(left_rmse)
    right_rmse = np.array(right_rmse)

    ldata = np.array(ldata)
    rdata = np.array(rdata[::-1])

    middle_rmse = (ldata[:, 1]+rdata[:, 1])/2
    mcut = np.argmin(middle_rmse)

    lcut = np.argmax(left_rmse <= tol*left_rmse.max())
    rcut = np.argmax(right_rmse <= tol*right_rmse.max())
    xcut = ldata[mcut, 0]

    left = x[n-lcut]
    right = x[rcut]

    return xcut, left, right, ldata, rdata, middle_rmse
