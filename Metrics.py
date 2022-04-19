#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101


"""
Calculate performance for CMAQ vs. Observations
~~~~~~~~
"""
import numpy as np
import xarray as xr
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap as lscmap

from _helper_functions_ import max_mult

class Metrics:
    def __init__(self, actual, predicted, axis=None, skipna=None):
        self.rmse = Metrics.rmse(actual, predicted, axis, skipna).to_dataframe('RMSE')
        self.mae = Metrics.mae(actual, predicted, axis, skipna).to_dataframe('MAE')
        self.max_err = Metrics.max_error(actual, predicted, axis, skipna).to_dataframe('MAXERR')
        self.mbe = Metrics.mbe(actual, predicted, axis, skipna).to_dataframe('MBE')
        self.df = pd.concat([self.rmse, self.mae['MAE'], self.max_err['MAXERR'], self.mbe['MBE']], 1, join='inner')


    @classmethod
    def rmse(cls, actual, predicted, axis=None, skipna=None):
        dif_sq = (actual - predicted) ** 2
        return np.sqrt(dif_sq.mean(axis=axis, skipna=skipna))

    @classmethod
    def mae(cls, actual, predicted, axis=None, skipna=None):
        difa = abs(actual - predicted)
        return difa.mean(axis=axis, skipna=skipna)

    @classmethod
    def max_error(cls, actual, predicted, axis=None, skipna=None):
        difa = abs(actual - predicted)
        return difa.max(axis=axis, skipna=skipna)

    @classmethod
    def mbe(cls, actual, predicted, axis=None, skipna=None):
        difa = actual - predicted
        return difa.mean(axis=axis, skipna=skipna)

    @classmethod
    def r_2(cls, actual, predicted, axis=None, skipna=None):
        difa = actual - predicted
        return difa.mean(axis=axis, skipna=skipna)

    @classmethod
    def cov(cls, actual, predicted, axis=None):

        cv = []
        for i, j in zip(actual, predicted):
            nas = np.logical_or(np.isnan(i), np.isnan(j))
            cv.append(np.cov(i[~nas], j[~nas]))
        [np.cov(i, j) for i, j in zip(actual, predicted)]
        return np.cov(actual, predicted, axis=axis, skipna=skipna)

    @classmethod
    def pearson(actual, predicted, axis=None):
        dim_str = None if axis is None else actual.dims[axis]
        return xr.corr(actual, predicted, dim=dim_str)


def split_xarr_by_pols(data, obs, pol_names):
    # if NA percent is greater than 25%, exclude data from analysis
    pols = {i: [] for i in pol_names}
    for i in range(obs.shape[1]):  # loop over pollutants
        p = obs[:, i]
        pn = p.coords['pol_name'].values.tolist()
        for j, s in enumerate(p):
            na_per = float(100 * np.isnan(s).sum() / len(s))
            if na_per < 25:
                o = data.isel(sta=j, pol_name=i, drop=False)
                o = o.expand_dims(['sta'], [0])
                pols[pol_names[i]].append(o)
    return {k: xr.concat(v, dim='sta') for k, v in pols.items()}



# load time-series data
ds = xr.open_dataarray('timeseries.nc').load()
pol_names = ds.coords['pol_name'].values.tolist()

actual = ds[:, 0]  # Observations

# if NA percent is greater than 25%, exclude data from analysis
pols = split_xarr_by_pols(ds, actual, pol_names)

actual = pols['co'][:, 0]
predicted = pols['co'][:, 1]

# Calculate Metrics as pandas dataframe
met = {k: Metrics(v[:, 0], v[:, 1], 1, True).df
                     for k, v in pols.items()}
metrics = pd.concat([Metrics(v[:, 0], v[:, 1], 1, True).df
                     for k, v in pols.items()])
metrics = metrics.rename(columns={'pol_name': 'Pollutant'})


