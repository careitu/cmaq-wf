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
        self.df = pd.concat([self.rmse, self.mae['MAE'], self.max_err['MAXERR']], 1, join='inner')

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


# load time-series data
ds = xr.open_dataarray('timeseries.nc').load()
pol_names = ds.coords['pol_name'].values.tolist()

actual = ds[:, 0]  # Observations

# if NA percent is greater than 25%, exclude data from analysis
pols = {i: [] for i in pol_names}
for i in range(actual.shape[1]):  # loop over pollutants
    p = actual[:, i]
    pn = p.coords['pol_name'].values.tolist()
    for j, s in enumerate(p):
        na_per = float(100 * np.isnan(s).sum() / len(s))
        if na_per < 25:
            o = ds.isel(sta=j, pol_name=i, drop=False)
            o = o.expand_dims(['sta'], [0])
            pols[pol_names[i]].append(o)
pols = {k: xr.concat(v, dim='sta') for k, v in pols.items()}

# plot time series
for k, v in pols.items():
    d = v[:, 0:2]
    n_sta = d.shape[0]
    cw = (n_sta, n_sta) if n_sta <= 5 else max_mult(n_sta)
    g = d.plot.line(row='sta', hue='project', col_wrap = min(cw),
                    aspect=2, size=3, sharey=False)
    g.fig.suptitle(k.upper(), y = 1.0)
    file_name = '_'.join(('plot', k)) + '.pdf'
    g.fig.savefig(file_name, bbox_inches='tight')
    plt.close(g.fig)
    print(file_name)

# Calculate Metrics as pandas dataframe
m = pd.concat([Metrics(v[:, 0], v[:, 1], 1, True).df for k, v in pols.items()])


