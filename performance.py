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

# plot scatter
for k, v in pols.items():
    vds = v.to_dataset(dim='project')
    g = vds.plot.scatter(row='sta', x='obs',y='cityair', marker='+', c='grey', col_wrap=3, sharey=False, sharex=False)
    # plt.xlim([0,300])
    # plt.ylim([0,300])
    # plt.axline((1, 1), slope=1)
    g.fig.suptitle(k.upper(), y = 1.0)
    file_name = '_'.join(('plot_scatter', k)) + '.pdf'
    g.fig.savefig(file_name, bbox_inches='tight')
    plt.close(g.fig)
    print(file_name)



# Calculate Metrics as pandas dataframe
met = {k: Metrics(v[:, 0], v[:, 1], 1, True).df
                     for k, v in pols.items()}
metrics = pd.concat([Metrics(v[:, 0], v[:, 1], 1, True).df
                     for k, v in pols.items()])
metrics = metrics.rename(columns={'pol_name': 'Pollutant'})

# -----------------------------------------------
# RELATIVE REDUCTION FACTOR
ds2 = xr.open_dataarray('reduction.nc').load()
pol_names = ds2.coords['pol_name'].values.tolist()

obs_2015 = ds2[:, 0]
r2015 = ds2[:, 1]
r2025 = ds2[:, 2]
rwama = ds2[:, 3]
rwamb = ds2[:, 4]

rrf1 = abs(r2025 - r2015)/r2015
rrf2 = abs(rwama - r2015)/r2015
rrf3 = abs(rwamb - r2015)/r2015

obs_2025 = obs_2015 * rrf1
obs_wama = obs_2015 * rrf2
obs_wamb = obs_2015 * rrf3

obs_2025 = xr.concat([r2015, obs_2025], dim = 'project')
obs_2025.coords['project'] = ['cityair', 'cityair_future']

obs_wama = xr.concat([r2015, obs_wama], dim = 'project')
obs_wama.coords['project'] = ['cityair', 'cityair_wama']

obs_wamb = xr.concat([r2015, obs_wamb], dim = 'project')
obs_wamb.coords['project'] = ['cityair', 'cityair_wamb']

p2025 = {i: [] for i in pol_names}
for i in range(obs_2025.shape[2]):  # loop over pollutants
    p = obs_2025[:, :, i]
    pn = p.coords['pol_name'].values.tolist()
    for j in range(p.shape[1]):  # loop over stations
        s = p[:, j]
        na_per = float(100 * np.isnan(s[1]).sum() / len(s[1]))
        if na_per < 25:
            o = obs_2025.isel(sta=j, pol_name=i, drop=False)
            o = o.expand_dims(['sta'], [0])
            p2025[pol_names[i]].append(o)
p2025 = {k: xr.concat(v, dim='sta') for k, v in p2025.items()}

# plot time series
for k, v in p2025.items():
    n_sta = v.shape[0]
    cw = (n_sta, n_sta) if n_sta <= 5 else max_mult(n_sta)
    g = v.plot.line(row='sta', hue='project', col_wrap = min(cw),
                    aspect=2, size=3, sharey=False)
    g.fig.suptitle(k.upper(), y = 1.0)
    file_name = '_'.join(('plot_2025', k)) + '.pdf'
    g.fig.savefig(file_name, bbox_inches='tight')
    plt.close(g.fig)
    print(file_name)

