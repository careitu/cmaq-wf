#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101


"""
Calculate performance for CMAQ vs. Observations
~~~~~~~~
"""
import numpy as np
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap as lscmap

from _helper_functions_ import max_mult

from settings import setting as _set

from DAL import *

import sklearn.metrics as mp
from math import sqrt

def rmse(actual, predicted, skipna=None):
    l = actual.shape[2]
    s = 100 * np.isnan(actual).sum(axis=2) / l
    n = int(.values)
    return np.sqrt(((actual - predicted) ** 2).mean(axis=2, skipna=skipna))

# load time-series data
ds = xr.open_dataarray('timeseries.nc').load()
pol_names = ds.coords['pol_name'].values.tolist()

actual = ds[:, 0]  # Observations

# if NA percent is greater than 25%, exclude data from analysis
selected = []
for i in range(ds.shape[0]):
    sta = ds[i]
    for j in range(ds.shape[2]):
        sta2 = sta[:, j]
        if float(100 * np.isnan(sta2[0]).sum() / len(sta2[0])) < 25:
            selected.append(ds.isel(sta=i, pol_name=j, drop=True))
selected = xr.concat(selected, dim='sta')

select = [float(100 * np.isnan(p).sum() / len(p)) < 25 for s in ds[:, 0] for p in s]

per = 100 * np.isnan(actual).sum(axis=2) / actual.shape[2]
actual.where(per > 25)
n = int(.values)

# plot time series
for p in pol_names:
    d = ds.loc[{'pol_name': p}][:, 0:2]
    new_ds = []
    for i in range(d.shape[0]):
        if not np.isnan(d[i, 0]).all().values.tolist():
            new_ds.append(d[i])
    d2 = xr.concat(new_ds, dim='sta')
    n_sta = d2.shape[0]
    cw = (n_sta, n_sta) if n_sta <= 5 else max_mult(n_sta)
    g = d2.plot.line(row='sta', hue='project', col_wrap = min(cw),
                     aspect=2, size=3, sharey=False)
    g.fig.suptitle(p.upper(), y = 1.0)
    # plt.tight_layout()
    file_name = '_'.join(('plot', p)) + '.pdf'
    g.fig.savefig(file_name, bbox_inches='tight')
    plt.close(g.fig)
    print(file_name)
    # plt.show()

# calculate metrics
actual = ds[:, 0]
predicted = ds[:, 1]
rm = rmse(actual, predicted, False)
dif = obs - model



istasyon =
alper = data.loc[dict(project=[‘obs’,‘cityair’], sta=‘alsancak’,pol_name=‘co’)]
dod =alper.values
dod2 = dod.transpose()
dod2 =dod2[~np.isnan(dod2).any(axis=1)]
y_actual = dod2[:,0]
y_predicted = dod2[:,1]*1.18
rms = sqrt(mp.mean_squared_error(y_actual, y_predicted))
max_error  = mp.max_error(y_actual,y_predicted)
r_2  = mp.r2_score(y_actual,y_predicted)
MAE = mp.mean_absolute_error(y_actual,y_predicted)
diff = (y_actual-y_predicted)
MB = diff.mean()
CO: 1.18*model
NOX: 1.95*model
O3: 2.03*model