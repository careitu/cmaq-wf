#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create datasets for analysis
~~~~~~~~
"""

import numpy as _np
import xarray as _xr

from DAL import Location
from DAL import PostProc
from DALObs import StaData
from DALObs import load_obs_data


# Conversion coefficients
coefs = {'co': 1.18, 'nox': 1.95, 'o3': 2.03}
# domain names
dom_names=['aegean', 'mediterranean',
           'south_central_anatolia', 'central_blacksea']
sta = StaData()  # load Stations

# load observations
obs = load_obs_data(sta, ['co', 'nox', 'o3', 'so2', 'pm10'],
                    simplify=False)
obs_pol_names = obs.coords['pol_name'].values.tolist()

# load time-series data from post-proc files
pp = PostProc(dom_names=dom_names, pol_names=obs_pol_names)
xd = [d.get_data_loc(loc=list(s.loc.values()), delta=0, layer_mean=False,
                     simplify=False) for s, d in zip(sta, pp.domains)]
doms = _xr.concat(xd, dim='sta')
for k, c in coefs.items():
    doms.loc[{'pol_name': k}] = doms.loc[{'pol_name': k}] * c

# Sample
# stations located in aegean
# locs=list(sta.aegean.loc.values())
# ngl = [pp.domains.aegean.nearest_grid_loc(lo) for lo in locs]
# daegean = pp.domains.aegean.get_data_loc(
#     loc=locs, delta=0, layer_mean=False, simplify=False)

# make sure observation dates are same with model dates
obs_dates = obs.coords['time'].values
model_dates = doms.coords['time'].values
actual = obs[:,:,:,[i in model_dates for i in obs_dates],:,:]

# if NA percent is greater than 25%, exclude data from analysis
# selected = []
# for i in range(actual.shape[0]):  # loop over stations
#     sta = actual[i]
#     for j in range(sta.shape[1]):  # loop over pollutants
#         sta2 = _np.squeeze(sta[:, j])
#         na_per = float(100 * _np.isnan(sta2).sum() / len(sta2))
#         if na_per < 25:
#             # o = actual.isel(sta=i, pol_name=j, drop=False)
#             o = actual[i, :, j]
#             o = o.expand_dims(['sta', 'pol_name'], [0, 2])
#             selected.append(o)
#         else:
#             print(sta2.coords['sta'].values.tolist(), ', ',
#                   sta2.coords['pol_name'].values.tolist(), ', ',
#                   f'na_per= %{na_per:.2f}', sep = '')
# sel = _xr.concat(selected, dim=('pol_name'))

# concat observations and model results
data = _xr.concat([actual, doms], dim='project')

# if you want you can drop 1-length dimensions
data = _np.squeeze(data)

# save time series data to file
data.to_netcdf('timeseries.nc')

# load mean data from post-proc files for reduction coefficients
xd2 = [d.get_data_loc(loc=list(s.loc.values()), delta=1, layer_mean=True,
                     simplify=False) for s, d in zip(sta, pp.domains)]
doms2 = _xr.concat(xd2, dim='sta')
for k, c in coefs.items():
    doms2.loc[{'pol_name': k}] = doms2.loc[{'pol_name': k}] * c
data2 = _xr.concat([actual, doms2], dim='project')
data2 = _np.squeeze(data2)
data2.to_netcdf('reduction.nc')
