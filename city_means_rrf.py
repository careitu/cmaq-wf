#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101,W0621,C0103


"""
Calculate performance for CMAQ vs. Observations
~~~~~~~~
"""
import pickle
import numpy as np
import xarray as xr
# import pandas as pd

# import matplotlib
import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# from matplotlib.colors import LinearSegmentedColormap as lscmap

from _helper_functions_ import max_mult


def split_xarr_by_pols(data, obs, pol_names):
    # if NA percent is greater than 25%, exclude data from analysis
    pols = {i: [] for i in pol_names}
    for i in range(obs.shape[1]):  # loop over pollutants
        p = obs[:, i]
        # pn = p.coords['pol_name'].values.tolist()
        for j, s in enumerate(p):
            na_per = float(100 * np.isnan(s).sum() / len(s))
            if na_per < 25:
                o = data.isel(sta=j, pol_name=i, drop=False)
                o = o.expand_dims(['sta'], [0])
                pols[pol_names[i]].append(o)
    return {k: xr.concat(v, dim='sta') for k, v in pols.items()}


# RELATIVE REDUCTION FACTOR DATA
ds = xr.open_dataarray('reduction.nc').load()
pol_names = ds.coords['pol_name'].values.tolist()

actual = ds[:, 0]  # Observations

# if NA percent is greater than 25%, exclude data from analysis
pols = split_xarr_by_pols(ds, actual, pol_names)

pols2 = {k: None for k in pol_names}
for k, v in pols.items():
    sta_names = v.coords['sta']
    sta_names = sta_names.to_index().tolist()
    cities = [s.split('-')[0].strip() for s in sta_names]
    cities = [c.split(' ')[0].strip() for c in cities]
    cities2 = list(set(cities))
    print(cities2)
    cit = []
    for c in cities2:
        city = v[[c == c2 for c2 in cities]].mean(axis=0, skipna=True)
        city = city.expand_dims(dict(zip(['city'], [1])))
        city = city.assign_coords(city=[c])
        cit.append(city)
    cit2 = xr.concat(cit, dim='city')
    pols2[k] = cit2

with open('city_means_rrf_dict.pickle', 'wb') as handle:
    pickle.dump(pols2, handle, protocol=pickle.HIGHEST_PROTOCOL)


obs_2015 = ds[:, 0]
r2015 = ds[:, 1]
r2025 = ds[:, 2]
rwama = ds[:, 3]
rwamb = ds[:, 4]

rrf1 = abs(r2025 - r2015)/r2015
rrf2 = abs(rwama - r2015)/r2015
rrf3 = abs(rwamb - r2015)/r2015

obs_2025 = obs_2015 * rrf1
obs_wama = obs_2015 * rrf2
obs_wamb = obs_2015 * rrf3

obs_2025 = xr.concat([r2015, obs_2025], dim='project')
obs_2025.coords['project'] = ['cityair', 'cityair_future']

obs_wama = xr.concat([r2015, obs_wama], dim='project')
obs_wama.coords['project'] = ['cityair', 'cityair_wama']

obs_wamb = xr.concat([r2015, obs_wamb], dim='project')
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
    g = v.plot.line(row='sta', hue='project', col_wrap=min(cw),
                    aspect=2, size=3, sharey=False)
    g.fig.suptitle(k.upper(), y=1.0)
    FILENAME = '_'.join(('plot_2025', k)) + '.pdf'
    g.fig.savefig(FILENAME, bbox_inches='tight')
    plt.close(g.fig)
    print(FILENAME)
