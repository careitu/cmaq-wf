#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Observations Module
~~~~~~~~
"""
import numpy as _np
import pandas as _pd
import xarray as _xr
from airdb import Database as _database
# from DAL import Location
# from DAL import PostProc

class StaData:
    def __init__(self, data=None, file='sta.csv'):
        if data is None:
            data = _pd.read_csv(file, skipinitialspace=True)
        if not isinstance(data, _pd.DataFrame):
            raise ValueError('data must be pandas DataFrame')

        data = {x: y for x, y in data.groupby('region', as_index=False)}
        if len(data) > 1:
            data = {k: StaData(v) for k, v in data.items()}
            self.__mem__ = data
            self.__dict__.update(data)
        else:
            data = list(data.values())[0]
            self.loc = {i[3]: Location(float(i[5]), float(i[4]), name=i[3],
                                       city=i[1], region=i[0])
                        for i in  data.values.tolist()}
            self.__mem__ = self.loc
            self.df = data

    def __repr__(self):
        return 'Stations:\n' + '\n'.join([f" * {k}" for k in self.__mem__])


    def __iter__(self):
        for i in self.__mem__.values():
            yield i

def load_obs_data(
    sta, pols=['co', 'nox', 'o3', 'so2', 'pm10', 'pm25'],
    date=['>2015-01-01', '<2015-04-01'], simplify=True):

    def load(sta, pols):
        dom_names = list(sta.__mem__.keys())
        for k, reg in sta.__mem__.items():
            dum1 = None
            for i, r2 in reg.df.iterrows():
                dum = None
                for pn2 in pols:
                    q2 = db.query(param=pn2, city=r2[1], sta=r2[2], 
                                  date=date)
                    if len(q2) > 0:
                        pn, r, q = pn2, r2, q2
                        obs = _xr.DataArray(q[5], dims=['t'], 
                                            coords={'time': (('t'), q[4])})
                        obs = obs.expand_dims(['y', 'x'],
                                              [len(obs.dims), len(obs.dims) + 1])
                        obs = obs.expand_dims(
                            dict(zip(['sta', 'project', 'pol_name'],
                                     [1, 1, 1, 1])))
                        obs = obs.assign_coords(
                            Latitude=('y', [r[4]]),
                            Longitude=('x', [r[5]]),
                            project=['obs'],
                            pol_name=[pn], sta=[r[3]])
                        if dum is None:
                            dum = obs
                        else:
                            dum = _xr.concat([dum, obs], dim='pol_name')
                if dum1 is None:
                    dum1 = dum
                else:
                    dum1 = _xr.concat([dum1, dum], dim='sta')
            yield dum1

    db = _database('air', return_type='long_list')
    obs = list(load(sta, pols))
    obs = _xr.concat(obs, dim='sta')
    if simplify:
        obs = _np.squeeze(obs)
    return obs

# --------------

# load Stations
sta = StaData()
dom_names=['aegean', 'mediterranean',
           'south_central_anatolia', 'central_blacksea']
pols = dict(zip(['co', 'nox', 'o3', 'so2', 'pm10', 'pm25'], 
                ['co', 'nox', 'o3', 'so2_ugm3', 'pm10', 'pm25_tot']))

# load observations
obs = load_obs_data(sta, simplify=False)
obs_dates = obs.coords['time'].values
obs_pol_names = obs.coords['pol_name'].values.tolist()
# obs.coords['pol_name'] = ['co', 'nox', 'o3', 'pm10', 'so2_ugm3']
# these are the pollutant names which are in station data.
selected_pols = [pols[p] for p in obs_pol_names]

# load model data
pp = PostProc(dom_names=dom_names, pol_names=obs_pol_names)

# stations located in aegean
locs=list(sta.aegean.loc.values())
# ngl = [pp.domains.aegean.nearest_grid_loc(lo) for lo in locs]
daegean = pp.domains.aegean.get_data_loc(
    loc=locs, delta=0, layer_mean=False, simplify=False)

# stations located in mediterranean
locs=list(sta.mediterranean.loc.values())
# ngl = [pp.domains.aegean.nearest_grid_loc(lo) for lo in locs]
dmed = pp.domains.mediterranean.get_data_loc(
    loc=locs, delta=0, layer_mean=False, simplify=False)

doms = _xr.concat([daegean, dmed], dim='sta')
# doms = _np.squeeze(doms)

model_dates = daegean.coords['time'].values

# make sure observation dates are same with model dates
sel_obs = obs[:,:,:,[i in model_dates for i in obs_dates],:,:]

data = _xr.concat([sel_obs, doms], dim='project')