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
from DAL import Location as _Location

class StaData:
    def __init__(self, data=None, file='sta.csv'):
        if data is None:
            data = _pd.read_csv(file, skipinitialspace=True,
                                skip_blank_lines=True)
        if not isinstance(data, _pd.DataFrame):
            raise ValueError('data must be pandas DataFrame')

        data = {x: y for x, y in data.groupby('region', as_index=False)}
        if len(data) > 1:
            data = {k: StaData(v) for k, v in data.items()}
            self.__mem__ = data
            self.__dict__.update(data)
        else:
            data = list(data.values())[0]
            self.loc = {i[3]: _Location(float(i[4]), float(i[5]), name=i[3],
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
                        obs = _xr.DataArray(q[5], dims=['time'], 
                                            coords={'time': (('time'), q[4])})
                        obs = obs.expand_dims(['y', 'x'],
                                              [len(obs.dims), len(obs.dims) + 1])
                        obs = obs.expand_dims(
                            dict(zip(['sta', 'project', 'pol_name'],
                                     [1, 1, 1])))
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
                    if dum is not None:
                        dum1 = _xr.concat([dum1, dum], dim='sta')
            dum1 = dum1.assign_coords(domain=k)
            yield dum1

    db = _database('air', return_type='long_list')
    obs = list(load(sta, pols))
    obs = _xr.concat(obs, dim='sta')
    if simplify:
        obs = _np.squeeze(obs)
    return obs

