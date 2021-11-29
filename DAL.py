#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import calendar
import numpy as np
import json as _json
import os as _os
import pandas as pd

from datetime import datetime as _dt

from json import JSONEncoder as _je
from os.path import join as _join
from pathlib import Path as _path
from warnings import warn as _warn

from netCDF4 import Dataset
import xarray as xr

from settings import setting as s

import itertools


def _expandgrid_(*itrs):
    product = list(itertools.product(*itrs))
    return [[x[i] for x in product] for i in range(len(itrs))]


def _get_latlon_from_cro_(nc_cro_file, lat_name='LAT', lon_name='LON'):
    """
    Get latitude and longitude values from CMAQ CRO files
    """
    nco = Dataset(nc_cro_file)
    lon = np.squeeze(nco.variables[lon_name][:])
    lat = np.squeeze(nco.variables[lat_name][:])
    nco.close()

    return lon, lat


def _get_data_(dom_names, proj, pol_names, years, months, tstep=True):
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'

    # check arg types
    if not isinstance(dom_names, list):
        dom_names = [dom_names]
    if not isinstance(pol_names, list):
        pol_names = [pol_names]
    if not isinstance(months, list):
        months = [months]
    if isinstance(proj, str):
        proj = s.get_proj_by_name(proj)

    for dn in dom_names:
        dom = proj.get_dom_by_name(dn)
        for y in years:
            for m in months:
                month_name = calendar.month_name[m].lower()
                DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size),
                                 dom.name, month_name + '_monthly')
                DIR_POST = proj.path.post

                CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name, y, m)
                COM_FILE = fmt_com.format(proj.cmaq_ver, proj.compiler,
                                          proj.name, dom.size, dom.name, y, m)

                lon, lat = _get_latlon_from_cro_(_join(DIR_MCIP, CRO_FILE))
                nco = Dataset(_join(DIR_POST, COM_FILE))

                for pol_name in pol_names:
                    TFLAG = np.squeeze(nco.variables['TFLAG'][:])
                    tf = TFLAG.reshape(TFLAG.shape[0] * TFLAG.shape[1], 2)
                    tf = np.unique(tf, axis=0)
                    dates = ['{} {:06d}'.format(i[0], i[1]) for i in tf]
                    dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
                    dates = pd.date_range(str(min(dates)), str(max(dates)),
                                          freq='H')
                    if tstep:
                        for i in range(len(dates)):
                            date = dates[i]
                            pol = np.asarray(nco.variables[pol_name][i][:])
                            yield xr.DataArray(
                                pol, dims=['t', 'y', 'x'],
                                coords={'time': (('t'), [date]),
                                        'Latitude': (('y', 'x'), lat),
                                        'Longitude': (('y', 'x'), lon)})
                    else:
                        pol = np.squeeze(nco.variables[pol_name][:])
                        yield xr.DataArray(
                            pol, dims=['t', 'y', 'x'],
                            coords={'time': (('t'), dates),
                                    'Latitude': (('y', 'x'), lat),
                                    'Longitude': (('y', 'x'), lon)})
                nco.close()


def _concat_(x, dim_names, recursive=True):
    if len(dim_names) != len(list(x.keys())[0]):
        msg = "length of x.keys() must be equal to length of dim_names"
        raise ValueError(msg)

    g2 = list(set(i[:-1] for i in x.keys()))
    g2 = [i[0] if len(i) == 1 else i for i in g2]
    doms = {i: list() for i in g2}
    for k, v in doms.items():
        for k2, d2 in x.items():
            k2 = k2[:-1]
            if isinstance(k2, tuple):
                if len(k2) == 1:
                    k2 = k2[0]
            if k == k2:
                doms[k].append(d2)

    # combine by last dim
    for k, v in doms.items():
        doms[k] = xr.concat(doms[k], dim_names[-1])

    if recursive:
        k = list(doms.keys())[0]
        if not isinstance(k, str):
            if len(k) > 1:
                return _concat_(doms, dim_names[:-1])

    return doms


class PostProc:
    def __init__(self, proj_names=None, dom_names=None, pol_names=None,
                 years=None, months=None):
        if proj_names is None:
            proj_names = s.get_proj_names()

        if not isinstance(proj_names, list):
            proj_names = [proj_names]

        _proj_ = s.get_proj_by_name(proj_names)

        if dom_names is None:
            dom_names = []
            for p in _proj_.values():
                dom_names = dom_names + p.get_dom_names()
            dom_names = list(set(dom_names))

        if not isinstance(dom_names, list):
            dom_names = [dom_names]

        for d in dom_names:
            for p in _proj_.values():
                if d not in p.get_dom_names():
                    str = "Domain {} is not belong to project {}"
                    raise ValueError(str.format(d, p.name))

        if years is None:
            years = []
            for p in _proj_.values():
                years = years + p.years
            years = list(set(years))

        if not isinstance(years, list):
            years = [years]

        if months is None:
            months = []
            for p in _proj_.values():
                months = months + p.months
            months = list(set(months))

        if not isinstance(months, list):
            months = [months]

        self.proj_names = proj_names
        self._proj_ = _proj_
        self.dom_names = dom_names
        self.pol_names = pol_names
        self.years = years
        self.months = months

    def get_data(self):
        # years = [2015]
        # months = [1, 2, 3]
        # proj_names = ['cityair', 'cityair_future', 'cityair_future_wama']
        # dom_names = ['aegean', 'central_blacksea', 'mediterranean']
        # pol_names = ['NOX', 'O3', 'CO', 'SO2_UGM3', 'PM10', 'PM25_TOT']

        dim_names = ['domain', 'project', 'pol_name']
        g = _expandgrid_(self.dom_names, self.proj_names, self.pol_names)
        g = list(map(tuple, zip(*g)))

        it = {i: _get_data_(i[0], i[1], i[2], self.years, self.months)
              for i in g}
        counter = 0
        while True:
            try:
                d = {k: next(v) for k, v in it.items()}
                for k, v in d.items():
                    d[k] = v.assign_coords(dict(zip(dim_names, k)))
                    d[k] = d[k].expand_dims(dict(zip(dim_names, [1, 1, 1])))

                d = _concat_(d, dim_names)
                yield d

            except StopIteration:
                break
            # for d2 in d.values():
            #     t = d2.coords['time'].values[0]
            #     break
            counter = counter + 1
            # print("{} - {}".format(counter, t))

