#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703
# pylint: disable=C0115,C0116,C0415,E1101

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import calendar
import numpy as np
import numpy.ma as _ma
import json as _json
import os as _os
import pandas as pd

from scipy.optimize import minimize
from collections import namedtuple as _nt
from datetime import datetime as _dt

from json import JSONEncoder as _je
from os.path import join as _join
from pathlib import Path as _path
from warnings import warn as _warn

from netCDF4 import Dataset
import xarray as xr

from copy import deepcopy

from settings import setting as s

import itertools

Location = _nt('Location', 'ilat ilon lat lon')


def _expandgrid_(*itrs):
    product = list(itertools.product(*itrs))
    return [[x[i] for x in product] for i in range(len(itrs))]


def _check_slice_(s, types=(int, float, str)):
    def check_type(x):
        typ = None
        if x is not None:
            if isinstance(x, types):
                typ = type(x)
            else:
                err_str = "slice value type can only be float or int"
                raise ValueError(err_str)
        return typ

    if s is None:
        s = slice(None, None)
    if not isinstance(s, slice):
        raise ValueError("{} argument must be a slice object")

    typ = list(set((check_type(s.start), check_type(s.stop))))
    typ = list(filter(None.__ne__, typ))
    ret = None
    if len(typ) == 0:
        ret = None
    elif len(typ) == 1:
        if typ[0] in types:
            ret = typ[0]
        else:
            err_str = f"slice value types can only be one of {types}"
            raise ValueError(err_str)
    else:
        raise ValueError("slice is not in proper format")
    return (s, ret)


def _get_data2_(dom_names, proj, pol_names, years, months,
                slice_dates=None, slice_ilats=None, slice_ilons=None,
                gd=None, iterate=True):
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'

    # check arg types
    if not isinstance(dom_names, list):
        dom_names = [dom_names]

    if not isinstance(pol_names, list):
        pol_names = [pol_names]

    if not isinstance(years, list):
        years = [years]

    if not isinstance(months, list):
        months = [months]

    if isinstance(proj, str):
        proj = s.get_proj_by_name(proj)

    slice_dates, _ = _check_slice_(slice_dates, (str, int))
    slice_ilats, type_ilats = _check_slice_(slice_ilats, (int, float))
    slice_ilons, type_ilons = _check_slice_(slice_ilons, (int, float))

    if (type_ilats, type_ilons) == (float, float):
        if gd is not None:
            l1 = gd.nearest_grid_loc(slice_ilats.start, slice_ilons.start)
            l2 = gd.nearest_grid_loc(slice_ilats.stop, slice_ilons.stop)
            ilats = (l1.ilat, l2.ilat)
            ilats = (min(ilats), max(ilats))
            if ilats[0] == ilats[1]:
                ilats = (min(ilats), max(ilats) + 1)
            ilons = (l1.ilon, l2.ilon)
            ilons = (min(ilons), max(ilons))
            if ilons[0] == ilons[1]:
                ilons = (min(ilons), max(ilons) + 1)
            slice_ilats = slice(ilats[0], ilats[1])
            slice_ilons = slice(ilons[0], ilons[1])
        else:
            err_msg = 'if lat/lon slices are float, gd (GridData) must be given'
            raise ValueError(err_msg)

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

                xy = GridData(_join(DIR_MCIP, CRO_FILE))
                lats = xy.lats[slice_ilats, slice_ilons]
                lons = xy.lons[slice_ilats, slice_ilons]

                dates = ncDates.from_proj(proj, dn, y, m)
                # dates = _get_dates_(proj, dn, y, m)
                dates = dates.loc(slice_dates)
                # dates = dates.loc[slice_dates]
                # date_values = dates.coords.to_index()
                # date_indices = dates.to_series().to_list()

                nco = Dataset(_join(DIR_POST, COM_FILE))

                for pol_name in pol_names:

                    if iterate:
                        for i, k in enumerate(dates.indices):
                            pol = nco.variables[pol_name]
                            pol = pol[k, :, slice_ilats, slice_ilons]
                            pol = xr.DataArray(
                                pol, dims=['t', 'l', 'y', 'x'],
                                coords={'time': (('t'), [dates.values[i]]),
                                        'layer': (('l'), [1]),
                                        'Latitude': (('y', 'x'), lats),
                                        'Longitude': (('y', 'x'), lons)})
                            pol = pol.squeeze('l')
                            pol = pol.drop(labels='layer')
                            yield pol
                    else:
                        if len(dates.indices) > 0:
                            pol = nco.variables[pol_name]
                            pol = pol[dates.indices, :, slice_ilats,
                                      slice_ilons]
                            pol = xr.DataArray(
                                pol, dims=['t', 'l', 'y', 'x'],
                                coords={'time': (('t'), dates.values),
                                        'layer': (('l'), [1]),
                                        'Latitude': (('y', 'x'), lats),
                                        'Longitude': (('y', 'x'), lons)})
                            pol = pol.squeeze('l')
                            pol = pol.drop(labels='layer')
                            yield pol
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


def _haversine_(lon1, lat1, lon2, lat2):
    from math import radians, cos, sin, asin, sqrt
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers
    return c * r


class ncDates:
    def __init__(self, xarray, nc_file=None):
        self.xarr = xarray
        self.values = xarray.coords.to_index()
        self.indices = xarray.to_series().to_list()
        self.nc_file = nc_file

    def __repr__(self):
        s = f'ncDates Object ({len(self.values)} dates):\n'
        s += f'  From : {self.values.min()}\n'
        s += f'  To   : {self.values.max()}'
        return s

    def loc(self, from_date=None, to_date=None):
        if from_date is None and to_date is None:
            raise ValueError('from_date and to_date cannot be both None')
        if not isinstance(from_date, slice):
            if from_date is not None and not isinstance(from_date, (int, str)):
                raise ValueError('from_date arg can be str or int')
            if to_date is not None and not isinstance(to_date, (int, str)):
                raise ValueError('to_date arg can be str or int')
            slice_dates = slice(from_date, to_date)
        else:
            slice_dates = from_date
        slice_dates, type_dates = _check_slice_(slice_dates, (str, int))
        if type_dates == int:
            sliced = self.xarr[slice_dates].coords['time'].values
            slice_dates = slice(str(np.datetime_as_string(sliced.min())),
                                str(np.datetime_as_string(sliced.max())))
        return ncDates(self.xarr.loc[slice_dates])

    @classmethod
    def _get_dates_from_nc_(cls, nc_file):
        nc = Dataset(nc_file)
        TFLAG = nc.variables['TFLAG'][:, 0, :]
        nc.close()
        dates = [f'{i[0]} {i[1]:06d}' for i in TFLAG]
        # dates = ['{} {:06d}'.format(i[0], i[1]) for i in TFLAG]
        dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
        xarr = xr.DataArray(range(0, len(dates)),
                            [("time", pd.to_datetime(dates))])
        return xarr

    @classmethod
    def from_nc_file(cls, nc_file):
        return cls(cls._get_dates_from_nc_(nc_file), nc_file)

    @classmethod
    def from_proj(cls, proj, dom_name, years, months):
        if not isinstance(years, list):
            years = [years]
        if not isinstance(months, list):
            months = [months]
        dom = proj.get_dom_by_name(dom_name)
        DIR_POST = proj.path.post
        fmt = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{{}}{{:02d}}.nc'
        fmt = fmt.format(proj.cmaq_ver, proj.compiler,
                         proj.name, dom.size, dom.name)
        fmt = _join(DIR_POST, fmt)

        list_of_nc_files = [fmt.format(y, m) for y in years for m in months]
        xarr = [cls._get_dates_from_nc_(f) for f in list_of_nc_files]
        xarr = xr.concat(xarr, 'time')
        return cls(xarr, list_of_nc_files)


class GridData:
    def __init__(self, nc_cro_file, lat_name='LAT', lon_name='LON'):
        from pyproj import Proj
        attributes = ('XORIG YORIG XCENT YCENT XCELL YCELL NROWS NCOLS ' +
                      'P_ALP P_BET GDTYP P_GAM').split(' ')
        variables = 'proj lons lats xlons ylats'.split(' ')

        nco = Dataset(nc_cro_file)
        atts = {k: nco.getncattr(k) for k in attributes}
        lons = np.squeeze(nco.variables[lon_name][:])
        lats = np.squeeze(nco.variables[lat_name][:])
        nco.close()

        # pylint: disable=C0209
        proj = Proj("+proj=lcc +lon_0={XCENT} +lat_0={YCENT} \
            lat_1={P_ALP} +lat_2={P_BET}".format(**atts))
        xlons, ylats = proj(lons, lats)  # pylint: disable=E0633
        data = dict(zip(variables, [proj, lons, lats, xlons, ylats]))
        self.nc_cro_file = nc_cro_file
        self.__dict__.update(data)
        self.__dict__.update(atts)

    def nearest_grid_loc(self, lat, lon):
        from math import floor
        x, y = self.proj(lon, lat)
        ilon, ilat = (floor((x - self.XORIG) / self.XCELL),
                      floor((y - self.YORIG) / self.YCELL))
        if ilon < 0 or ilon >= self.NCOLS:
            raise ValueError('lon is out of domain bounds')
        if ilat < 0 or ilat >= self.NROWS:
            raise ValueError('lat is out of domain bounds')
        return Location(ilat, ilon,
                        self.lats[ilat, ilon],
                        self.lons[ilat, ilon])

    def nearest_grid_hav(self, lat, lon):
        # lon, lat = 27.1633, 38.4217 - İzmir
        lons = self.lons
        lats = self.lats
        h = _ma.empty(lats.shape)
        for i in range(lats.shape[0]):
            for j in range(lats.shape[1]):
                h[i, j] = _haversine_(lon, lat, lons[i, j], lats[i, j])
        g = _ma.where(h == h.min())
        return Location(g[0][0], g[1][0], lats[g][0], lons[g][0])

    @classmethod
    def from_dom(cls, proj, dom):
        fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'

        if isinstance(proj, str):
            proj = s.get_proj_by_name(proj)

        if isinstance(dom, str):
            dom = proj.get_dom_by_name(dom)

        m = proj.months[0]
        month_name = calendar.month_name[m].lower()
        DIR_MCIP = _join(proj.path.mcip, f'{dom.size}km',
                         dom.name, f'{month_name}_monthly')
        CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name,
                                  proj.years[0], m)
        return GridData(_join(DIR_MCIP, CRO_FILE))


class Domains:  # pylint: disable=R0903
    def __init__(self, dic):
        self.__dict__.update(dic)

    def __repr__(self):
        return 'Domains:\n' + '\n'.join([f" * {k}" for k in self.__dict__])


class Domain:
    def __init__(self, name, gd, pp):
        self.name = name
        self.gd = gd
        self._pp_ = deepcopy(pp)
        self._pp_._gd_ = gd
        self._pp_.dom_names = [name]

    def __repr__(self):
        s = 'Domain:\n'
        s += f'  Name: {self.name}\n'
        s += f'  Size: {self.gd.XCELL} x {self.gd.YCELL}\n'
        s += f'  ncol: {self.gd.NCOLS}\n'
        s += f'  nrow: {self.gd.NROWS}\n'
        return s

    def get_data(self, slice_dates=None, slice_ilats=None, slice_ilons=None):
        return self._pp_.get_data(slice_dates, slice_ilats, slice_ilons)

    def iterate(self, slice_dates=None, slice_ilats=None, slice_ilons=None):
        yield self._pp_._iterate_(slice_dates, slice_ilats, slice_ilons)


class PostProc:

    def __init__(self, proj_names=None, dom_names=None, pol_names=None):
        years, months = None, None
        if proj_names is None:
            proj_names = s.get_proj_names()
        if not isinstance(proj_names, list):
            proj_names = [proj_names]
        _proj_ = s.get_proj_by_name(proj_names)

        if dom_names is None:
            dom_names = [p.get_dom_names() for p in _proj_.values()]
            dom_names = list(set(sum(dom_names, [])))
        if not isinstance(dom_names, list):
            dom_names = [dom_names]

        # check domain was defined in the project
        for d in dom_names:
            for p in _proj_.values():
                if d not in p.get_dom_names():
                    err_str = f"Domain {d} is not belong to project {p.name}"
                    raise ValueError(err_str)

        # check domain sizes are consistent among projects
        p0 = next(iter(_proj_.values()))
        doms = {d: p0.get_dom_by_name(d) for d in dom_names}
        dom_size = {k: (d.nrow, d.ncol) for k, d in doms.items()}
        for p in _proj_.values():
            ds = {k: (d.nrow, d.ncol) for k, d in
                  {d: p.get_dom_by_name(d) for d in dom_names}.items()}
            for k, d in dom_size.items():
                if d != ds[k]:
                    err_str = "Incompatible domain size for {}.{} and {}.{}"
                    raise ValueError(err_str.format(p0.name, d, p.name, ds[k]))

        if pol_names is None:
            pol_names = list(p0.pols.keys())
        if not isinstance(pol_names, list):
            pol_names = [pol_names]

        # check pol names are defined in project settings
        for p in pol_names:
            if p.lower() not in p0.pols.keys():
                err_str = "Pollutant {} can not be found in settings"
                raise ValueError(err_str.format(p))

        pols = {p: p0.pols[p] for p in pol_names}

        if years is None:
            years = [p.years for p in _proj_.values()]
            years = list(set(sum(years, [])))
        if not isinstance(years, list):
            years = [years]

        if months is None:
            months = [p.months for p in _proj_.values()]
            months = list(set(sum(months, [])))
        if not isinstance(months, list):
            months = [months]

        self.proj_names = proj_names
        self.pol_names = pol_names
        self._proj_ = _proj_
        self._gd_ = None
        self.dom_names = dom_names
        self.pols = pols
        self.years = years
        self.months = months
        self.dates = ncDates.from_proj(p0, dom_names[0], years, months)
        doms2 = {k: Domain(k, GridData.from_dom(p0, doms[k]), self)
                 for k in dom_names}
        self.domains = Domains(doms2)

    def get_data(self, slice_dates=None, slice_ilats=None, slice_ilons=None):
        ret = list(self._iterate_(slice_dates, slice_ilats,
                                  slice_ilons, iterate=False))
        if len(ret) > 1:
            ret = {k: [r[k] for r in ret] for k in ret[0].keys()}
            ret = {k: xr.concat(v, 't') for k, v in ret.items()}
        else:
            ret = ret[0]
        return ret

    def iterate(self, slice_dates=None, slice_ilats=None, slice_ilons=None):
        yield self._iterate_(slice_dates, slice_ilats, slice_ilons)

    def _iterate_(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
                  iterate=True):

        slice_dates, type_dates = _check_slice_(slice_dates, (str, int))
        if type_dates == int:
            sliced = self.dates.xarr[slice_dates].coords['time'].values
            slice_dates = slice(str(np.datetime_as_string(sliced.min())),
                                str(np.datetime_as_string(sliced.max())))

        pol_names = [p.name for p in self.pols.values()]
        dim_names = ['domain', 'project', 'pol_name']
        g = _expandgrid_(self.dom_names, self.proj_names, pol_names)
        g = list(map(tuple, zip(*g)))

        it = {i: _get_data2_(i[0], i[1], i[2], self.years, self.months,
                             slice_dates, slice_ilats, slice_ilons,
                             gd=self._gd_, iterate=iterate)
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
