#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

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

from settings import setting as s

import itertools

Location = _nt('Location', 'ilat ilon lat lon')


def _expandgrid_(*itrs):
    product = list(itertools.product(*itrs))
    return [[x[i] for x in product] for i in range(len(itrs))]


def _get_latlon_from_cro_(nc_cro_file, lat_name='LAT', lon_name='LON'):
    """
    Get XY coordinate values from CMAQ CRO files
    """
    import pyproj
    attr_str = 'XORIG YORIG XCENT YCENT XCELL YCELL NROWS NCOLS P_ALP P_BET'
    attr = _nt('Attributes', attr_str)
    XYData = _nt('XYData', 'attr proj lons lats xlons ylats')
    nco = Dataset(nc_cro_file)
    # GDTYP = nco.getncattr('GDTYP')
    # P_GAM = nco.getncattr('P_GAM')
    P_ALP, P_BET = nco.getncattr('P_ALP'), nco.getncattr('P_BET')
    XCENT, YCENT = nco.getncattr('XCENT'), nco.getncattr('YCENT')
    XORIG, YORIG = nco.getncattr('XORIG'), nco.getncattr('YORIG')
    XCELL, YCELL = nco.getncattr('XCELL'), nco.getncattr('YCELL')
    NROWS, NCOLS = nco.getncattr('NROWS'), nco.getncattr('NCOLS')
    proj_str = "+proj=lcc +lon_0={} +lat_0={} +lat_1={} +lat_2={}"
    proj = pyproj.Proj(proj_str.format(XCENT, YCENT, P_ALP, P_BET))
    lons = np.squeeze(nco.variables[lon_name][:])
    lats = np.squeeze(nco.variables[lat_name][:])
    nco.close()
    xlons, ylats = proj(lons, lats)
    attr = attr(XORIG, YORIG, XCENT, YCENT, XCELL, YCELL,
                NROWS, NCOLS, P_ALP, P_BET)
    return XYData(attr, proj, lons, lats, xlons, ylats)


def _get_closest_grid_loc_(lat, lon, xydata):
    from math import floor
    a = xydata.attr
    x, y = xydata.proj(lon, lat)
    ilon, ilat = floor((x - a.XORIG) / a.XCELL), floor((y - a.YORIG) / a.YCELL)
    if ilon < 0 or ilon >= a.NCOLS:
        raise ValueError('lon is out of domain bounds')
    if ilat < 0 or ilat >= a.NROWS:
        raise ValueError('lat is out of domain bounds')
    return Location(ilat, ilon,
                    xydata.lats[ilat, ilon],
                    xydata.lons[ilat, ilon])


def _get_latlon_(proj, dom):
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'

    if isinstance(proj, str):
        proj = s.get_proj_by_name(proj)

    if isinstance(dom, str):
        dom = proj.get_dom_by_name(dom)

    m = proj.months[0]
    month_name = calendar.month_name[m].lower()
    DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size),
                     dom.name, month_name + '_monthly')
    CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name, proj.years[0], m)
    return _get_latlon_from_cro_(_join(DIR_MCIP, CRO_FILE))


def _get_dates_from_nc_(nc_file):
    nco = Dataset(nc_file)
    TFLAG = nco.variables['TFLAG'][:, 0, :]
    nco.close()
    dates = ['{} {:06d}'.format(i[0], i[1]) for i in TFLAG]
    dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
    dt = xr.DataArray(range(0, len(dates)),
                      [("time", pd.to_datetime(dates))])
    return dt


def _get_dates_(proj, dom_name, years, months):
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

    dates = [_get_dates_from_nc_(fmt.format(y, m))
             for y in years for m in months]
    dates = xr.concat(dates, 'time')
    return dates


def _get_data_(dom_names, proj, pol_names, years, months, tstep=True):
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

                xy = _get_latlon_from_cro_(_join(DIR_MCIP, CRO_FILE))
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
                                        'Latitude': (('y', 'x'), xy.lats),
                                        'Longitude': (('y', 'x'), xy.lons)})
                    else:
                        pol = np.squeeze(nco.variables[pol_name][:])
                        yield xr.DataArray(
                            pol, dims=['t', 'y', 'x'],
                            coords={'time': (('t'), dates),
                                    'Latitude': (('y', 'x'), xy.lats),
                                    'Longitude': (('y', 'x'), xy.lons)})
                nco.close()


def _get_data2_(dom_names, proj, pol_names, years, months,
                slice_dates=None, slice_lats=None, slice_lons=None,
                iterate=True):
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

                xy = _get_latlon_from_cro_(_join(DIR_MCIP, CRO_FILE))

                dates = _get_dates_(proj, dn, y, m)
                dates = dates.loc[slice_dates]
                date_values = dates.coords.to_index()
                date_indices = dates.to_series().to_list()

                nco = Dataset(_join(DIR_POST, COM_FILE))

                for pol_name in pol_names:

                    if iterate:
                        for i, k in enumerate(date_indices):
                            pol = np.asarray(nco.variables[pol_name][k][:])
                            yield xr.DataArray(
                                pol, dims=['t', 'y', 'x'],
                                coords={'time': (('t'), [date_values[i]]),
                                        'Latitude': (('y', 'x'), xy.lats),
                                        'Longitude': (('y', 'x'), xy.lons)})
                    else:
                        if len(date_indices) > 0:
                            pol = nco.variables[pol_name][date_indices][:]
                            pol = np.squeeze(pol)
                            yield xr.DataArray(
                                pol, dims=['t', 'y', 'x'],
                                coords={'time': (('t'), date_values),
                                        'Latitude': (('y', 'x'), xy.lats),
                                        'Longitude': (('y', 'x'), xy.lons)})
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


class Domain:
    def __init__(self, name, xydata):
        self.name = name
        self.xydata = xydata

    def __repr__(self):
        s = 'Domain:\n'
        s += '  Name: {}\n'.format(self.name)
        s += '  Size: {} x {}\n'.format(self.xydata.attr.XCELL,
                                        self.xydata.attr.YCELL)
        s += '  ncol: {}\n'.format(self.xydata.attr.NCOLS)
        s += '  nrow: {}\n'.format(self.xydata.attr.NROWS)
        return s

    def closest_grid(self, lat, lon):
        return _get_closest_grid_loc_(lat, lon, self.xydata)

    def closest_grid_hav(self, lat, lon):
        # lon, lat = 27.1633, 38.4217 - İzmir
        lons = self.xydata.lons
        lats = self.xydata.lats
        h = _ma.empty(lats.shape)
        for i in range(lats.shape[0]):
            for j in range(lats.shape[1]):
                h[i, j] = _haversine_(lon, lat, lons[i, j], lats[i, j])
        g = _ma.where(h == h.min())
        result = Location(g[0][0], g[1][0], lats[g][0], lons[g][0])
        return result


class PostProc:
    def __init__(self, proj_names=None, dom_names=None, pol_names=None,
                 years=None, months=None):
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
                    str = "Domain {} is not belong to project {}"
                    raise ValueError(str.format(d, p.name))

        # check domain sizes are consistent among projects
        p0 = next(iter(_proj_.values()))
        doms = {d: p0.get_dom_by_name(d) for d in dom_names}
        dom_size = {k: (d.nrow, d.ncol) for k, d in doms.items()}
        for p in _proj_.values():
            ds = {k: (d.nrow, d.ncol) for k, d in
                  {d: p.get_dom_by_name(d) for d in dom_names}.items()}
            for k, d in dom_size.items():
                if d != ds[k]:
                    str = "Incompatible domain size for {}.{} and {}.{}"
                    raise ValueError(str.format(p0.name, d, p.name, ds[k]))

        doms2 = {k: Domain(k, _get_latlon_(p0, doms[k])) for k in dom_names}

        if pol_names is None:
            pol_names = list(p0.pols.keys())
        if not isinstance(pol_names, list):
            pol_names = [pol_names]

        # check pol names are defined in project settings
        for p in pol_names:
            if p.lower() not in p0.pols.keys():
                str = "Pollutant {} can not be found in settings"
                raise ValueError(str.format(p))

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

        self._dates_ = _get_dates_(p0, dom_names[0], years, months)

        self.proj_names = proj_names
        self._proj_ = _proj_
        # self.doms = doms2
        self.__dict__.update(doms2)
        self.dom_names = dom_names
        self.pols = pols
        self.years = years
        self.months = months
        self.dates = self._dates_.coords['time'].values

    def get_data(self, slice_dates=None):
        ret = []
        for i in self._iterate_(slice_dates, iterate=False):
            ret.append(i)
        if len(ret) > 0:
            keys = list(ret[0].keys())
            ret = {k: [r[k] for r in ret] for k in keys}
            ret = {k: xr.concat(v, 't') for k, v in ret.items()}
        return ret

    def iterate(self, slice_dates=None):
        yield self._iterate_(slice_dates)

    def _iterate_(self, slice_dates=None, iterate=True):
        # years = [2015]
        # months = [1, 2, 3]
        # proj_names = ['cityair', 'cityair_future', 'cityair_future_wama']
        # dom_names = ['aegean', 'central_blacksea', 'mediterranean']
        # pol_names = ['NOX', 'O3', 'CO', 'SO2_UGM3', 'PM10', 'PM25_TOT']

        if slice_dates is None:
            slice_dates = slice(None, None)
        if not isinstance(slice_dates, slice):
            raise ValueError("{} argument must be a slice object")

        # dt = self._dates_
        # dt = dt.loc[slice_dates]
        # idates = dt.to_series().to_list()
        # idates = slice(min(idates), max(idates))

        pol_names = [p.name for p in self.pols.values()]
        dim_names = ['domain', 'project', 'pol_name']
        g = _expandgrid_(self.dom_names, self.proj_names, pol_names)
        g = list(map(tuple, zip(*g)))

        it = {i: _get_data2_(i[0], i[1], i[2], self.years, self.months,
                             slice_dates, iterate=iterate)
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

