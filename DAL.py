#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703
# pylint: disable=C0115,C0116,C0415,E1101
# pylint: disable=R0913

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import calendar as _cal
from datetime import datetime as _dt
from os.path import join as _join
from collections import namedtuple as _nt

import numpy as _np
import numpy.ma as _ma
import pandas as _pd

import netCDF4 as _nc
import xarray as _xr

from _helper_functions_ import check_slice as _check_slice_
from settings import setting as _set


def _get_data2_(dom_names, proj, pol_names,
                slice_dates=None, slice_ilats=None, slice_ilons=None,
                layer_mean=False, gd=None, iterate=True):
    # TODO: This function need to simplify
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'

    # check arg types
    if not isinstance(dom_names, list):
        dom_names = [dom_names]

    if isinstance(proj, str):
        proj = _set.get_proj_by_name(proj)

    if not isinstance(pol_names, list):
        pol_names = [pol_names]

    # Get years/months from project
    # We don't have any reason to supply years/months.
    # We already have slice_dates argument.
    nc_dates = ncDates.from_proj(proj, dom_names[0])
    slice_dates, _ = _check_slice_(slice_dates, (str, int))
    nc_dates = nc_dates.loc(slice_dates)

    slice_ilats, type_ilats = _check_slice_(slice_ilats, (int, float))
    slice_ilons, type_ilons = _check_slice_(slice_ilons, (int, float))

    if (type_ilats, type_ilons) == (float, float):
        if gd is not None:
            mid_lat = (slice_ilats.start + slice_ilats.stop) / 2
            mid_lon = (slice_ilons.start + slice_ilons.stop) / 2
            l1 = gd.nearest_grid_loc(slice_ilats.start, slice_ilons.start)
            l2 = gd.nearest_grid_loc(slice_ilats.stop, slice_ilons.stop)
            ilats = (l1.ilat, l2.ilat)
            ilats = (min(ilats), max(ilats))
            if ilats[0] == ilats[1]:
                ilats = (min(ilats), max(ilats) + 1)
                mid_lat = slice_ilats.start

            ilons = (l1.ilon, l2.ilon)
            ilons = (min(ilons), max(ilons))
            if ilons[0] == ilons[1]:
                ilons = (min(ilons), max(ilons) + 1)
                mid_lon = slice_ilons.start
            slice_ilats = slice(ilats[0], ilats[1])
            slice_ilons = slice(ilons[0], ilons[1])
        else:
            err_msg = 'if lat/lon slices are float, \
            gd (GridData) must be given'
            raise ValueError(err_msg)

    for dn in dom_names:
        dom = proj.get_dom_by_name(dn)
        for y in nc_dates.years:
            for m in nc_dates.months:
                month_name = _cal.month_name[m].lower()
                DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size),
                                 dom.name, month_name + '_monthly')
                DIR_POST = proj.path.post

                CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name, y, m)
                COM_FILE = fmt_com.format(proj.cmaq_ver, proj.compiler,
                                          proj.name, dom.size, dom.name, y, m)

                xy = GridData(_join(DIR_MCIP, CRO_FILE))
                lats = xy.lats[slice_ilats, slice_ilons]
                lons = xy.lons[slice_ilats, slice_ilons]
                # get middle
                k1 = 0 if lats.shape[0] == 2 else (lats.shape[0] // 2)
                k2 = 0 if lats.shape[1] == 2 else (lats.shape[1] // 2)
                mid_lat, mid_lon =  lats[k1, k2], lons[k1, k2]

                dates = ncDates.from_proj(proj, dn, y, m)
                dates = dates.loc(slice_dates)

                nco = _nc.Dataset(_join(DIR_POST, COM_FILE))

                for pol_name in pol_names:

                    if iterate:
                        for i, k in enumerate(dates.indices):
                            pol = nco.variables[pol_name]
                            pol = pol[k, :, slice_ilats, slice_ilons]
                            pol = _xr.DataArray(
                                pol, dims=['t', 'y', 'x'],
                                coords={'time': (('t'), [dates.values[i]]),
                                        'Latitude': (('y', 'x'), lats),
                                        'Longitude': (('y', 'x'), lons)})
                            if layer_mean:
                                pol = pol.mean(['y', 'x'])
                                pol = pol.expand_dims(['y', 'x'],
                                                      [len(pol.dims),
                                                      len(pol.dims) + 1])
                                pol = pol.assign_coords(
                                    Latitude=('y', [mid_lat]),
                                    Longitude=(['x'], [mid_lon]))
                            yield pol
                    else:
                        if len(dates.indices) > 0:
                            pol = nco.variables[pol_name]
                            pol = pol[dates.indices, :, slice_ilats,
                                      slice_ilons]
                            pol = _xr.DataArray(
                                pol, dims=['t', 'l', 'y', 'x'],
                                coords={'time': (('t'), dates.values),
                                        'layer': (('l'), [1]),
                                        'Latitude': (('y', 'x'), lats),
                                        'Longitude': (('y', 'x'), lons)})
                            pol = pol.squeeze('l')
                            pol = pol.drop(labels='layer')
                            if layer_mean:
                                pol = pol.mean(['y', 'x'])
                                pol = pol.expand_dims(['y', 'x'],
                                                      [len(pol.dims),
                                                      len(pol.dims) + 1])
                                pol = pol.assign_coords(
                                    Latitude=('y', [mid_lat]),
                                    Longitude=(['x'], [mid_lon]))
                            yield pol
                nco.close()


def _concat_(x, dim_names, recursive=True):
    if len(dim_names) != len(list(x.keys())[0]):
        msg = "length of x.keys() must be equal to length of dim_names"
        raise ValueError(msg)

    g2 = list(set(i[:-1] for i in x.keys()))
    g2 = [i[0] if len(i) == 1 else i for i in g2]
    doms = {i: [] for i in g2}
    for k in doms:
        for k2, d2 in x.items():
            k2 = k2[:-1]
            if isinstance(k2, tuple):
                if len(k2) == 1:
                    k2 = k2[0]
            if k == k2:
                doms[k].append(d2)

    # combine by last dim
    for k in doms:
        doms[k] = _xr.concat(doms[k], dim_names[-1])

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


class Location:
    def __init__(self, lat, lon, ilat=None, ilon=None, name=None,
                 city=None, region=None):
        # if lat is not None:
        #     lat = float(lat)
        # if lon is not None:
        #     lon = float(lon)
        # if ilat is not None:
        #     ilat = int(ilat)
        # if ilon is not None:
        #     ilon = int(ilon)
        d = {'lat': lat, 'lon': lon, 'ilat': ilat,
             'ilon': ilon, 'name': name, 'city': city, 'region': region}
        self.__dict__.update(d)

    def __repr__(self):
        name = f'({self.name})' if self.name is not None else ''
        latlon = ' '.join([f"{k}: {v:.5f}" for k, v in self.__dict__.items()
                         if v is not None and k in ['lat', 'lon']])
        ilatlon = ' '.join([f"{k}: {v}" for k, v in self.__dict__.items()
                         if v is not None and k in ['ilat', 'ilon']])
        return f'Location:{name} {{' + latlon + ' ' + ilatlon + '}'

    def is_in(self, x):
        if not isinstance(x, (Domain, GridData)):
            raise ValueError('x must be a Domain or GridData object')
        return x.contains(self)


class ncDates:
    def __init__(self, xarray, nc_file=None):
        self.xarr = xarray
        self.values = xarray.coords.to_index()
        self.indices = xarray.to_series().to_list()
        self.years = list(set(i.year for i in self.values))
        self.months = list(set(i.month for i in self.values))
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
            slice_dates = slice(str(_np.datetime_as_string(sliced.min())),
                                str(_np.datetime_as_string(sliced.max())))
        return ncDates(self.xarr.loc[slice_dates])

    @classmethod
    def _get_dates_from_dates_(cls, dates):
        xarr = _xr.DataArray(range(0, len(dates)),
                             [("time", _pd.to_datetime(dates))])
        return xarr

    @classmethod
    def _get_dates_from_nc_(cls, nc_file):
        nc = _nc.Dataset(nc_file)
        TFLAG = nc.variables['TFLAG'][:, 0, :]
        nc.close()
        dates = [f'{i[0]} {i[1]:06d}' for i in TFLAG]
        dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
        return cls._get_dates_from_dates_(dates)

    @classmethod
    def from_dates(cls, dates):
        return cls(cls._get_dates_from_dates_(dates))

    @classmethod
    def from_nc_file(cls, nc_file):
        return cls(cls._get_dates_from_nc_(nc_file), nc_file)

    @classmethod
    def from_proj(cls, proj, dom_name, years=None, months=None):
        if years is None:
            years = proj.years
        if months is None:
            months = proj.months
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

        list_of_nc_files = [fmt.format(y, m) for y in years 
                            for m in months]
        xarr = [cls._get_dates_from_nc_(f) for f in list_of_nc_files]
        xarr = _xr.concat(xarr, 'time')
        return cls(xarr, list_of_nc_files)


class GridData:
    def __init__(self, nc_cro_file, lat_name='LAT', lon_name='LON'):
        from pyproj import Proj
        attributes = ('XORIG YORIG XCENT YCENT XCELL YCELL NROWS NCOLS ' +
                      'P_ALP P_BET GDTYP P_GAM').split(' ')
        variables = 'proj lons lats xlons ylats'.split(' ')

        nco = _nc.Dataset(nc_cro_file)
        atts = {k: nco.getncattr(k) for k in attributes}
        lons = _np.squeeze(nco.variables[lon_name][:])
        lats = _np.squeeze(nco.variables[lat_name][:])
        nco.close()

        # pylint: disable=C0209
        proj = Proj("+proj=lcc +lon_0={XCENT} +lat_0={YCENT} \
            lat_1={P_ALP} +lat_2={P_BET}".format(**atts))
        xlons, ylats = proj(lons, lats)  # pylint: disable=E0633
        data = dict(zip(variables, [proj, lons, lats, xlons, ylats]))
        self.nc_cro_file = nc_cro_file
        self.__dict__.update(data)
        self.__dict__.update(atts)
        Bounds = _nt('Bounds', ['lat', 'lon'])
        self.bounds = Bounds([float(self.lats.min()), float(self.lats.max())],
                             [float(self.lons.min()), float(self.lons.max())])
        self.ibounds = Bounds([0, int(self.NROWS - 1)],
                              [0, int(self.NCOLS - 1)])

    def __repr__(self):
        s = 'Grid:\n'
        s += f'  SIZE                : {self.XCELL} x {self.YCELL}\n'
        s += f'  NCOL x NROW         : {self.NCOLS} x {self.NROWS}\n'
        s += f'  XORIG, YORIG        : {self.XORIG}, {self.YORIG}\n'
        s += f'  XCENT, YCENT        : {self.XCENT}, {self.YCENT}\n'
        s += f'  P_ALP, P_BET, P_GAM : {self.P_ALP}, {self.P_BET},'
        s += f' {self.P_GAM}\n'
        s += f'  GDTYP               : {self.GDTYP}\n'
        s += 'BOUNDS:\n'
        s += '        MIN    MAX\n'
        s += f'  LAT: {self.bounds.lat[0]:.4f} {self.bounds.lat[1]:.4f}\n'
        s += f'  LON: {self.bounds.lon[0]:.4f} {self.bounds.lon[1]:.4f}\n'
        return s

    def get_lat_at(self, i, j):
        return self.lats[i, j]

    def get_lon_at(self, i, j):
        return self.lons[i, j]

    def get_latlon_at(self, i, j):
        return (self.lats[i, j], self.lons[i, j])

    def nearest_grid_loc(self, lat, lon, name=None, raise_error=True):
        from math import floor
        x, y = self.proj(lon, lat)
        ilon, ilat = (floor((x - self.XORIG) / self.XCELL),
                      floor((y - self.YORIG) / self.YCELL))
        if raise_error:
            if ilon < 0 or ilon >= self.NCOLS:
                raise ValueError('lon is out of domain bounds')
            if ilat < 0 or ilat >= self.NROWS:
                raise ValueError('lat is out of domain bounds')
        else:
            ilat = max(ilat, 0)
            ilat = min(ilat, self.NROWS - 1)
            ilon = max(ilon, 0)
            ilon = min(ilon, self.NCOLS - 1)
        # lat = float(self.lats[ilat, ilon])
        # lon = float(self.lons[ilat, ilon])
        lat = self.lats[ilat, ilon]
        lon = self.lons[ilat, ilon]
        return Location(lat, lon, ilat, ilon, name)

    def nearest_grid_hav(self, lat, lon, name=None):
        # lon, lat = 27.1633, 38.4217 - İzmir
        lons = self.lons
        lats = self.lats
        h = _ma.empty(lats.shape)
        for i in range(lats.shape[0]):
            for j in range(lats.shape[1]):
                h[i, j] = _haversine_(lon, lat, lons[i, j], lats[i, j])
        g = _ma.where(h == h.min())
        return Location(lats[g][0], lons[g][0], g[0][0], g[1][0], name)

    def contains(self, loc):
        if isinstance(loc, list):
            return [self.contains(lo) for lo in loc]
        try:
            self.nearest_grid_loc(loc.lat, loc.lon)
            return True
        except ValueError:
            return False

    @classmethod
    def from_dom(cls, proj, dom):
        fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'

        if isinstance(proj, str):
            proj = _set.get_proj_by_name(proj)

        if isinstance(dom, str):
            dom = proj.get_dom_by_name(dom)

        m = proj.months[0]
        month_name = _cal.month_name[m].lower()
        DIR_MCIP = _join(proj.path.mcip, f'{dom.size}km',
                         dom.name, f'{month_name}_monthly')
        CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name,
                                  proj.years[0], m)
        return GridData(_join(DIR_MCIP, CRO_FILE))


class Domains:
    def __init__(self, dic):
        self.__mem__ = dic
        self.__dict__.update(dic)

    def __repr__(self):
        return 'Domains:\n' + '\n'.join([f" * {k}" for k in self.__mem__])

    def __iter__(self):
        for i in self.__mem__.values():
            yield i


class Domain:
    def __init__(self, name, gd, pp):
        from copy import deepcopy
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

    def contains(self, loc):
        return self.gd.contains(loc)

    def nearest_grid_loc(self, loc, name=None):
        return self.gd.nearest_grid_loc(loc.lat, loc.lon, name,
                                        raise_error=False)

    def get_data_loc(self, slice_dates=None, loc=None, delta=0,
                     layer_mean=False, simplify=True):
        ret = None
        if not isinstance(delta, int):
            raise ValueError('delta must be integer')
        if delta < 0:
            raise ValueError('delta must be positive integer')
        if isinstance(loc, list):
            for i, lo in enumerate(loc):
                if lo.name is None:
                    lo.name = f'S{i}'
            data = [self.get_data_loc(slice_dates, lo, delta, layer_mean,
                                      simplify) for lo in loc]
            if len(data) > 1:
                data = _xr.concat(data, dim='sta')
            else:
                data = data[0]
            ret = data
        if isinstance(loc, Location):
            if loc.name is None:
                import random
                n = random.randint(0, 99999)
                loc.name = f'S{n:05d}'
            slice_ilats = slice(loc.lat, loc.lat)
            slice_ilons = slice(loc.lon, loc.lon)
            if (type(loc.lat), type(loc.lon)) == (float, float):
                if self.gd is not None:
                    l1 = self.nearest_grid_loc(loc)
                    ilat_start, ilat_end = l1.ilat - delta, l1.ilat + delta + 1
                    ilat_start = max(ilat_start, 0)
                    ilat_start = min(ilat_start, self.gd.ibounds.lat[1])
                    ilat_end = max(ilat_end, 0)
                    ilat_end = min(ilat_end, self.gd.ibounds.lat[1])
                    if ilat_start == ilat_end:
                        ilat_end += 1
                    slice_ilats = slice(ilat_start, ilat_end)

                    ilon_start, ilon_end = l1.ilon - delta, l1.ilon + delta + 1
                    ilon_start = max(ilon_start, 0)
                    ilon_start = min(ilon_start, self.gd.ibounds.lon[1])
                    ilon_end = max(ilon_end, 0)
                    ilon_end = min(ilon_end, self.gd.ibounds.lon[1])
                    if ilon_start == ilon_end:
                        ilon_end += 1
                    slice_ilons = slice(ilon_start, ilon_end)
                else:
                    err_msg = 'if loc.lat/loc.lon are float, \
                    gd (GridData) must be given'
                    raise ValueError(err_msg)
            ret = self.get_data(slice_dates, slice_ilats,
                                slice_ilons, layer_mean, simplify)
            ret = ret.squeeze(dim='domain')
            if simplify and isinstance(ret, _xr.DataArray):
                ret = _np.squeeze(ret)
            ret = ret.expand_dims({'sta': 1})
            ret = ret.assign_coords(sta=[loc.name])
        if ret is not None:
            return ret
        err_msg = 'loc must be Location object or list of Location objects'
        raise ValueError(err_msg)

    def get_data(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
                 layer_mean=False, simplify=True):
        return self._pp_.get_data(slice_dates, slice_ilats, slice_ilons,
                                  layer_mean, simplify)

    def iter(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
             layer_mean=False, simplify=True):
        for i in self._pp_.iter(slice_dates, slice_ilats, slice_ilons,
                                layer_mean, simplify):
            yield i


class PostProc:

    def __init__(self, proj_names=None, dom_names=None, pol_names=None):
        years, months = None, None
        if proj_names is None:
            proj_names = _set.get_proj_names()
        if not isinstance(proj_names, list):
            proj_names = [proj_names]
        _proj_ = _set.get_proj_by_name(proj_names)

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
            proj_pols = [p.lower() for p in list(p0.pols.keys())]
            if p.lower() not in proj_pols:
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

    def __iter__(self):
        for i in self.iter():
            yield i

    def get_data(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
                 layer_mean=False, simplify=True):
        ret = list(self._iter_(slice_dates, slice_ilats,
                               slice_ilons, layer_mean, simplify,
                               iterate=False))
        if len(ret) > 1:
            ret = {k: [r[k] for r in ret] for k in ret[0].keys()}
            ret = {k: _xr.concat(v, 't') for k, v in ret.items()}
        else:
            ret = ret[0]
        if isinstance(ret, dict) and len(ret) == 1:
            ret = list(ret.values())[0]
        return ret

    def iter(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
             layer_mean=False, simplify=True):
        for i in self._iter_(slice_dates, slice_ilats, slice_ilons,
                             layer_mean, simplify):
            yield i

    def _iter_(self, slice_dates=None, slice_ilats=None, slice_ilons=None,
               layer_mean=False, simplify=True, iterate=True):
        def expandgrid(*itrs):
            import itertools
            product = list(itertools.product(*itrs))
            return [[x[i] for x in product] for i in range(len(itrs))]

        pol_names = [p.name for p in self.pols.values()]
        dim_names = ['domain', 'project', 'pol_name']
        g = expandgrid(self.dom_names, self.proj_names, pol_names)
        g = list(map(tuple, zip(*g)))

        it = {i: _get_data2_(i[0], i[1], i[2],
                             slice_dates, slice_ilats, slice_ilons,
                             layer_mean, self._gd_, iterate)
              for i in g}
        counter = 0
        while True:
            try:
                d = {k: next(v) for k, v in it.items()}
                for k, v in d.items():
                    d[k] = v.assign_coords(dict(zip(dim_names, k)))
                    d[k] = d[k].expand_dims(dict(zip(dim_names, [1, 1, 1])))

                d = _concat_(d, dim_names)
                if simplify:
                    for k in d:
                        d[k] = _np.squeeze(d[k])
                yield d

            except StopIteration:
                break
            counter = counter + 1
