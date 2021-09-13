#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create plots for CMAQ
~~~~~~~
Python script to create plots
"""
import warnings
from os import makedirs as _mkdir
from os.path import join as _join
import calendar
from datetime import datetime as _dt
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr
from xarray.plot import pcolormesh as pcm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
from copy import copy
import string
import cmocean
from settings import setting as s
import matplotlib
matplotlib.use('Agg')

proj = s.get_active_proj()

proj = s.get_active_proj()


def get_latlon_from_cro(nc_cro_file, lat_name='LAT', lon_name='LON'):
    """
    Get latitude and longitude values from CMAQ CRO files
    """
    nco = Dataset(nc_cro_file)
    lon = np.squeeze(nco.variables[lon_name][:])
    lat = np.squeeze(nco.variables[lat_name][:])
    nco.close()

    return lon, lat


def lambert_ticks(ax, ticks, axis='x'):
    """Draw ticks on the left or bottom axis of a Lamber Conformal
       projection."""

    def _lambert_ticks_(ax, ticks, tick_location, line_constructor,
                        tick_extractor):
        """Get the tick locations and labels for an axis of a Lambert
           Conformal projection."""

        def find_side(ls, side):
            """ Given a shapely LineString which is assumed to be rectangular,
            return the line corresponding to a given side of the rectangle. """

            minx, miny, maxx, maxy = ls.bounds
            points = {'left': [(minx, miny), (minx, maxy)],
                      'right': [(maxx, miny), (maxx, maxy)],
                      'bottom': [(minx, miny), (maxx, miny)],
                      'top': [(minx, maxy), (maxx, maxy)]}
            return sgeom.LineString(points[side])

        outline_patch = sgeom.LineString(
            ax.outline_patch.get_path().vertices.tolist())
        axis = find_side(outline_patch, tick_location)
        n_steps = 30
        extent = ax.get_extent(ccrs.PlateCarree())
        _ticks = []
        for t in ticks:
            xy = line_constructor(t, n_steps, extent)
            proj_xyz = ax.projection.transform_points(ccrs.Geodetic(),
                                                      xy[:, 0], xy[:, 1])
            xyt = proj_xyz[..., :2]
            ls = sgeom.LineString(xyt.tolist())
            locs = axis.intersection(ls)
            if not locs:
                tick = [None]
            else:
                tick = tick_extractor(locs.xy)
            _ticks.append(tick[0])

        ticklabels = copy(ticks)
        while True:
            try:
                index = _ticks.index(None)
            except ValueError:
                break
            _ticks.pop(index)
            ticklabels.pop(index)
        return _ticks, ticklabels

    def te(xy):
        return xy[0] if axis == 'x' else xy[1]

    def lc(t, n, b):
        npz = np.zeros(n) + t
        if axis == 'x':
            npv = (npz, np.linspace(b[2], b[3], n))
        else:
            npv = (np.linspace(b[0], b[1], n), npz)
        return np.vstack(npv).T

    side = 'bottom' if axis == 'x' else 'left'
    ticks, ticklabels = _lambert_ticks_(ax, ticks, side, lc, te)
    axiss = getattr(ax, '{}axis'.format(axis))
    tick_side = getattr(axiss, 'tick_{}'.format(side))
    tick_side()
    set_ticks = getattr(ax, 'set_{}ticks'.format(axis))
    set_ticks(ticks)
    set_ticklabels = getattr(ax, 'set_{}ticklabels'.format(axis))
    set_ticklabels([axiss.get_major_formatter()(tick) for tick in ticklabels])


def plot_map(doms, path, suffix='', cmap='twilight', rast_zorder=None,
             cb_limits=None):
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    for dom_name, d in doms.items():
        for i, a in enumerate(d.transpose('pol_name', ...)):
            pol_name = a.coords['pol_name'].values.tolist()
            x = a.drop(labels=['pol_name'])

            _, _, month = (i[1].values for i in x.coords.items())
            facet_labels = ['{}'.format(calendar.month_name[i].title())
                            for i in month]

            # central_longitude=24, central_latitude=45,
            ccrs_proj = ccrs.LambertConformal(
                central_longitude=24, central_latitude=45,
                false_easting=400000, false_northing=400000,
                standard_parallels=(46, 49))
            cbar_kws = dict(label=' '.join((POL_LABELS[i], POL_UNITS[i])),
                            pad=0.02, shrink=0.8)
            res = '10m'
            grid_interval = 1
            if dom_name == 'tr':
                grid_interval = 2
            elif dom_name == 'eu':
                grid_interval = 3
            xticks = list(np.arange(-180, 180, grid_interval))
            yticks = list(np.arange(-90, 90, grid_interval))

            # with warnings.catch_warnings():
            #     warnings.simplefilter("ignore")
            if cb_limits is None:
                p = x.plot(x='Longitude', y='Latitude', col='month',
                           robust=True, col_wrap=3, size=4,
                           aspect=x.shape[2] / x.shape[1],
                           cbar_kwargs=cbar_kws, cmap=cmap,
                           subplot_kws=dict(projection=ccrs_proj),
                           transform=ccrs.PlateCarree(), zorder=0,
                           alpha=1.0)
            else:
                vmin = min(cb_limits[pol_name])
                vmax = max(cb_limits[pol_name])
                p = x.plot(x='Longitude', y='Latitude', col='month',
                           robust=True, col_wrap=3, size=4,
                           aspect=x.shape[2] / x.shape[1],
                           cbar_kwargs=cbar_kws, cmap=cmap,
                           vmin=vmin, vmax=vmax,
                           subplot_kws=dict(projection=ccrs_proj),
                           transform=ccrs.PlateCarree(), zorder=0,
                           alpha=1.0)

            for i, ax in enumerate(p.axes.flat):
                ax.set_title(facet_labels[i])
                trans = mtransforms.ScaledTranslation(
                    10 / 72, -5 / 72, ax.figure.dpi_scale_trans)
                ax.text(-0.02, 1.005, '{})'.format(string.ascii_lowercase[i]),
                        transform=ax.transAxes + trans,
                        fontsize='large', verticalalignment='top',
                        fontfamily='serif',
                        bbox=dict(facecolor='none', edgecolor='none',
                                  pad=2.5))
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.coastlines(resolution=res, alpha=0.6)
                # ax.add_feature(cfeature.BORDERS.with_scale(res),
                #                linestyle=':', alpha=1)
                # ax.add_feature(cfeature.OCEAN.with_scale(res),
                #                facecolor=("lightblue"))
                # ax.add_feature(cfeature.LAND.with_scale(res),
                #                facecolor=("peachpuff"))
                p.fig.canvas.draw()

                gl = ax.gridlines(xlocs=xticks, ylocs=yticks,
                                  dms=True, color='indigo', alpha=0.5,
                                  linestyle='--')
                # gl = ax.gridlines(xlocs=xticks, ylocs=yticks, linewidth=1.0,
                #                   color='black', alpha=0.5, linestyle='--')
                # gl.top_labels = gl.right_labels = False
                # gl.rotate_labels = False

                ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
                ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
                # lambert_xticks(ax, xticks)
                # lambert_yticks(ax, yticks)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    lambert_ticks(ax, xticks)
                    lambert_ticks(ax, yticks, 'y')
                if rast_zorder is not None:
                    ax.set_rasterization_zorder(rast_zorder)
                # ax.set_rasterization_zorder(1)
            # p.fig.subplots_adjust(hspace=0, wspace=0)
            # cbaxes = xplot.fig.add_axes([0.8, 0.1, 0.03, 0.8])
            # cbar = plt.colorbar(xplot.fig, label='new cbar')
            # p.fig.show()

            _mkdir(path, exist_ok=True)
            file_name = '_'.join(('plot', dom_name, pol_name, suffix)) + '.pdf'
            p.fig.savefig(_join(path, file_name), bbox_inches='tight')
            plt.close(p.fig)
            print(file_name)


def calc_stat(dom_names, pol_names, months, stat_day='mean', stat_mon='mean'):
    assert stat_day in ['mean', 'max', 'min'], \
        'stat_day arg must be mean, max or min.'
    assert stat_mon in ['mean', 'max', 'min'], \
        'stat_mon arg must be mean, max or min.'
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'
    doms = {i: None for i in dom_names}
    for dn in dom_names:
        pols = {i: None for i in pol_names}
        dom = proj.get_dom_by_name(dn)
        for m in months:
            month_name = calendar.month_name[m].lower()
            DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size), dom.name,
                             month_name + '_monthly')
            DIR_POST = proj.path.post

            CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name, year, m)
            COM_FILE = fmt_com.format(proj.cmaq_ver, proj.compiler, proj.name,
                                      dom.size, dom.name, year, m)

            nc_cro_file = _join(DIR_MCIP, CRO_FILE)
            nc_com_file = _join(DIR_POST, COM_FILE)
            lon, lat = get_latlon_from_cro(nc_cro_file)
            nco = Dataset(nc_com_file)

            for pol_name in pol_names:
                TFLAG = np.squeeze(nco.variables['TFLAG'][:])
                poll = np.squeeze(nco.variables[pol_name][:])
                tf = TFLAG.reshape(TFLAG.shape[0] * TFLAG.shape[1], 2)
                tf = np.unique(tf, axis=0)
                dates = ['{} {:06d}'.format(i[0], i[1]) for i in tf]
                dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
                dates = pd.date_range(str(min(dates)), str(max(dates)),
                                      freq='H')

                x = xr.DataArray(np.squeeze(poll[:, :, :]),
                                 dims=['t', 'y', 'x'],
                                 coords={'time': (('t'), dates),
                                         'Latitude': (('y', 'x'), lat),
                                         'Longitude': (('y', 'x'), lon)})

                xd = getattr(x.groupby('time.day'), stat_day)()

                days = [_dt.strptime(str(i[0]), '%Y%j') for i in tf]
                days = pd.date_range(min(days), max(days))
                xd = xd.assign_coords(day=days)
                xd.attrs['long_name'] = pol_name

                xm = getattr(xd.groupby('day.month'), stat_mon)()

                if pols[pol_name] is None:
                    pols[pol_name] = xm
                else:
                    pols[pol_name] = xr.concat([pols[pol_name], xm],
                                               dim='month')
            nco.close()
        doms[dn] = xr.concat([i for i in pols.values()],
                             pd.Index(list(pols.keys()), name="pol_name"))
        doms[dn].attrs['units'] = 'ugm-3'
    return doms


POL_NAMES = ['NOX', 'O3', 'CO', 'SO2_UGM3', 'PM10', 'PM25_TOT']
DOM_NAMES = ['tr', 'aegean', 'central_blacksea', 'mediterranean',
             'south_central_anatolia']
POL_LABELS = ['$NO_x$', '$O_3$', '$CO$', '$SO_2$', r'$PM_{10}$',
              r'$PM_{2.5}$']
POL_UNITS = ['$(\\mu g/m^3)$', '$(\\mu g/m^3)$',
             '$(\\mu g/m^3)$', '$(\\mu g/m^3)$', '$(\\mu g/m^3)$',
             '$(\\mu g/m^3)$']
year = 2015
months = [1, 2, 3]
cmap = cmocean.cm.thermal_r
# cmap = 'twilight'
NO_LIMIT = False

cmap_str = cmap if isinstance(cmap, str) else cmap.name
no_limit_str = 'no_limit' if NO_LIMIT else ''
pth = _join('plots_single', '_'.join(('plots', cmap_str, no_limit_str)))

if NO_LIMIT:
    CB_LIMITS_mean = None
    CB_LIMITS_daily_max = None
    CB_LIMITS_hourly_max = None
else:
    # CB_LIMITS_mean = {'CO': [0, 100],
    #                   'NOX': [0, 10],
    #                   'O3': [0, 60],
    #                   'PM10': [0, 80],
    #                   'PM25_TOT': [0, 40],
    #                   'SO2_UGM3': [0, 20]}
    # CB_LIMITS_daily_max = {'CO': [0, 200],
    #                        'NOX': [0, 20],
    #                        'O3': [0, 80],
    #                        'PM10': [0, 200],
    #                        'PM25_TOT': [0, 100],
    #                        'SO2_UGM3': [0, 30]}
    # CB_LIMITS_hourly_max = {'CO': [0, 400],
    #                         'NOX': [0, 40],
    #                         'O3': [0, 80],
    #                         'PM10': [0, 300],
    #                         'PM25_TOT': [0, 150],
    #                         'SO2_UGM3': [0, 50]}
    CB_LIMITS_mean = {'CO': [70, 100],
                      'NOX': [0, 10],
                      'O3': [0, 60],
                      'PM10': [0, 80],
                      'PM25_TOT': [0, 40],
                      'SO2_UGM3': [0, 20]}
    CB_LIMITS_daily_max = {'CO': [0, 200],
                           'NOX': [0, 20],
                           'O3': [0, 80],
                           'PM10': [0, 200],
                           'PM25_TOT': [0, 100],
                           'SO2_UGM3': [0, 30]}
    CB_LIMITS_hourly_max = {'CO': [0, 400],
                            'NOX': [0, 40],
                            'O3': [0, 80],
                            'PM10': [0, 300],
                            'PM25_TOT': [0, 150],
                            'SO2_UGM3': [0, 50]}


doms = calc_stat(DOM_NAMES, POL_NAMES, months, stat_day='mean',
                 stat_mon='mean')
plot_map(doms, pth, 'mean', cmap, cb_limits=CB_LIMITS_mean)

doms = calc_stat(DOM_NAMES, POL_NAMES, months, stat_day='mean', stat_mon='max')
plot_map(doms, pth, 'daily_max', cmap, cb_limits=CB_LIMITS_daily_max)

doms = calc_stat(DOM_NAMES, POL_NAMES, months, stat_day='max', stat_mon='max')
plot_map(doms, pth, 'hourly_max', cmap, cb_limits=CB_LIMITS_hourly_max)

