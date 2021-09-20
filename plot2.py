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
import string
from copy import copy
import numpy as np
import pandas as pd
from netCDF4 import Dataset  # pylint: disable=E0611
import xarray as xr
# from xarray.plot import pcolormesh as pcm
import cartopy.crs as ccrs
# import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely.geometry as sgeom
import cmocean
from settings import setting as s

import matplotlib
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap as lscmap
matplotlib.use('Agg')


proj = s.get_active_proj()


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def get_latlon_from_cro(nc_cro_file, lat_name='LAT', lon_name='LON'):
    """
    Get latitude and longitude values from CMAQ CRO files
    """
    nco = Dataset(nc_cro_file)
    lon = np.squeeze(nco.variables[lon_name][:])
    lat = np.squeeze(nco.variables[lat_name][:])
    nco.close()

    return lon, lat


def get_lat_lon(dom, month):
    """
    Get lat lon from domain
    """
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    month_name = calendar.month_name[month].lower()
    DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size),
                     dom.name, month_name + '_monthly')
    CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name,
                              year, month)
    nc_cro_file = _join(DIR_MCIP, CRO_FILE)
    return get_latlon_from_cro(nc_cro_file)


def get_dates_from_nco(nco):
    """ Get dates from NCO """
    TFLAG = np.squeeze(nco.variables['TFLAG'][:])
    tf = TFLAG.reshape(TFLAG.shape[0] * TFLAG.shape[1], 2)
    tf = np.unique(tf, axis=0)
    dates = ['{} {:06d}'.format(i[0], i[1]) for i in tf]
    dates = [_dt.strptime(i, '%Y%j %H%M%S') for i in dates]
    return pd.date_range(str(min(dates)), str(max(dates)),
                         freq='H')


def get_days_from_nco(nco):
    """ Get days from NCO file """
    TFLAG = np.squeeze(nco.variables['TFLAG'][:])
    tf = TFLAG.reshape(TFLAG.shape[0] * TFLAG.shape[1], 2)
    tf = np.unique(tf, axis=0)
    days = [_dt.strptime(str(i[0]), '%Y%j') for i in tf]
    return pd.date_range(min(days), max(days))


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


def plot_map(doms, path, cmap='twilight', cb_levels=None,
             rast_zorder=None, cb_limits=None):
    """
    Plot datasets in doms object.
    """
    import matplotlib.pyplot as plt  # pylint: disable=C0415
    import matplotlib.transforms as mtransforms  # pylint: disable=C0415
    for dom_name, d in doms.items():
        for i, a in enumerate(d.transpose('pol_name', ...)):
            pol_name = a.coords['pol_name'].values.tolist()
            x = a.drop(labels=['pol_name'])

            _, _, month, stat = (i[1].values for i in x.coords.items())

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
            levels = None if cb_levels is None else cb_levels[pol_name]
            if cb_limits is None:
                p = x.plot(x='Longitude', y='Latitude', col='month',
                           row='stat', robust=True, size=4,
                           aspect=x.shape[3] / x.shape[2],
                           cbar_kwargs=cbar_kws, cmap=cmap, levels=levels,
                           subplot_kws=dict(projection=ccrs_proj),
                           transform=ccrs.PlateCarree(), zorder=0,
                           alpha=1.0)
            else:
                vmin = min(cb_limits[pol_name])
                vmax = max(cb_limits[pol_name])
                p = x.plot(x='Longitude', y='Latitude', col='month',
                           row='stat', robust=True, size=4,
                           aspect=x.shape[3] / x.shape[2],
                           cbar_kwargs=cbar_kws, cmap=cmap, levels=levels,
                           vmin=vmin, vmax=vmax,
                           subplot_kws=dict(projection=ccrs_proj),
                           transform=ccrs.PlateCarree(), zorder=0,
                           alpha=1.0)

            for i, st in enumerate(stat):
                p.row_labels[i].set_text(st)
            for i, m in enumerate(month):
                p.col_labels[i].set_text(calendar.month_name[m].title())

            for i, ax in enumerate(p.axes.flat):
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

                # pylint: disable=C0415
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
            file_name = '_'.join(('plot', dom_name, pol_name)) + '.pdf'
            p.fig.savefig(_join(path, file_name), bbox_inches='tight')
            plt.close(p.fig)
            print(file_name)


def calc_stat(dom_names, pol_names, year, months, stats=None):
    """
    Calculate statistics for domains.
    """
    STATS = {'mean': {'day': 'mean', 'mon': 'mean'},
             'daily_max': {'day': 'mean', 'mon': 'max'},
             'hourly_max': {'day': 'max', 'mon': 'max'}}
    stats = stats or ['mean', 'daily_max', 'hourly_max']
    if not isinstance(stats, list):
        stats = list(stats)
    for st in stats:
        assert st in list(STATS.keys()), \
            'stats arg must be mean, daily_max or hourly_max.'

    statistics = ['mean', 'max', 'min']
    for k in stats:
        assert STATS[k]['day'] in statistics, \
            'stat_day arg must be mean, max or min.'
        assert STATS[k]['mon'] in statistics, \
            'stat_mon arg must be mean, max or min.'

    def get_monthly_stats(x, stats, days):
        xm2 = None
        for k in stats:
            stat_day = STATS[k]['day']
            stat_mon = STATS[k]['mon']
            xd = getattr(x.groupby('time.day'), stat_day)()
            xd = xd.assign_coords(day=days)
            xm = getattr(xd.groupby('day.month'), stat_mon)()
            xm2 = xm if xm2 is None else xr.concat([xm2, xm], dim='stat')
        return xm2.assign_coords(stat=stats)

    def get_monthly_stats_for_file(dom_name, pol_names, months):
        fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'
        pols = {i: None for i in pol_names}
        dom = proj.get_dom_by_name(dom_name)
        for m in months:
            nco = Dataset(_join(proj.path.post,
                          fmt_com.format(
                            proj.cmaq_ver, proj.compiler,
                            proj.name, dom.size, dom.name,
                            year, m)))
            lon, lat = get_lat_lon(dom, m)
            dates = get_dates_from_nco(nco)
            days = get_days_from_nco(nco)

            for pn in pols.keys():
                x = xr.DataArray(np.squeeze(nco.variables[pn][:]),
                                 dims=['t', 'y', 'x'],
                                 coords={'time': (('t'), dates),
                                         'Latitude': (('y', 'x'), lat),
                                         'Longitude': (('y', 'x'), lon)})
                x = get_monthly_stats(x, stats, days)
                pols[pn] = x if pols[pn] is None else xr.concat([pols[pn], x],
                                                                dim='month')

            nco.close()
        return pols

    doms = {i: None for i in dom_names}
    for dn in doms.keys():
        pols = get_monthly_stats_for_file(dn, pol_names, months)
        doms[dn] = xr.concat(list(pols.values()),
                             pd.Index(list(pols.keys()), name="pol_name"))
        doms[dn].attrs['units'] = 'ugm-3'
    return doms


def print_stats(doms):
    """
    Print statistics for domains
    """
    for dom_name, x in doms.items():
        print(dom_name)
        for _, x in enumerate(x.transpose('pol_name', ...)):
            pol_name = x.coords['pol_name'].values.tolist()
            print(pol_name)
            x = x.drop(labels=['pol_name'])
            for _, x in enumerate(x.transpose('stat', ...)):
                stat = x.coords['stat'].values.tolist()
                print(' ' + stat)
                x = x.drop(labels=['stat'])
                for _, x in enumerate(x.transpose('month', ...)):
                    month = x.coords['month'].values.tolist()
                    print('  ' + str(month))
                    x = x.drop(labels=['month'])
                    print('   Min  : {:.2f}'.format(float(x.min())))
                    print('   Mean : {:.2f}'.format(float(x.mean())))
                    print('   Max  : {:.2f}'.format(float(x.max())))


# c = mcolors.ColorConverter().to_rgb
# rvb = make_colormap(
#     [c('red'), c('violet'), 0.33, c('violet'), c('blue'), 0.66, c('blue')])

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
# cmap = cmocean.cm.thermal_r  # pylint: disable=E1101
# cmap = 'twilight'
cmap = lscmap.from_list("",
                        ['#E6E1E7', '#d9ca92', '#77c0ab', '#7581C0',
                         '#633A90', '#3B2044', '#7F3660', '#b266a9',
                         '#928149', '#D3B8A5', '#E6E1E7'])
cmap.name = 'mycmap'

cmap = matplotlib.colors.ListedColormap(['#7dede6', '#71c9ac', '#efe661',
                                         '#ee5c57', '#8b1a34', '#742c7d'])
cmap.name = 'discrete'
NO_LIMIT = True
LEVELS = True

cmap_str = cmap if isinstance(cmap, str) else cmap.name
no_limit_str = 'no_limit' if NO_LIMIT else ''
pth = _join('plots_combined', '_'.join(('plots', cmap_str, no_limit_str)))

if NO_LIMIT:
    CB_LIMITS = None
else:
    CB_LIMITS = {'CO': [0, 100],
                 'NOX': [0, 10],
                 'O3': [0, 60],
                 'PM10': [0, 80],
                 'PM25_TOT': [0, 40],
                 'SO2_UGM3': [0, 20]}

if LEVELS:
    CB_LEVELS = {'CO': [0, 10, 20, 25, 50, 75, 800],
                 'NOX': [0, 40, 90, 120, 230, 340, 1000],
                 'O3': [0, 50, 100, 130, 240, 380, 800],
                 'PM10': [0, 20, 40, 50, 100, 150, 1200],
                 'PM25_TOT': [0, 10, 20, 25, 50, 75, 800],
                 'SO2_UGM3': [0, 100, 200, 350, 500, 750, 1250]}
else:
    CB_LEVELS = None

doms = calc_stat(DOM_NAMES, POL_NAMES, year, months)
plot_map(doms, pth, cmap, CB_LEVELS, cb_limits=CB_LIMITS)
