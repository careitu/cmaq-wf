#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create plots for CMAQ
~~~~~~~
Python script to create plots
"""
import calendar
import numpy as np
import warnings

from datetime import datetime as _dt
from os import makedirs as _mkdir
from os.path import isfile as _isfile
from os.path import join as _join

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import shapely.geometry as sgeom
import string
import xarray as xr

from cartopy.mpl.gridliner import LATITUDE_FORMATTER
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER
from copy import copy
from netCDF4 import Dataset
from xarray.plot import pcolormesh as pcm
# import cmocean

import matplotlib
# matplotlib.use('Agg')
from settings import setting as s

import matplotlib.animation as animation
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation


proj_name1 = 'cityair'
proj_name2 = 'cityair_future'
proj1 = s.projects[proj_name1]
proj2 = s.projects[proj_name2]


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


def plot_map_single(x, dom_name, pol_name, cbar_label='', cmap='twilight',
                    cb_limits=None, ax=None):
    # import matplotlib.pyplot as plt
    # import matplotlib.transforms as mtransforms
    ccrs_proj = ccrs.LambertConformal(
        central_longitude=24, central_latitude=45,
        false_easting=400000, false_northing=400000,
        standard_parallels=(46, 49))
    cbar_kws = dict(label=cbar_label, pad=0.02, shrink=0.8)
    res = '10m'
    grid_interval = 1
    if dom_name == 'tr':
        grid_interval = 2
    elif dom_name == 'eu':
        grid_interval = 3

    xticks = list(np.arange(-180, 180, grid_interval))
    yticks = list(np.arange(-90, 90, grid_interval))
    if ax is None:
        if cb_limits is None:
            p = x.plot(x='Longitude', y='Latitude',
                       robust=True,
                       cbar_kwargs=cbar_kws, cmap=cmap,
                       subplot_kws=dict(projection=ccrs_proj),
                       transform=ccrs.PlateCarree(), zorder=0,
                       alpha=1.0)
        else:
            vmin = min(cb_limits)
            vmax = max(cb_limits)
            p = x.plot(x='Longitude', y='Latitude',
                       robust=True,
                       cbar_kwargs=cbar_kws, cmap=cmap,
                       vmin=vmin, vmax=vmax,
                       subplot_kws=dict(projection=ccrs_proj),
                       transform=ccrs.PlateCarree(), zorder=0,
                       alpha=1.0)
    else:
        if cb_limits is None:
            p = x.plot(x='Longitude', y='Latitude',
                       robust=True, ax=ax,
                       cbar_kwargs=cbar_kws, cmap=cmap,
                       transform=ccrs.PlateCarree(), zorder=0,
                       alpha=1.0)
        else:
            vmin = min(cb_limits)
            vmax = max(cb_limits)
            p = x.plot(x='Longitude', y='Latitude',
                       robust=True, ax=ax,
                       cbar_kwargs=cbar_kws, cmap=cmap,
                       vmin=vmin, vmax=vmax,
                       transform=ccrs.PlateCarree(), zorder=0,
                       alpha=1.0)
    ax = p.axes
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.coastlines(resolution=res, alpha=0.6)
    p.figure.canvas.draw()
    ax.gridlines(xlocs=xticks, ylocs=yticks,
                 dms=True, color='indigo', alpha=0.5,
                 linestyle='--')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lambert_ticks(ax, xticks)
        lambert_ticks(ax, yticks, 'y')

    return p

    # _mkdir(path, exist_ok=True)
    # file_name = '_'.join(('plot', dom_name, pol_name, suffix)) + '.pdf'
    # p.fig.savefig(_join(path, file_name), bbox_inches='tight')
    # plt.close(p.fig)
    # print(file_name)


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


def get_data(dom_names, pol_names, months, proj, tstep=True):
    fmt_cro = 'GRIDCRO2D_{}_{}km_{}_{}{:02d}01.nc'
    fmt_com = 'COMBINE_ACONC_v{}_{}_{}_{}km_{}_{}{:02d}.nc'

    # check arg types
    if not isinstance(dom_names, list):
        dom_names = [dom_names]
    if not isinstance(pol_names, list):
        pol_names = [pol_names]
    if not isinstance(months, list):
        months = [months]

    for dn in dom_names:
        dom = proj.get_dom_by_name(dn)
        for m in months:
            month_name = calendar.month_name[m].lower()
            DIR_MCIP = _join(proj.path.mcip, '{}km'.format(dom.size), dom.name,
                             month_name + '_monthly')
            DIR_POST = proj.path.post

            CRO_FILE = fmt_cro.format(proj.name, dom.size, dom.name, year, m)
            COM_FILE = fmt_com.format(proj.cmaq_ver, proj.compiler, proj.name,
                                      dom.size, dom.name, year, m)

            lon, lat = get_latlon_from_cro(_join(DIR_MCIP, CRO_FILE))
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


def calc_stat(dom_names, pol_names, months, stat_day='mean', stat_mon='mean'):
    assert stat_day in ['mean', 'max', 'min'], \
        'stat_day arg must be mean, max or min.'
    assert stat_mon in ['mean', 'max', 'min'], \
        'stat_mon arg must be mean, max or min.'
    doms = {i: None for i in dom_names}
    for dn in dom_names:
        pols = {i: None for i in pol_names}
        for m in months:
            for pol_name in pol_names:
                for x in get_data(dn, pol_name, m, tstep=False):
                    xd = getattr(x.groupby('time.day'), stat_day)()
                    days = pd.to_datetime(x.coords['time'].values)
                    days = pd.date_range(min(days), max(days))
                    xd = xd.assign_coords(day=days)
                    xd.attrs['long_name'] = pol_name

                    xm = getattr(xd.groupby('day.month'), stat_mon)()

                    if pols[pol_name] is None:
                        pols[pol_name] = xm
                    else:
                        pols[pol_name] = xr.concat([pols[pol_name], xm],
                                                   dim='month')
        doms[dn] = xr.concat([i for i in pols.values()],
                             pd.Index(list(pols.keys()), name="pol_name"))
        doms[dn].attrs['units'] = 'ugm-3'
    return doms


POL_NAMES = ['NOX', 'O3', 'CO', 'SO2_UGM3', 'PM10', 'PM25_TOT']
DOM_NAMES = ['aegean', 'central_blacksea', 'mediterranean',
             'south_central_anatolia']
# DOM_NAMES = ['aegean']
POL_LABELS = dict(zip(POL_NAMES,
                      ['$NO_x$', '$O_3$', '$CO$', '$SO_2$', r'$PM_{10}$',
                       r'$PM_{2.5}$']))

POL_UNITS = dict(zip(POL_NAMES,
                     ['$(\\mu g/m^3)$', '$(\\mu g/m^3)$',
                      '$(\\mu g/m^3)$', '$(\\mu g/m^3)$', '$(\\mu g/m^3)$',
                      '$(\\mu g/m^3)$']))

year = 2015
months = [1, 2, 3]
# cmap = cmocean.cm.thermal_r
cmap = 'twilight'
NO_LIMIT = False

cmap_str = cmap if isinstance(cmap, str) else cmap.name
no_limit_str = 'no_limit' if NO_LIMIT else ''
pth = _join('plots_single_anom', '_'.join(('plots', cmap_str, no_limit_str)))

CB_LIMITS = None if NO_LIMIT else {'CO': [-50, 50],
                                   'NOX': [-50, 50],
                                   'O3': [-50, 50],
                                   'PM10': [-50, 50],
                                   'PM25_TOT': [-50, 50],
                                   'SO2_UGM3': [-50, 50]}

for dn in DOM_NAMES:
    for pn in POL_NAMES:
        file_name = '_'.join((proj_name1, 'vs', proj_name2, dn, pn)) + '.mp4'
        if not _isfile(file_name):
            print(file_name)
            cbar_label = ' '.join((POL_LABELS[pn], POL_UNITS[pn]))
            dat1 = get_data(dn, pn, months, proj1, tstep=True)
            dat2 = get_data(dn, pn, months, proj2, tstep=True)
            d1 = next(dat1)
            d2 = next(dat2)
            x = 100 * (d2 - d1) / d1
            cax = plot_map_single(x, dn, pn, cbar_label, cmap, CB_LIMITS[pn])

            fig = cax.figure
            ax = cax.axes

            def update(i):
                d1 = next(dat1)
                d2 = next(dat2)
                x = 100 * (d2 - d1) / d1
                plt.clf()
                cax = plot_map_single(x, dn, pn, cbar_label, cmap,
                                      CB_LIMITS[pn])
                # cax.set_array(x.values.flatten())
                # ax.set_title("Time = " + str(x.coords['time'].values))
                return cax

            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=10, metadata=dict(artist='isezen'),
                            bitrate=192000)
            anim = FuncAnimation(fig, update, frames=2060, interval=5)

            anim.save(file_name, writer=writer,
                      dpi=100, savefig_kwargs={'facecolor': '#ffffff'})

            plt.close(fig)
