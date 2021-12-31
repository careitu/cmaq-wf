#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101


"""
Plot station coordinates
~~~~~~~~
"""

import warnings
import numpy as np
from DAL import GridData
from DALObs import StaData
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shapely.geometry as sgeom
import cartopy.feature as cfeature
from adjustText import adjust_text
from scipy import interpolate
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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
            print(ls.bounds)
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

def plot_station(locs, g):
    lons = np.array([l.lon for l in locs])
    lats = np.array([l.lat for l in locs])
    sta_names = [l.name for l in locs]
    title = locs[0].region
    lonb = [min(lons), max(lons)]
    latb = [min(lats), max(lats)]

    geo = ccrs.Geodetic()
    ccrs_proj = ccrs.LambertConformal(
        central_longitude=g.XCENT, central_latitude=g.YCENT,
        false_easting=400000, false_northing=400000,
        standard_parallels=(46, 49))
    # ccrs_proj = ccrs.PlateCarree()
    points = ccrs_proj.transform_points(geo, lons, lats)

    res = '10m'
    grid_interval = 0.5
    xticks = list(np.arange(-180, 180, grid_interval))
    yticks = list(np.arange(-90, 90, grid_interval))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs_proj)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.scatter(lons, lats, edgecolors='b', transform=geo)
    ax.coastlines(resolution=res, alpha=0.7)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAKES, alpha=0.7)
    gl = ax.gridlines(xlocs=xticks, ylocs=yticks,
                      dms=True, color='indigo', alpha=0.5,
                      linestyle='--')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    # ax.set_extent([lonb[0], lonb[1],
    #                latb[0], latb[1]])
    plt.title(title)
    ax.set_extent([g.bounds.lon[0], g.bounds.lon[1],
                   g.bounds.lat[0], g.bounds.lat[1]])
    texts = []
    for x, y, s in zip(points[:, 0], points[:, 1], sta_names):
        texts.append(plt.text(x, y, s))
    f = interpolate.interp1d(points[:, 0], points[:, 1])
    x = np.arange(min(points[:, 0]), max(points[:, 0]), 100)
    y = f(x)
    adjust_text(texts, x=x, y=y, autoalign='y',
                only_move={'points':'y', 'texts':'y'}, force_points=100,
                arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    plt.savefig(f'plot_{title}.pdf', bbox_inches='tight')
    plt.close(ax.figure)


sta = StaData()  # load Stations

for r in sta:
    locs = [v for _, v in r.loc.items()]
    g = GridData.from_dom('cityair', locs[0].region)
    plot_station(locs, g)
