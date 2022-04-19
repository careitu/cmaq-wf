#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101


"""
Plot stations and grid points
~~~~~~~~
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shapely.geometry as sgeom
import cartopy.feature as cfeature
from adjustText import adjust_text
from scipy import interpolate
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from DAL import GridData
from DALObs import StaData

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

def plot_station(locs, g, plot_name=None, ext=None, ng_method='loc'):
    geo = ccrs.Geodetic()
    # ccrs_proj = ccrs.PlateCarree()
    ccrs_proj = ccrs.LambertConformal(
        central_longitude=g.XCENT, central_latitude=g.YCENT,
        false_easting=400000, false_northing=400000,
        standard_parallels=(46, 49))

    # glons, glats = g.xlons, g.ylats
    glons, glats = g.lons, g.lats
    lons = np.array([l.lon for l in locs])
    lats = np.array([l.lat for l in locs])

    ng_method = 'loc' if ng_method == 'loc' else 'hav'
    ng_method = getattr(g, f'nearest_grid_{ng_method}')
    ngl = [ng_method(lat, lon) for lat, lon in zip(lats, lons)]
    # ngl = [g.nearest_grid_loc(loc) for loc in locs]
    grid_lons = np.array([l.lon for l in ngl])
    grid_lats = np.array([l.lat for l in ngl])

    points = ccrs_proj.transform_points(geo, lons, lats)
    # tlons, tlats = points[:, 0], points[:, 1]
    sta_names = [l.name for l in locs]
    title = plot_name if plot_name is not None else locs[0].region
    lonb = [min(lons), max(lons)]
    latb = [min(lats), max(lats)]
    
    res = '10m'
    grid_interval = 0.5
    xticks = list(np.arange(-180, 180, grid_interval))
    yticks = list(np.arange(-90, 90, grid_interval))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs_proj)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
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

    ax.scatter(glons, glats, marker='x', transform=geo)
    ax.scatter(lons, lats, edgecolors='black', transform=geo)

    for i in range(len(locs)):
        ax.plot([lons[i], grid_lons[i]], [lats[i], grid_lats[i]],
                transform=geo)

    delta = 2000
    if ext is None:
        ext = [lonb[0], lonb[1], latb[0], latb[1]]
    
    ext = np.array(ext)
    p = ccrs_proj.transform_points(geo, ext[0:2], ext[2:4])
    ext = [p[0, 0] - delta, p[1, 0] + delta,
           p[0, 1] - delta, p[1, 1] + delta]
    ax.set_extent(ext, crs=ccrs_proj)

    plt.title(title)

    texts = []
    for x, y, s in zip(points[:, 0], points[:, 1], sta_names):
        if ext[0] < x < ext[1] and ext[2] < y < ext[3]:
            texts.append(plt.text(x, y, s))
    adjust_text(texts, autoalign='y')
    plot_name = plot_name if plot_name is not None else title
    plt.savefig(f'plot_{plot_name.lower()}.pdf', bbox_inches='tight')
    plt.close(ax.figure)


sta = StaData()  # load Stations

# Aegean
g = GridData.from_dom('cityair', 'aegean')
locs = [v for k, v in sta.aegean.loc.items()]
plot_station(locs, g, 'aegean1', ext=[26.83, 27.25, 38.23, 38.55])
locs = [v for k, v in sta.aegean.loc.items() if not k.startswith('Izmir')]
plot_station(locs, g, 'aegean2', ext=[27.3, 27.7, 38.6, 39.2])
locs = [v for v in locs if not v.name.startswith('Manisa')]
plot_station(locs, g, 'aegean3', ext=[27.7, 28.6, 37.2, 37.8])

# South Central Anatolia
g = GridData.from_dom('cityair', 'south_central_anatolia')
locs = [v for k, v in sta.south_central_anatolia.loc.items()]
plot_station(locs, g, 'sca1', ext=[32.45, 32.55, 37.86, 37.93])
plot_station(locs, g, 'sca2', ext=[35.3, 35.55, 38.67, 38.76], ng_method='hav')

kay_melikgazi = locs[7]
lat = kay_melikgazi.lat
lon = kay_melikgazi.lon
loc = g.nearest_grid_loc(lat, lon)
hav = g.nearest_grid_hav(lat, lon)

# In [231]: ilon, ilat
# Out[231]: (133, 88)

# x, y = (968376.9539606982, -606138.3500973004)

# x, y = (1404376.078960698, -1564501.1625973005)

ilon, ilat = (floor((x - 2 * self.XORIG) / self.XCELL),
              floor((y - 2 * self.YORIG) / self.YCELL))