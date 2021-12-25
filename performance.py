#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=E1101


"""
Calculşate performance for CMAQ vs. Observations
~~~~~~~~
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap as lscmap

from settings import setting as _set

matplotlib.use('Agg')

from settings import setting as _set
from DAL import *

proj = _set.get_proj_by_name('cityair')
dom_names=['aegean', 'mediterranean',
           'south_central_anatolia',
           'central_blacksea']
pol_names = list(proj.pols.keys())
pols = {p: proj.pols[p] for p in pol_names}
pol_names = [p.name for p in pols.values()]

# proj_names='cityair',
pp = PostProc(dom_names=dom_names)

d = pp.domains.aegean.get_data(slice(0, 640))

lats, lons = [36.2, 38.4217, 39.0, 37.0], [25.8, 27.1633, 28.0, 29.0]
# lats, lons = [38.4217, 39.0, 37.0], [27.1633, 28.0, 29.0]
locs = [Location(lat, lon) for (lat, lon) in zip(lats, lons)]

d2 = pp.domains.aegean.get_data_loc(slice('2015-01-05', '2015-01-10'),
                                    loc=locs, delta=1, layer_mean=True,
                                    simplify=False)

d2 = pp.domains.aegean.get_data_loc(loc=locs, delta=1, layer_mean=True,
                                    simplify=False)


d3 = pp.domains.aegean.get_data_loc(slice('2015-01-05', '2015-01-10'),
                                    locs, delta=0, layer_mean=False)
