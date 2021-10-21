#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Compare pickles
~~~~~~~
"""


import os
import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# aegean central_blacksea mediterranean south_central_anatolia
reg = 'south_central_anatolia'
pth = '/mnt/ssd2/APPS/programs/CreateEmis'
pth1 = os.path.join(pth, 'outputs/CMAQ_ReadyEmissions_city_air_{}.pickle')
pth2 = os.path.join(pth, 'outputs_old/CMAQ_ReadyEmissions_city_air_{}.pickle')
f1 = open(pth1.format(reg), 'rb')
f2 = open(pth2.format(reg), 'rb')

p1 = pickle.load(f1)
p2 = pickle.load(f2)

for i in range(0, 13):
    for j in range(0, 7):
        d = p1[:, :, i, j] - p2[:, :, i, j]
        file_name = 'plot_{}_{:02d}_{:02d}.png'.format(reg, i, j)
        ax = plt.subplot()
        plt.title('{} (New - Old) Pickle [{}, {}]'.format(reg, i, j))
        im = ax.imshow(d, interpolation='nearest')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
