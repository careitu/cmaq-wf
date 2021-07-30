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
import numpy as np

reg = 'mediterranean'
pth = '/mnt/ssd2/APPS/programs/CreateEmis'
pth1 = os.path.join(pth, 'outputs/CMAQ_ReadyEmissions_city_air_{}.pickle')
pth2 = os.path.join(pth, 'outputs_old/CMAQ_ReadyEmissions_city_air_{}.pickle')
f1 = open(pth1.format(reg), 'rb')
f2 = open(pth2.format(reg), 'rb')

p1 = pickle.load(f1)
p2 = pickle.load(f2)

for i in range(0, 13):
    for j in range(0, 7):
        file_name = 'plot_{:02d}_{:02d}.png'.format(i, j)
        plt.subplot(1, 2, 1)
        plt.imshow(p1[:, :, i, j], interpolation='nearest')
        plt.title('New Pickle ({}, {})'.format(i, j))
        plt.subplot(1, 2, 2)
        plt.imshow(p2[:, :, i, j], interpolation='nearest')
        plt.title('Old Pickle ({}, {})'.format(i, j))
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
