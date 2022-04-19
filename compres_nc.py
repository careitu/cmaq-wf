#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
ACF plot for ensembles
~~~~~~~
Python script to create plots
"""

# from os import remove
# from os.path import join
# from os import makedirs as mkdir
import xarray as xr
import numpy as np
# import pandas as pd

encoding = {'zlib': True, 'shuffle': True, 'complevel': 5,
            'dtype': np.dtype('float32')}
FILE_NAME = 'CCTM_CONC_v532_gcc_cityair_future_wamd_2015_4km_20150102.nc'
ds = xr.open_dataset(FILE_NAME, cache=True)

for k, v in ds.variables.items():
    ds[k] = np.round(v, 1)
    ds[k].encoding = encoding

ds.to_netcdf('deneme.nc')
del ds

ds1 = xr.open_dataset(FILE_NAME, cache=True)
ds2 = xr.open_dataset('deneme.nc', cache=True)


np.min(ds1['PRES'] - ds2['PRES'])
# array(0.0078125)

np.max(np.abs(ds1['PRES'] - ds2['PRES']))
np.mean(np.abs(ds1['PRES'] - ds2['PRES']))

# ds1['PRES'][0, 0, 0, 2].values.tolist()
# 101863.1953125

# ds2['PRES'][0, 0, 0, 2].values.tolist()
# 101863.203125

# -0.0078125

# {'zlib': True,
#  'shuffle': True,
#  'complevel': 5,
#  'fletcher32': False,
#  'contiguous': False,
#  'chunksizes': (5, 7, 20, 22),
#  'source': '/home/isezen/cmaq-wf/CCTM_CONC_v532_gcc_cityair_future_wamd_2015_4km_20150102.nc',
#  'original_shape': (25, 34, 94, 103),
#  'dtype': dtype('float32')}
