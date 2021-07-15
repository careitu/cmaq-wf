#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Change year in CMAQ netcdf files
~~~~~~~
Python script to change year in CMAQ netcdf files
"""

from netCDF4 import Dataset
from settings import setting as s

proj = s.get_active_proj()


def change_year(nc_file, new_year, var_name='TFLAG'):
    nco = Dataset(nc_file, 'r+')
    v = nco.variables[var_name][:]
    s = v.shape
    for i in range(s[0]):
        for j in range(s[1]):
            val = v[i][j][0]
            nco[var_name][i, j, 0] = val + (new_year - val // 1000) * 1000
    nco.close()


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'Change year in netcdf file\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -y 2025 file_name.nc\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('files', metavar='FILES', type=str, nargs='+',
                   help='File names')
    p.add_argument('-y', '--year', type=int, default=2025,
                   help='default is 2025.')
    return p.parse_args()


var_name = 'TFLAG'
nc_file = 'METCRO2D_cityair_4km_aegean_20150102.nc'
new_year = 2025


if __name__ == "__main__":
    a = _parse_args_()

    for f in a.files:
        change_year(f, a.year)
