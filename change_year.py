#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Change year in CMAQ netcdf files
~~~~~~~
Python script to change year in CMAQ netcdf files
"""
import os
from netCDF4 import Dataset
from settings import setting as s

proj = s.get_active_proj()


def getListOfFiles(dirName):
    listOfFile = os.listdir(dirName)
    for entry in listOfFile:
        fullPath = os.path.join(dirName, entry)
        if os.path.isdir(fullPath):
            yield from getListOfFiles(fullPath)
        else:
            yield fullPath


def change_year(nc_file, new_year, var_name='TFLAG'):
    nco = Dataset(nc_file, 'r+')
    v = nco.variables[var_name][:]
    s = v.shape
    for i in range(s[0]):
        for j in range(s[1]):
            val = v[i][j][0]
            # nco[var_name][i, j, 0] = val + (new_year - val // 1000) * 1000
            nco[var_name][i, j, 0] = new_year * 1000 + val % 1000
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
    p.add_argument('-d', '--dir', type=str, default=os.getcwd(),
                   help='Directory to search. Default is current directory')
    p.add_argument('file_pattern', metavar='FILE_PATTERN', type=str, nargs='+',
                   help='File patterns to search')
    p.add_argument('-y', '--year', type=int, default=2025,
                   help='default is 2025.')
    return p.parse_args()


var_name = 'TFLAG'
nc_file = 'METCRO2D_cityair_4km_aegean_20150102.nc'
new_year = 2025


if __name__ == "__main__":
    a = _parse_args_()
    if a.dir is None:
        d = os.getcwd()
    else:
        if a.dir in proj.path.__dict__.keys():
            d = proj.path.__dict__[a.dir]
        else:
            d = a.dir

    print('\nDirectory is "{}"'.format(d))
    user_input = input('Confirm? [Y/N] ')

    if user_input.lower() in ('y', 'yes'):
        for f in getListOfFiles(d):
            for fp in a.file_pattern:
                if fp.lower() in f.lower():
                    print(f)



    # for f in a.files:
    #     if f in proj.path.__dict__.keys():
    #         f = proj.path.__dict__[f]
    #         files = getListOfFiles(f)
    #         for k in files:
    #             print(k)
    #     else:
    #         print(f)
        # change_year(f, a.year)
