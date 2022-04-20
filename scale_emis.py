#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703


"""
Scale emission files
~~~~~~~~
Python script to scale emissions by a factor of a number
"""

import logging
import os
import sys
from pathlib import Path
from netCDF4 import Dataset

from os.path import join as _join
from os.path import abspath as _abspath
from os.path import isdir as _isdir
from timeit import default_timer as timer
import shutil

import calendar

from settings import setting as s

proj = s.get_active_proj()

# ----------------------------------
log = logging.getLogger('scaling')
log.setLevel(logging.DEBUG)
# create formatter
fmt_str = '%(asctime)s.%(msecs)03d - %(name)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
# create file handler
fh = logging.FileHandler(s.log.file)
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
log.addHandler(fh)
# ----------------------------------

_emis_vars_ = ['TFLAG', 'AACD', 'ACET', 'ALD2', 'ALDX', 'APIN', 'BENZ',
               'CH4', 'CO', 'ETH', 'ETHA', 'ETHY', 'ETOH', 'FACD', 'FORM',
               'IOLE', 'ISOP', 'IVOC', 'KET', 'MEOH', 'NAPH', 'NH3', 'NO',
               'NO2', 'NVOL', 'OLE', 'PAL', 'PAR', 'PCA', 'PCL', 'PEC',
               'PFE', 'PH2O', 'PK', 'PMC', 'PMG', 'PMN', 'PMOTHR', 'PNA',
               'PNCOM', 'PNH4', 'PNO3', 'POC', 'PRPA', 'PSI', 'PSO4',
               'PTI', 'SO2', 'SULF', 'TERP', 'TOL', 'UNR', 'XYLMN']


def scale_emis_file(filename, factor, force=False):
    nco = Dataset(filename, 'r+')
    log.info(f"Scaling Emission File '{filename}' by {factor}")
    variables = list(nco.variables.keys())
    if not force:
        for v in variables:
            if v not in _emis_vars_:
                msg = f"Variable '{v}' is not in emission file. "
                msg += f"{filename} may not be a emission file. "
                msg += "Execution stopped. Use '--force' option to force "
                msg += "scale this file."
                raise ValueError(msg)
    if 'TFLAG' in variables:
        variables.remove('TFLAG')  # except TFLAG
    for v in variables:
        nco[v][:] *= factor
        log.debug(f'{v} variable scaled by {factor}')
    nco.close()


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'Scale emission files by factor of a floating number\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -s 2\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-f', '--force', action='store_true', default=False,
                   help='Force scaling even if not emission file')
    p.add_argument('-o', '--overwrite', action='store_true', default=False,
                   help='Overwrite if file exist')
    p.add_argument('FACTOR', type=float, help='Scaling Factor')
    p.add_argument('SOURCE', help='Source Directory')
    p.add_argument('TARGET', help='Target Directory')
    return p.parse_args()


if __name__ == "__main__":
    from _helper_functions_ import ExitHelper
    from _helper_functions_ import ScriptError
    from _helper_functions_ import run_script_combine
    from _helper_functions_ import expandgrid
    from _helper_functions_ import get_days

    a = _parse_args_()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    factor, source, target = a.FACTOR, a.SOURCE, a.TARGET
    force, overwrite = a.force, a.overwrite
    # factor, source, target = 0.3, "/mnt/disk2/projects/cityair_future/emis/", "/mnt/disk2/projects/cityair_future_wamc/emis/"
    # force, overwrite = False, False

    if not _isdir(source):
        raise ValueError(f"{source} is not a directory")

    # if not _isdir(target):
    #     raise ValueError(f"{target} is not a directory")

    source, target = _abspath(source), _abspath(target)

    # log.info(f'Scaling Emission files at {a.SOURCE}')
    # Make sure FACTOR is a numeric value
    # print(factor)
    # print(source)
    # print(target)

    # f1 = list(Path(source).rglob("*.[nN][cC]"))[0]
    for f1 in Path(source).rglob("*.[nN][cC]"):
        f2 = Path(str(f1).replace(source, target))
        dir2, _ = os.path.split(f2)
        os.makedirs(dir2, exist_ok=True)
        if not overwrite and f2.exists():
            raise FileExistsError(f"{f2} is exist. Use '--overwrite' to overwrite.")
        shutil.copy2(f1, f2)
        scale_emis_file(f2, factor, force)
        # base = os.path.basename(f)
        print(dir2, f2)


