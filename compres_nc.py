#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
ACF plot for ensembles
~~~~~~~
Python script to create plots
"""

import logging

import sys
from os import getcwd as _getcwd
from os import makedirs as _makedirs
from os import remove as _remove
from os.path import split as _split
from os.path import isdir as _isdir
from os.path import isfile as _isfile
from os.path import abspath as _abspath
from pathlib import Path
from filecmp import cmp as _cmp

import shutil
import numpy as np
import xarray as xr

from settings import setting as s

proj = s.get_active_proj()

# ----------------------------------
log = logging.getLogger('compress_nc')
log.setLevel(logging.DEBUG)
# create formatter
fmt_str = '%(asctime)s.%(msecs)03d - %(name)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
# create file handler
fh = logging.FileHandler(s.log.file)
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
log.addHandler(fh)


def common_start(sa, sb):
    """returns the longest common substring from the beginning of sa and sb"""
    def _iter():
        for a, b in zip(sa, sb):
            if a == b:
                yield a
            else:
                return

    return ''.join(_iter())


def _is_netcdf(f):
    """ Check file is a valied netcdf file """
    try:
        xr.open_dataset(f, engine='netcdf4')
    except OSError:
        return False
    return True


def _update_is_required(f, dec_prec=2, complevel=5):
    """ Check NetCDF file is compressed """
    dp_str = 'decimal_precision'
    compression_changed = False
    with xr.open_dataset(f, engine='netcdf4') as ds:
        for _, v in ds.variables.items():
            if v.encoding['complevel'] != complevel:
                compression_changed = True
                break
        if dp_str in ds.attrs.keys():
            if ds.attrs[dp_str] is not None and dec_prec is not None:
                if dec_prec >= ds.attrs[dp_str]:
                    dec_prec = None
        if dec_prec is not None:
            compression_changed = True
    return compression_changed


def _get_source_dir(source):
    """ Get source_dir from file names """
    source_dir = ''
    if len(source) > 1:
        source_dir = source[0]
        for i in source[1:]:
            source_dir = common_start(source_dir, i)
            if source_dir == '':
                break
        if not _isdir(source_dir):
            source_dir, _ = _split(source_dir)
    return source_dir


def _get_file_args(paths):
    """ Get file/path arguments """
    if len(paths) == 1:
        source, target = paths[0], paths[0]
    elif len(paths) == 2:
        source, target = paths[0], paths[1]
    else:
        target = paths.pop() if _isdir(paths[-1]) else _getcwd()
        source = paths

    if not isinstance(source, list):
        if _isdir(source):
            source = list(Path(source).rglob("*.*"))
        elif _isfile(source):
            source = [source]
        else:
            raise ValueError('Error in source files')

    source, target = [_abspath(i) for i in source], _abspath(target)
    source_dir = _get_source_dir(source)

    if _isfile(target):
        if len(source) > 1:
            raise ValueError('Target must be a directory')

    return source, source_dir, target


def compress_nc(from_file, to_file=None, dec_prec=None, complevel=5):
    """ Compress NetCDF File """
    dp_str = 'decimal_precision'
    pm_str = 'precision_MAE'
    if to_file is None:
        to_file = from_file
    encoding = {'zlib': True, 'shuffle': True, 'complevel': complevel,
                'dtype': np.dtype('float32')}

    modify = from_file == to_file
    if modify:
        to_file += '.tmp'

    glob_mae = 0
    with xr.open_dataset(from_file) as ds:
        log.debug(f'Compressing {from_file}')
        for k, v in ds.variables.items():
            if k != 'TFLAG':
                if dec_prec is not None:
                    v = v.astype('float128')
                    rounded = np.round(v, dec_prec).astype('float32')
                    ds[k] = rounded
                    mae = np.max(abs(v - rounded)).values.tolist()
                    mae = mae.astype('float32')
                    ds[k].attrs[pm_str] = mae
                    glob_mae = max(glob_mae, mae)
                    log.debug(f'Processing {k}: MAE: {mae}')
            ds[k].encoding = encoding

        if dec_prec is not None:
            ds.attrs.update({dp_str: dec_prec, pm_str: glob_mae})

    ds.to_netcdf(to_file)
    if modify:
        _remove(from_file)
        shutil.move(to_file, from_file)
        log.debug(f'tmp file renamed to {from_file}')
    log.info(f'Compressed: {from_file}')


def _parse_args_():
    from _helper_functions_ import _create_argparser_  # pylint: disable=C0415
    DESCRIPTION = 'Compress cmaq files\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -s 2\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-d', '--dlevel', type=int, required=False, default=5,
                   help="Compression/deflate level between 0-9")
    p.add_argument('-q', '--quantize', type=int, required=False, default=2,
                   help="Truncate data in variables to a given decimal \n\
precision, e.g. -q 2.")
    p.add_argument('-o', '--overwrite', action='store_true', default=False,
                   help='Overwrite if file exist')
    p.add_argument('PATHS', nargs='*', default=_getcwd(),
                   help='Files/Paths to compress from/to.')
    return p.parse_args()


if __name__ == "__main__":
    a = _parse_args_()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    dlevel, dec_precision, overwrite = a.dlevel, a.quantize, a.overwrite
    source, source_dir, target = _get_file_args(a.PATHS)

    for f1 in source:
        f2 = target if source_dir == '' else \
            Path(str(f1).replace(source_dir, target))
        dir2, _ = _split(f2)
        _makedirs(dir2, exist_ok=True)
        compress = False
        if _isfile(f1):
            if not _is_netcdf(f1):
                if _isfile(f2):
                    if _cmp(f1, f2):
                        msg = f"{f1} is identical to target."
                        log.debug(msg)
                else:
                    shutil.copy2(f1, f2)
                    msg = f"{f1} is not a valid netcdf file. Copied to target."
                    log.debug(msg)
                compress = False
                continue

        if _isfile(f2):
            if _is_netcdf(f2):
                compress = _update_is_required(f2, dec_precision, dlevel)
                if not compress:
                    log.info(f"{f2} is already compressed.")
                elif not overwrite:
                    raise FileExistsError('Use --overwrite to overwrite files')
            else:
                shutil.copy2(f1, f2)
                msg = f"{f1} is not a valid netcdf file. Copied to target."
                log.debug(msg)
        else:
            compress = True

        if compress:
            compress_nc(f1, f2, dec_prec=dec_precision, complevel=dlevel)


# f2 = 'tmp/CCTM_CONC_v532_gcc_cityair_future_wamd_2015_4km_20150102_yedek.nc'
# compress_nc(FILE_NAME, 'deneme.nc', 2)

# ds1 = xr.open_dataset(FILE_NAME, cache=True)
# ds2 = xr.open_dataset('deneme.nc', cache=True)


# dif = ds1 - ds2

# np.min(dif)

# np.min(ds1['PRES'] - ds2['PRES'])
# # array(0.0078125)

# np.max(np.abs(ds1['PRES'] - ds2['PRES']))
# np.mean(np.abs(ds1['PRES'] - ds2['PRES']))

# {'zlib': True,
#  'shuffle': True,
#  'complevel': 5,
#  'fletcher32': False,
#  'contiguous': False,
#  'chunksizes': (5, 7, 20, 22),
#  'source': '/home/isezen/cmaq-wf/CCTM_CONC_v532_gcc_cityair_future_wamd_2015\
#             _4km_20150102.nc',
#  'original_shape': (25, 34, 94, 103),
#  'dtype': dtype('float32')}
