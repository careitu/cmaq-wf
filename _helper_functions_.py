#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import argparse as _ap
import os as _os
import sys as _sys
from argparse import RawTextHelpFormatter as _rtformatter

__version__ = '0.0.1.dev'
__author__ = 'Ismail SEZEN'
__email__ = 'sezenismail@gmail.com'
__license__ = 'AGPL v3.0'
__year__ = '2021'


def _create_argparser_(description, epilog):
    """ Create an argparser object """
    file_py = _os.path.basename(_sys.argv[0])
    p = _ap.ArgumentParser(description=description,
                           epilog=epilog.format(file_py),
                           formatter_class=_rtformatter)
    p.add_argument('-v', '--version', help="Version", action="version",
                   version='{} {}\n{} (c) {} {}'.format(file_py, __version__,
                                                        __license__, __year__,
                                                        __author__))
    p.add_argument('-p', '--print', action='store_true', default=False,
                   help='Verbose output')
    return p
