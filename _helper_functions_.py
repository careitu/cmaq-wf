#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703

"""
Create/get CMAQ-WF settings
~~~~~~~~
"""
import argparse as _ap
import os as _os
import signal
import sys as _sys

from argparse import RawTextHelpFormatter as _rtformatter
from os.path import join as _join

import itertools

from collections.abc import Iterable
from datetime import datetime as _dt
from datetime import timezone as _tz

__version__ = '0.0.1.dev'
__author__ = 'Ismail SEZEN'
__email__ = 'sezenismail@gmail.com'
__license__ = 'AGPL v3.0'
__year__ = '2021'


class ScriptError(Exception):
    """ Script Error """

    def __init__(self, message, detail=None):
        self.message = message
        self.detail = detail

    def __str__(self):
        return str(self.message)


class ExitHelper():
    """ Class to exit from script gracefully """

    def __init__(self):
        self._state = False
        signal.signal(signal.SIGINT, self.change_state)

    def change_state(self, a, b):
        """ Change statie of exit status """
        print("\nStopping...", a, b)
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        self._state = True

    @property
    def exit(self):
        """ Get exit state """
        return self._state


def max_mult(n):
    """ find maximum multipliers of an integer """
    n = int(n)
    import math
    n_org = n
    i = 2
    fc = 1
    dif_prev = 0
    dif = math.inf
    prev = None
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            fc *= i
            dif_prev = dif
            dif = abs(n - fc)
            if dif > dif_prev:
                return prev
            else:
                prev = (n, fc)
    if prev is None:
        return(max_mult(n_org - 1))
    return prev


def _create_argparser_(description, epilog):
    """ Create an argparser object """
    file_py = _os.path.basename(_sys.argv[0])
    p = _ap.ArgumentParser(description=description,
                           epilog=epilog.format(file_py), add_help=False,
                           formatter_class=_rtformatter)
    p.add_argument('--help', action='help', default=_ap.SUPPRESS,
                   help='Show this help message and exit.')
    p.add_argument('--version', help="Version", action="version",
                   version='{} {}\n{} (c) {} {}'.format(file_py, __version__,
                                                        __license__, __year__,
                                                        __author__))
    p.add_argument('-p', '--print', action='store_true', default=False,
                   help='Verbose output')
    return p


def expandgrid(*itrs):
    """ expand iterables """
    v = [x if isinstance(x, Iterable) else [x] for i, x in enumerate(itrs)]
    product = list(itertools.product(*v))
    x = list({'Var{}'.format(i + 1): [x[i] for x in product]
              for i in range(len(v))}.values())
    return list(map(tuple, zip(*x)))


def get_days(year, month, day=list(range(1, 32))):
    """ Return days in specific year and month """
    def date_is_ok(year, month, day):
        """ Check (year, month, day) is a correct day """
        try:
            date_str = '{}-{}-{}'.format(year, month, day)
            _dt.strptime(date_str, '%Y-%m-%d').replace(tzinfo=_tz.utc)
            return True
        except:  # noqa: E722
            pass
        return False
    days = []
    for i in expandgrid(year, month, day):
        if date_is_ok(i[0], i[1], i[2]):
            days.append(i[2])
    return days


def del_files(dir, files):
    import glob
    files = _join(dir, files)
    for f in glob.glob(files):
        _os.remove(f)


def is_in_file(file_name, search_string):
    import mmap
    with open(file_name, 'rb', 0) as f, mmap.mmap(
      f.fileno(), 0, access=mmap.ACCESS_READ) as s:
        return s.find(search_string) != -1


def run_script(script, dom, str_date, run_name='mcip',
               success_str=b'NORMAL TERMINATION'):
    from _helper_functions_ import ScriptError
    import tempfile as tf
    import subprocess as subp
    from settings import setting as s

    proj = s.get_active_proj()
    log_dir = _join(proj.path.logs, run_name)
    _os.makedirs(log_dir, exist_ok=True)

    tmp_name = next(tf._get_candidate_names())
    fmt = '{}_{}_{}_{}'.format(run_name, dom.size, dom.name, str_date)
    file_script = '{}_{}.csh'.format(fmt, tmp_name)
    tmp_dir = tf.gettempdir()
    _os.chdir(tmp_dir)
    file_script = _join(tmp_dir, file_script)
    file_log = _join(log_dir, '{}.log'.format(fmt))
    file_err = _join(log_dir, '{}.err'.format(fmt))

    with open(file_script, 'w') as f:
        f.write("#!/bin/csh -f\n")
        f.write(script)

    subp.call(['chmod', '+x', file_script])
    with open(file_err, 'w') as fe:
        with open(file_log, 'w') as fl:
            subp.call([file_script], stdout=fl, stderr=fe)

    _os.remove(file_script)

    if not is_in_file(file_log, success_str):
        msg = 'Error running {}'.format(run_name)
        raise ScriptError(msg, 'See: {}'.format(file_log))


def run_script_mcip(script, dom, str_date):
    run_script(script, dom, str_date)


def run_script_icon(script, dom, str_date):
    run_script(script, dom, str_date, 'icon',
               b'Program  ICON completed successfully')


def run_script_bcon(script, dom, str_date):
    run_script(script, dom, str_date, 'bcon', b'BCON completed successfully')


def run_script_cctm(script, dom, str_date):
    run_script(script, dom, str_date, 'cctm', b'CMAQ TIMING REPORT')


def run_script_combine(script, dom, str_date):
    run_script(script, dom, str_date, 'cctm', b'API Shut down')
