#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703


"""
Run mcip
~~~~~~~~
Python script to run mcip
"""

import os
import sys
import itertools
import logging
from os.path import join
from timeit import default_timer as timer

import calendar
from collections.abc import Iterable
from datetime import datetime as _dt
from datetime import timedelta as _td
from datetime import timezone as _tz

from settings import setting as s

proj = s.get_active_proj()

# ----------------------------------
log = logging.getLogger('{}.icon'.format(proj.name))
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


def get_script(year, month, day, dom, type='auto'):
    type = 'regrid' if type == 'auto' and dom.__parent__ is not None \
        else 'profile'
    dom_parent_size = None if dom.__parent__ is None else dom.__parent__.size
    mn = calendar.month_name[month].lower()
    script = """
setenv compiler {}
source /mnt/ssd2/APPS/CMAQ/config_cmaq.csh ${{compiler}}

if ( ! -e $CMAQ_DATA ) then
    echo "   $CMAQ_DATA path does not exist"
    exit 1
endif
echo " "; echo " Input data path, CMAQ_DATA set to $CMAQ_DATA"; echo " "

set year = {:04}
set month = {:02d}
set month_name = {}
set day = {:02d}
set dom_size = {}km
set dom_size_par = {}km
set proj_name = {}
set dom_name = {}

set cmaq_home = {}
set OUTDIR = {}
set dir_mcip = {}
set parent_dom_cctm_path = {}/${{dom_size_par}}

set VRSN     = v{}
set ICTYPE   = {}

set dir_mcip1 = ${{dir_mcip}}/${{dom_size_par}}/${{dom_name}}/${{month_name}}
set dir_mcip2 = ${{dir_mcip}}/${{dom_size}}/${{dom_name}}/${{month_name}}

set path_bld = ${{cmaq_home}}/PREP/icon/scripts
set BLD      = ${{path_bld}}/BLD_ICON_${{VRSN}}_${{compiler}}
set EXEC     = ICON_${{VRSN}}.exe
cat $BLD/ICON_${{VRSN}}.cfg; echo " "; set echo

setenv GRID_NAME ${{dom_size}}
setenv GRIDDESC  ${{dir_mcip2}}/GRIDDESC
setenv IOAPI_ISPH 20

setenv IOAPI_LOG_WRITE F
setenv IOAPI_OFFSET_64 YES
setenv EXECUTION_ID ${{EXEC}}

setenv ICON_TYPE ` echo $ICTYPE | tr "[A-Z]" "[a-z]" `

set ymd  = ${{year}}${{month}}${{day}}
set APPL = ${{proj_name}}_${{dom_size}}_${{dom_name}}_${{ymd}}

setenv INIT_CONC_1 "${{OUTDIR}}/ICON_${{VRSN}}_${{ICON_TYPE}}_${{APPL}} -v"
setenv MET_CRO_3D_FIN ${{dir_mcip2}}/METCRO3D_${{APPL}}.nc

if ( $ICON_TYPE == regrid ) then
  set cctm_file = CCTM_CONC_${{VRSN}}_${{compiler}}_${{proj_name}}_${{year}}
  set cctm_file = ${{cctm_file}}_${{dom_size_par}}_${{ymd}}.nc
  setenv CTM_CONC_1 ${{parent_dom_cctm_path}}/${{cctm_file}}
endif

if ( $ICON_TYPE == profile ) then
  set av = avprofile_cb6r3m_ae7_kmtbr_hemi2016_v53beta2_m3dry_col051_row068.csv
  setenv IC_PROFILE $BLD/profiles/$av
endif

if ( ! -d "${{OUTDIR}}" ) mkdir -p ${{OUTDIR}}

ls -l $BLD/$EXEC; size $BLD/$EXEC
unlimit
limit

time $BLD/$EXEC

exit()""".format(proj.compiler, year, month, mn, day, dom.size,
                 dom_parent_size, proj.name, dom.name, proj.path.cmaq_app,
                 proj.path.icon, proj.path.mcip, proj.path.cctm,
                 proj.cmaq_ver, type)
    return script


def expandgrid(*itrs):
    """ expand iterables """
    v = [x if isinstance(x, Iterable) else [x] for i, x in enumerate(itrs)]
    product = list(itertools.product(*v))
    x = list({'Var{}'.format(i + 1): [x[i] for x in product]
              for i in range(len(v))}.values())
    return list(map(tuple, zip(*x)))


def date_is_ok(year, month, day):
    """ Check (year, month, day) is a correct day """
    try:
        date_str = '{}-{}-{}'.format(year, month, day)
        _dt.strptime(date_str, '%Y-%m-%d').replace(tzinfo=_tz.utc)
        return True
    except:  # noqa: E722
        pass
    return False


def get_days(year, month, day=list(range(1, 32))):
    """ Return days in specific year and month """
    days = []
    for i in expandgrid(year, month, day):
        if date_is_ok(i[0], i[1], i[2]):
            days.append(i[2])
    return days


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'icon script\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -d 4 5 6\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-n', '--domain', nargs='+', type=int,
                   default=proj.get_dom_ids(),
                   help="domain Id(s). Default is all domains in the project.")
    p.add_argument('-t', '--type', default='auto',
                   choices=['auto', 'profile', 'regrid'],
                   help="default is 'auto'.")
    p.add_argument('-y', '--years', nargs='+', type=int, default=proj.years,
                   help='default is all years in config file.')
    p.add_argument('-m', '--months', nargs='+', type=int, default=proj.months,
                   help='default is all months in config file.')
    p.add_argument('-d', '--days', nargs='+', type=int, default=1,
                   help='default is 1.')
    return p.parse_args()


if __name__ == "__main__":
    from _helper_functions_ import ExitHelper
    from _helper_functions_ import ScriptError
    from _helper_functions_ import run_script_icon
    a = _parse_args_()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    year, month, day = a.years, a.months, a.days
    doms = [proj.get_dom_by_id(i) for i in a.domain]
    ym = expandgrid(year, month)  # Year and months

    log_dir = os.path.join(proj.path.logs, 'icon')
    os.makedirs(log_dir, exist_ok=True)

    log.info('Creating icon files')
    log.info('log dir: {}'.format(log_dir))

    # make sure output dir exist
    dir_out = join(proj.path.proj, 'icon')
    os.makedirs(dir_out, exist_ok=True)
    log.info('Output Dir: {}'.format(dir_out))

    time_fmt = '[Time: {:.2f} secs]'
    flag = ExitHelper()
    icon_timer_start = timer()
    for dom in doms:
        dom_timer_start = timer()
        dom_str = 'Dom: {}-{}'.format(dom.size, dom.name)
        log.info(dom_str + ' starting...')

        for y, m in ym:
            month_str = 'Month: {}-{:02d}'.format(y, m)
            log.info('Processing Month {}-{:02d}'.format(y, m))
            days = get_days(y, m, day)

            month_timer_start = timer()
            for d in days:
                if flag.exit:
                    log.info('User stopped execution')
                    sys.exit()
                str_date = '{:04d}-{:02d}-{:02d}'.format(y, m, d)
                day_str = 'Day: {}'.format(str_date)
                script = get_script(y, m, d, dom, a.type)

                day_timer_start = timer()
                err = None
                try:
                    run_script_icon(script, dom, str_date)
                except ScriptError as error:
                    err = error
                day_timer_end = timer()

                el = time_fmt.format(day_timer_end - day_timer_start)
                msg = dom_str + ', ' + day_str
                if err is None:
                    log.info(msg + ' ' + el)
                else:
                    log.error(msg + ', ' + err.detail + ' ' + el)
                    err = None

            el = time_fmt.format(day_timer_end - month_timer_start)
            msg = dom_str + ', ' + month_str
            log.info(msg + ' ' + el)

        msg = '{} completed '.format(dom_str)
        el = day_timer_end - dom_timer_start
        log.info(msg + time_fmt.format(el))

    el = day_timer_end - icon_timer_start
    log.info('icon completed ' + time_fmt.format(el))
