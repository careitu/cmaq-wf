#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703


"""
Run post combine
~~~~~~~~
Python script to run post combine
"""

import logging
import os
import sys

from os.path import join as _join
from timeit import default_timer as timer

import calendar

from settings import setting as s

proj = s.get_active_proj()

# ----------------------------------
log = logging.getLogger('{}.post'.format(proj.name))
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


def get_script(year, month, start_day, end_day, dom):
    fmt = '{:04d}-{:02d}-{:02d}'
    START_DATE = fmt.format(year, month, start_day)
    END_DATE = fmt.format(year, month, end_day)
    mn = calendar.month_name[month].lower()
    script = """
setenv compiler {}

source /mnt/ssd2/APPS/CMAQ/config_cmaq.csh ${{compiler}}

set year = {}
set month_name = {}
set proj_name = {}
set dom_name = {}
set dom_size = {}km
set mcip_dir = {}
set cctm_dir = {}
set post_dir = {}

set VRSN = v{}
set MECH = cb6r3_ae7_aq
set APPL = ${{proj_name}}_${{year}}_${{dom_size}}

setenv RUNID  ${{VRSN}}_${{compilerString}}_${{APPL}}

if ( ! $?BINDIR ) then
  set bin_path = ${{CMAQ_HOME}}/POST/combine/scripts
  setenv BINDIR ${{bin_path}}/BLD_combine_${{VRSN}}_${{compilerString}}
endif

setenv EXEC combine_${{VRSN}}.exe

setenv REPO_HOME  ${{CMAQ_REPO}}

setenv METDIR     ${{mcip_dir}}/${{dom_size}}/${{dom_name}}/${{month_name}}_monthly
setenv CCTMOUTDIR ${{cctm_dir}}/${{dom_size}}/${{dom_name}}
setenv POSTDIR    ${{post_dir}}

if ( ! -e $POSTDIR ) then
  mkdir $POSTDIR
endif

set START_DATE = "{}"
set END_DATE   = "{}"

set spec_def_path = $REPO_HOME/POST/combine/scripts/spec_def_files
setenv SPEC_CONC $spec_def_path/SpecDef_${{MECH}}.txt
setenv SPEC_DEP  $spec_def_path/SpecDef_Dep_${{MECH}}.txt

setenv GENSPEC N

setenv SPECIES_DEF $SPEC_CONC

set TODAYG = ${{START_DATE}}
set TODAYJ = `date -ud "${{START_DATE}}" +%Y%j`
set STOP_DAY = `date -ud "${{END_DATE}}" +%Y%j`

while ($TODAYJ <= $STOP_DAY )

  set YYYY = `date -ud "${{TODAYG}}" +%Y`
  set YY = `date -ud "${{TODAYG}}" +%y`
  set MM = `date -ud "${{TODAYG}}" +%m`
  set DD = `date -ud "${{TODAYG}}" +%d`
  setenv OUTFILE ${{POSTDIR}}/COMBINE_ACONC_${{RUNID}}_$YYYY$MM.nc

  appl_mcip = ${{proj_name}}_${{dom_size}}_${{dom_name}}_$YYYY$MM01.nc
  setenv INFILE1 $CCTMOUTDIR/CCTM_ACONC_${{RUNID}}_$YYYY$MM$DD.nc
  setenv INFILE2 $METDIR/METCRO3D_${{appl_mcip}}
  setenv INFILE3 $CCTMOUTDIR/CCTM_APMDIAG_${{RUNID}}_$YYYY$MM$DD.nc
  setenv INFILE4 $METDIR/METCRO2D_${{appl_mcip}}

  ${{BINDIR}}/${{EXEC}}

  set TODAYG = `date -ud "${{TODAYG}}+1days" +%Y-%m-%d`
  set TODAYJ = `date -ud "${{TODAYG}}" +%Y%j`

end

setenv SPECIES_DEF $SPEC_DEP

set TODAYG = ${{START_DATE}}
set TODAYJ = `date -ud "${{START_DATE}}" +%Y%j`
set STOP_DAY = `date -ud "${{END_DATE}}" +%Y%j`

while ($TODAYJ <= $STOP_DAY )

  set YYYY = `date -ud "${{TODAYG}}" +%Y`
  set YY = `date -ud "${{TODAYG}}" +%y`
  set MM = `date -ud "${{TODAYG}}" +%m`
  set DD = `date -ud "${{TODAYG}}" +%d`
  setenv OUTFILE ${{POSTDIR}}/COMBINE_DEP_${{RUNID}}_$YYYY$MM

  appl_mcip = ${{proj_name}}_${{dom_size}}_${{dom_name}}_$YYYY$MM01.nc
  setenv INFILE1 $CCTMOUTDIR/CCTM_DRYDEP_${{RUNID}}_$YYYY$MM$DD.nc
  setenv INFILE2 $CCTMOUTDIR/CCTM_WETDEP1_${{RUNID}}_$YYYY$MM$DD.nc
  setenv INFILE3 $METDIR/METCRO2D_${{appl_mcip}}
  setenv INFILE4

  ${{BINDIR}}/${{EXEC}}

  set TODAYG = `date -ud "${{TODAYG}}+1days" +%Y-%m-%d`
  set TODAYJ = `date -ud "${{TODAYG}}" +%Y%j`

end

exit()""".format(proj.compiler, year, mn, proj.name, dom.name, dom.size,
                 proj.path.mcip, proj.path.cctm, proj.path.post,
                 proj.cmaq_ver, START_DATE, END_DATE)
    return script


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'Combine post script\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -s 2\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-n', '--domain', nargs='+', type=int,
                   default=proj.get_dom_ids(),
                   help="domain Id(s). Default is all domains in the project.")
    p.add_argument('-y', '--years', nargs='+', type=int, default=proj.years,
                   help='default is all years in config file.')
    p.add_argument('-m', '--months', nargs='+', type=int, default=proj.months,
                   help='default is all months in config file.')
    p.add_argument('-s', '--start_day', type=int, default=1,
                   help='default is 1.')
    p.add_argument('-e', '--end_day', type=int, default=31,
                   help='default is 31.')
    return p.parse_args()


if __name__ == "__main__":
    from _helper_functions_ import ExitHelper
    from _helper_functions_ import ScriptError
    from _helper_functions_ import run_script_cctm
    from _helper_functions_ import expandgrid
    from _helper_functions_ import get_days

    a = _parse_args_()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    year, month, start_day, end_day = a.years, a.months, a.start_day, a.end_day
    day = [start_day] + [end_day + 1]
    doms = [proj.get_dom_by_id(i) for i in a.domain]
    ym = expandgrid(year, month)  # Year and months

    log_dir = _join(proj.path.logs, 'post')
    os.makedirs(log_dir, exist_ok=True)

    log.info('Creating combine files')
    log.info('log dir: {}'.format(log_dir))

    # make sure output dir exist
    dir_out = proj.path.post
    os.makedirs(dir_out, exist_ok=True)
    log.info('Output Dir: {}'.format(dir_out))

    time_fmt = '[Time: {:.2f} secs]'
    flag = ExitHelper()
    post_timer_start = timer()
    for dom in doms:
        dom_timer_start = timer()
        dom_str = 'Dom: {}-{}'.format(dom.size, dom.name)
        log.info(dom_str + ' starting...')

        for y, m in ym:
            if flag.exit:
                log.info('User stopped execution')
                sys.exit()
            month_str = 'Month: {}-{:02d}'.format(y, m)
            log.info('Processing Month {}-{:02d}'.format(y, m))
            days = get_days(y, m, day=list(range(day[0], day[1])))
            start_day = days[0]
            end_day = days[len(days) - 1]
            str_month = '{:04d}-{:02d}'.format(y, m)

            script = get_script(y, m, start_day, end_day, dom)

            month_timer_start = timer()

            err = None
            try:
                run_script_cctm(script, dom, str_month)
            except ScriptError as error:
                err = error
            month_timer_end = timer()

            el = time_fmt.format(month_timer_end - month_timer_start)
            msg = dom_str + ', ' + month_str
            if err is None:
                log.info(msg + ' ' + el)
            else:
                log.error(msg + ', ' + err.detail + ' ' + el)
                err = None

        msg = '{} completed '.format(dom_str)
        el = month_timer_end - dom_timer_start
        log.info(msg + time_fmt.format(el))

    el = month_timer_end - post_timer_start
    log.info('post combine completed ' + time_fmt.format(el))
