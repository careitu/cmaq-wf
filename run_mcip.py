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
from os.path import join as _join
from timeit import default_timer as timer

import calendar
from collections.abc import Iterable
from datetime import datetime as _dt
from datetime import timedelta as _td
from datetime import timezone as _tz

from settings import setting as s

proj = s.get_active_proj()
logging_dir = _join(os.path.expanduser("~"), '.config/cwf')
os.makedirs(logging_dir, exist_ok=True)

# ----------------------------------
log = logging.getLogger('{}.mcip'.format(proj.name))
log.setLevel(logging.DEBUG)
# create formatter
fmt_str = '%(asctime)s.%(msecs)03d - %(name)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
# create file handler
fh = logging.FileHandler(_join(logging_dir, 'cwf.log'))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
log.addHandler(fh)
# ----------------------------------

dir_prog = _join(proj.path.cmaq_app, 'PREP', 'mcip', 'src')
wrfout_fmt = '${{InMetDir}}' \
             '/wrfout_${{dom_num}}_${{year}}-${{month}}-{:02d}_00:00:00'
dir_in_geo = proj.path.wps
dir_out_fmt = _join(proj.path.mcip, '{}km', '{}', '{}')
dir_in_met_fmt = _join(proj.path.wrf, '{}')


def get_script(year, month, day, dom, proj_name, dir_in_met, dir_in_geo,
               dir_out, dir_prog, in_met_files, monthly=False, compiler='gcc'):
    fmt = '%Y-%m-%d-%H:%M:%S'
    if monthly:
        day = 1
        next_day = calendar.monthrange(year, month)[1]
        date_str = '{}-{}-{}-00:00:00'.format(year, month, day)
        day_after = _dt.strptime(date_str, fmt).replace(tzinfo=_tz.utc)
        day_after = day_after + _td(days=next_day) - _td(hours=1)
        day_after = day_after.strftime(fmt)
        mcip_start = '{:04d}-{:02d}-{:02d}-01:00:00'.format(year, month, day)
        mcip_end = day_after
    else:
        next_day = 1
        date_str = '{}-{}-{}-00:00:00'.format(year, month, day)
        day_after = _dt.strptime(date_str, fmt).replace(tzinfo=_tz.utc)
        day_after = day_after + _td(days=1)
        day_after = day_after.strftime(fmt)
        mcip_start = '{:04d}-{:02d}-{:02d}-00:00:00'.format(year, month, day)
        mcip_end = day_after
    script = """
source /mnt/ssd2/APPS/CMAQ/config_cmaq.csh {}

set year = {}
set month = {:02d}
set day = {:02d}
set dom_size = {}km
set dom_name = {}
set dom_num = d{:02d}
set project_name = {}

set ymd        = ${{year}}${{month}}${{day}}
set APPL       = ${{project_name}}_${{dom_size}}_${{dom_name}}_${{ymd}}
set CoordName  = LambertConformal
set GridName   = ${{dom_size}}

set DataPath   = $CMAQ_DATA
set InMetDir   = {}
set InGeoDir   = {}

set OutDir     = {}
set ProgDir    = {}
set WorkDir    = ${{OutDir}}

{}

set IfGeo      = "T"
set InGeoFile  = ${{InGeoDir}}/geo_em.${{dom_num}}.nc

set LPV     = 0
set LWOUT   = 0
set LUVBOUT = 1

set MCIP_START={}.0000
set MCIP_END={}.0000

set INTVL      = 60

set IOFORM = 1

set BTRIM = 0

set X0    =  1
set Y0    =  1
set NCOLS =  {}
set NROWS =  {}

set LPRT_COL = 0
set LPRT_ROW = 0

set WRF_LC_REF_LAT = -999.0

set PROG = mcip

date

if ( ! -d $InMetDir ) then
  echo "No such input directory $InMetDir"
  exit 1
endif

if ( ! -d $OutDir ) then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if ( $status != 0 ) then
    echo "Failed to make output directory, $OutDir"
    exit 1
  endif
endif

if ( ! -d $ProgDir ) then
  echo "No such program directory $ProgDir"
  exit 1
endif

if ( $IfGeo == "T" ) then
  if ( ! -f $InGeoFile ) then
    echo "No such input file $InGeoFile"
    exit 1
  endif
endif

foreach fil ( $InMetFiles )
  if ( ! -f $fil ) then
    echo "No such input file $fil"
    exit 1
  endif
end

if ( ! -f $ProgDir/$PROG.exe ) then
  echo "Could not find $PROG.exe"
  exit 1
endif

if ( ! -d $WorkDir ) then
  mkdir -p $WorkDir
  if ( $status != 0 ) then
    echo "Failed to make work directory, $WorkDir"
    exit 1
  endif
endif

cd $WorkDir

if ( $IfGeo == "T" ) then
  if ( -f $InGeoFile ) then
    set InGeo = $InGeoFile
  else
    set InGeo = "no_file"
  endif
else
  set InGeo = "no_file"
endif

set FILE_GD  = $OutDir/GRIDDESC

set MACHTYPE = `uname`
if ( ( $MACHTYPE == "AIX" ) || ( $MACHTYPE == "Darwin" ) ) then
  set Marker = "/"
else
  set Marker = "&END"
endif

cat > $WorkDir/namelist.$PROG << !

 &FILENAMES
  file_gd    = "$FILE_GD"
  file_mm    = "$InMetFiles[1]",
!

if ( $#InMetFiles > 1 ) then
  @ nn = 2
  while ( $nn <= $#InMetFiles )
    cat >> $WorkDir/namelist.$PROG << !
               "$InMetFiles[$nn]",
!
    @ nn ++
  end
endif

if ( $IfGeo == "T" ) then
cat >> $WorkDir/namelist.$PROG << !
  file_geo   = "$InGeo"
!
endif

cat >> $WorkDir/namelist.$PROG << !
  ioform     =  $IOFORM
 $Marker

 &USERDEFS
  lpv        =  $LPV
  lwout      =  $LWOUT
  luvbout    =  $LUVBOUT
  mcip_start = "$MCIP_START"
  mcip_end   = "$MCIP_END"
  intvl      =  $INTVL
  coordnam   = "$CoordName"
  grdnam     = "$GridName"
  btrim      =  $BTRIM
  lprt_col   =  $LPRT_COL
  lprt_row   =  $LPRT_ROW
  wrf_lc_ref_lat = $WRF_LC_REF_LAT
 $Marker

 &WINDOWDEFS
  x0         =  $X0
  y0         =  $Y0
  ncolsin    =  $NCOLS
  nrowsin    =  $NROWS
 $Marker

!

rm fort.*
if ( -f $FILE_GD ) rm -f $FILE_GD

ln -s $FILE_GD                   fort.4
ln -s $WorkDir/namelist.$PROG  fort.8

set NUMFIL = 0
foreach fil ( $InMetFiles )
  @ NN = $NUMFIL + 10
  ln -s $fil fort.$NN
  @ NUMFIL ++
end

setenv IOAPI_CHECK_HEADERS  T
setenv EXECUTION_ID         $PROG

setenv GRID_BDY_2D $OutDir/GRIDBDY2D_$APPL.nc
setenv GRID_CRO_2D $OutDir/GRIDCRO2D_$APPL.nc
setenv GRID_DOT_2D $OutDir/GRIDDOT2D_$APPL.nc
setenv MET_BDY_3D  $OutDir/METBDY3D_$APPL.nc
setenv MET_CRO_2D  $OutDir/METCRO2D_$APPL.nc
setenv MET_CRO_3D  $OutDir/METCRO3D_$APPL.nc
setenv MET_DOT_3D  $OutDir/METDOT3D_$APPL.nc
setenv LUFRAC_CRO  $OutDir/LUFRAC_CRO_$APPL.nc
setenv SOI_CRO     $OutDir/SOI_CRO_$APPL.nc
setenv MOSAIC_CRO  $OutDir/MOSAIC_CRO_$APPL.nc

if ( -f $GRID_BDY_2D ) rm -f $GRID_BDY_2D
if ( -f $GRID_CRO_2D ) rm -f $GRID_CRO_2D
if ( -f $GRID_DOT_2D ) rm -f $GRID_DOT_2D
if ( -f $MET_BDY_3D  ) rm -f $MET_BDY_3D
if ( -f $MET_CRO_2D  ) rm -f $MET_CRO_2D
if ( -f $MET_CRO_3D  ) rm -f $MET_CRO_3D
if ( -f $MET_DOT_3D  ) rm -f $MET_DOT_3D
if ( -f $LUFRAC_CRO  ) rm -f $LUFRAC_CRO
if ( -f $SOI_CRO     ) rm -f $SOI_CRO
if ( -f $MOSAIC_CRO  ) rm -f $MOSAIC_CRO

if ( -f $OutDir/mcip.nc      ) rm -f $OutDir/mcip.nc
if ( -f $OutDir/mcip_bdy.nc  ) rm -f $OutDir/mcip_bdy.nc

$ProgDir/$PROG.exe

if ( $status == 0 ) then
  rm fort.*
  exit 0
else
  echo "Error running $PROG"
  exit 1
endif""".format(compiler, year, month, day, dom.size, dom.name, dom.id2,
                proj_name, dir_in_met, dir_in_geo, dir_out, dir_prog,
                in_met_files, mcip_start, mcip_end, dom.ncol, dom.nrow)
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


def get_InMetFiles(days):
    """ get input meteorology file paths as list """
    fmt = wrfout_fmt + ' \\'
    return [fmt.format(d) for d in days]


def create_InMetFiles(days):
    """ Create input meteorology file paths string """
    return 'set InMetFiles = ( {})'.format(
        '\n\t'.join(get_InMetFiles(days))[:-1])


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'mcip script\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -d 4 5 6\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-n', '--domain', nargs='+', type=int,
                   default=proj.get_dom_ids(),
                   help="domain Id(s). Default is all domains in the project.")
    p.add_argument('--monthly', action='store_true',
                   help='Run mcip as monthly')
    p.add_argument('-y', '--years', nargs='+', type=int, default=proj.years)
    p.add_argument('-m', '--months', nargs='+', type=int, default=proj.months)
    p.add_argument('-d', '--days', nargs='+', type=int, default=proj.days)
    return p.parse_args()


if __name__ == "__main__":
    from _helper_functions_ import ExitHelper
    from _helper_functions_ import ScriptError
    from _helper_functions_ import run_script_mcip
    from _helper_functions_ import del_files

    a = _parse_args_()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    year, month, day = a.years, a.months, a.days
    doms = [proj.get_dom_by_id(i) for i in a.domain]
    run_name = 'monthly' if a.monthly else 'daily'
    if a.monthly:
        day = [1]
    ym = expandgrid(year, month)  # Year and months

    log_dir = _join(proj.path.logs, 'mcip')
    os.makedirs(log_dir, exist_ok=True)

    log.info('Creating {} mcip files'.format(run_name))
    log.info('log dir: {}'.format(log_dir))

    time_fmt = '[Time: {:.2f} secs]'
    mcip_timer_start = timer()
    flag = ExitHelper()
    for dom in doms:
        dom_timer_start = timer()
        dom_str = 'Dom: {}-{}'.format(dom.size, dom.name)
        log.info(dom_str + ' starting...')
        for y, m in ym:
            month_str = 'Month: {}-{:02d}'.format(y, m)
            if not a.monthly:
                log.info('Processing {}'.format(month_str))

            days = get_days(y, m, day)

            if a.monthly:
                days_tmp = get_days(y, m, list(range(1, 32)))
                in_met_files = create_InMetFiles(days_tmp)
            else:
                in_met_files = create_InMetFiles(days)

            mn = calendar.month_name[m].lower()
            dir_out = dir_out_fmt.format(dom.size, dom.name, mn)
            if a.monthly:
                dir_out += '_monthly'
            os.makedirs(dir_out, exist_ok=True)
            log.info('Output Dir: {}'.format(dir_out))
            dir_in_met = dir_in_met_fmt.format(mn)

            month_timer_start = timer()
            for d in days:
                if flag.exit:
                    log.info('User stopped execution')
                    sys.exit()
                str_date = '{:04d}-{:02d}-{:02d}'.format(y, m, d)
                day_str = 'Day: {}'.format(str_date)
                script = get_script(y, m, d, dom, proj.name, dir_in_met,
                                    dir_in_geo, dir_out, dir_prog,
                                    in_met_files, a.monthly,
                                    proj.compiler)

                day_timer_start = timer()
                err = None
                try:
                    run_script_mcip(script, dom, str_date)
                except ScriptError as error:
                    err = error
                day_timer_end = timer()

                if not a.monthly:
                    elapsed = day_timer_end - day_timer_start
                    msg = dom_str + ', ' + day_str
                    if err is not None:
                        msg = msg + ', ' + err.detail
                    msg = msg + ' ' + time_fmt.format(elapsed)
                    if err is None:
                        log.info(msg)
                    else:
                        log.error(msg)
                        err = None

            del_files(dir_out, 'fort.*')

            elapsed = day_timer_end - month_timer_start
            msg = dom_str + ', ' + month_str
            if err is not None:
                msg = msg + ', ' + err.detail
            msg = msg + ' ' + time_fmt.format(elapsed)
            if err is None:
                log.info(msg)
            else:
                log.error(msg)
                err = None

        msg = '{} completed '.format(dom_str)
        elapsed = day_timer_end - dom_timer_start
        log.info(msg + time_fmt.format(elapsed))

    elapsed = day_timer_end - mcip_timer_start
    log.info('mcip completed ' + time_fmt.format(elapsed))
