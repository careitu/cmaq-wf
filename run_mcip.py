#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703


"""
Run mcip
~~~~~~~~
Python script to run mcip
"""

import os
import itertools
import logging
from os.path import join

import calendar
from collections.abc import Iterable
from datetime import datetime as _dt
from datetime import timedelta as _td
from datetime import timezone as _tz

from settings import setting as s
from _helper_functions_ import _create_argparser_

proj = s.get_active_project()
logging_dir = os.path.join(os.path.expanduser("~"), '.config/cwf')
os.makedirs(logging_dir, exist_ok=True)

# ----------------------------------
log = logging.getLogger('mcip')
log.setLevel(logging.DEBUG)
# create formatter
fmt_str = '%(asctime)s.%(msecs)03d - %(name)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(fmt_str, datefmt='%Y-%m-%dT%H:%M:%S')
# create file handler
fh = logging.FileHandler(os.path.join(logging_dir, 'cwf.log'))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
log.addHandler(fh)

# ----------------------------------

dir_proj = join(proj.path, proj.name)
dir_prog = join(proj.path_cmaq, 'PREP/mcip/src')
wrfout_fmt = '${{InMetDir}}' \
             '/wrfout_${{dom_num}}_${{year}}-${{month}}-{:02d}_00:00:00'
dir_in_geo = join(dir_proj, 'WPS')
dir_out_fmt = join(dir_proj, 'mcip/{}km/{}/{}')
dir_in_met_fmt = join(dir_proj, 'wrf/{}')


def get_script(year, month, day, dom, proj_name, region, dir_in_met,
               dir_in_geo, dir_out, dir_prog, in_met_files, monthly=False,
               compiler='gcc'):
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
set region = {}

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
endif""".format(compiler, year, month, day, dom.size, dom.name, dom.id,
                proj_name, region, dir_in_met, dir_in_geo, dir_out, dir_prog,
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


def create_InMetFiles(days):
    """ Create input meteorology file paths """
    fmt = wrfout_fmt + ' \\'
    list_of_files = [fmt.format(d) for d in days]
    str_files = '\n\t'.join(list_of_files)[:-1]
    return 'set InMetFiles = ( {})'.format(str_files)


if __name__ == "__main__":
    DESCRIPTION = 'Run mcip script\n' + \
                  'Project: {}\nProject Dir: {}\nCMAQ Dir: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path, proj.path_cmaq)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -d \n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-n', '--domain', help="domain Id(s)", nargs='+',
                   type=int, default=[d.id for d in proj.doms])
    p.add_argument('--monthly', action='store_true',
                   help='foo the bars before frobbling')
    p.add_argument('-y', '--years', nargs='+', type=int, default=proj.years)
    p.add_argument('-m', '--months', nargs='+', type=int, default=proj.months)
    p.add_argument('-d', '--days', nargs='+', type=int, default=proj.days)
    a = p.parse_args()

    verbose = a.print
    if verbose:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    year, month, day = a.years, a.months, a.days
    doms = [proj.doms[i - 1] for i in a.domain]
    run_name = 'monthly' if a.monthly else 'daily'
    if a.monthly:
        day = [1]
    ym = expandgrid(year, month)  # Year and months

    log_dir = os.path.join(dir_proj, 'logs', 'mcip')
    os.makedirs(log_dir, exist_ok=True)

    import tempfile as tf
    import subprocess as subp
    for dom in doms:
        log.info('Creating {} mcip files for domain: {}'.format(run_name,
                                                                dom.name))
        for y, m in ym:
            log.info('Processing Year: {}, Month: {}'.format(y, m))
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

            for d in days:
                str_date = '{:04d}-{:02d}-{:02d}'.format(y, m, d)
                if not a. monthly:
                    log.info('Processing day: {}'.format(str_date))
                # subp.call('rm ' + os.path.join(dir_out, 'fort*'), shell=True)
                # pylint: disable=W0212
                tmp_name = next(tf._get_candidate_names())
                fmt = 'mcip_{}_{}'.format(dom.name, str_date)
                file_script = '{}_{}.csh'.format(fmt, tmp_name)
                file_script = os.path.join(tf.gettempdir(), file_script)
                file_log = os.path.join(log_dir, '{}.log'.format(fmt))
                file_err = os.path.join(log_dir, '{}.err'.format(fmt))
                script = get_script(y, m, d, dom, proj.name, dom.name,
                                    dir_in_met, dir_in_geo, dir_out,
                                    dir_prog, in_met_files, a.monthly,
                                    proj.compiler)
                with open(file_script, 'w') as f:
                    f.write("#!/bin/csh -f\n")
                    f.write(script)

                subp.call(['chmod', '+x', file_script])
                with open(file_err, 'w') as fe:
                    with open(file_log, 'w') as fl:
                        subp.call([file_script], stdout=fl, stderr=fe)

                os.remove(file_script)
            # subp.call('rm ' + os.path.join(dir_out, 'fort*'), shell=True)
