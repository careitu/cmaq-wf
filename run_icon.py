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

if ( ! -e $CMAQ_DATA ) then
    echo "   $CMAQ_DATA path does not exist"
    exit 1
endif
echo " "; echo " Input data path, CMAQ_DATA set to $CMAQ_DATA"; echo " "

set year = {}
set month = {:02d} # march
set domainsize = {}km
set domainsize_larger = {}km # one upper level domain size
set domain_num = d{:02d}
set project_name = {}
set region = {}
set proj_mcip_path = ${{proj_path}}/mcip
set proj_path = /mnt/ssd2/APPS/CMAQ/data/projects/CityAir
set upper_dom_path = ${{proj_path}}/${domainsize_larger}

set ymd      = ${{year}}${{month}}${{day}}
set APPL     = ${{project_name}}_${{dom_size}}_${{dom_name}}_${{ymd}}
set VRSN     = v532
set ICTYPE   = regrid

set BLD      = ${CMAQ_HOME}/PREP/icon/scripts/BLD_ICON_${VRSN}_${compiler}
set EXEC     = ICON_${VRSN}.exe
cat $BLD/ICON_${VRSN}.cfg; echo " "; set echo

setenv GRID_NAME 4km
setenv GRIDDESC  ${{proj_mcip_path}}/${domainsize}/${region}/${month}/GRIDDESC
setenv IOAPI_ISPH 20

setenv IOAPI_LOG_WRITE F
setenv IOAPI_OFFSET_64 YES
setenv EXECUTION_ID $EXEC

setenv ICON_TYPE ` echo $ICTYPE | tr "[A-Z]" "[a-z]" `

set OUTDIR   = ${{proj_path}}/icon

set DATE = "2015-03-02"
set YYYYJJJ  = `date -ud "${DATE}" +%Y%j`   #> Convert YYYY-MM-DD to YYYYJJJ
set YYMMDD   = `date -ud "${DATE}" +%y%m%d` #> Convert YYYY-MM-DD to YYMMDD
set YYYYMMDD = `date -ud "${DATE}" +%Y%m%d` #> Convert YYYY-MM-DD to YYYYMMDD

if ( $ICON_TYPE == regrid ) then
    setenv CTM_CONC_1 ${{upper_dom_path}}/CCTM_CONC_v532_gcc_CityAir_2015_${domainsize_larger}_${YYYYMMDD}.nc
    setenv MET_CRO_3D_CRS ${{proj_mcip_path}}/${domainsize_larger}/${month}/METCRO3D_CityAir_${domainsize_larger}_${YYYYMMDD}.nc
    setenv MET_CRO_3D_FIN ${{proj_mcip_path}}/${domainsize}/${region}/${month}/METCRO3D_CityAir_${domainsize}_${YYYYMMDD}.nc
    setenv INIT_CONC_1    "$OUTDIR/ICON_${VRSN}_${APPL}_${ICON_TYPE}_${YYYYMMDD} -v"
endif

if ( $ICON_TYPE == profile ) then
    setenv IC_PROFILE $BLD/profiles/avprofile_cb6r3m_ae7_kmtbr_hemi2016_v53beta2_m3dry_col051_row068.csv
    setenv MET_CRO_3D_FIN /ssd2/programs/cmaq/data/mcip/36km/March_2015/METCRO3D_cityair_march_2015_36km.nc
    setenv INIT_CONC_1    "$OUTDIR/ICON_${VRSN}_${APPL}_${ICON_TYPE}_${YYYYMMDD} -v"
endif

if ( ! -d "$OUTDIR" ) mkdir -p $OUTDIR

ls -l $BLD/$EXEC; size $BLD/$EXEC
unlimit
limit

time $BLD/$EXEC

exit()""".format(compiler, year, month, day, dom.size, dom.name, dom.id,
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
