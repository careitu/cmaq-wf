#! /usr/bin/env python
# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0621,W0702,W0703


"""
Run cctm
~~~~~~~~
Python script to run cctm
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
log = logging.getLogger('{}.cctm'.format(proj.name))
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


def get_script(year, month, start_day, end_day, dom, proc_type='mpi',
               new_start=True, ncores='auto', icon_type='auto',
               bcon_type='auto'):
    import psutil
    from _helper_functions_ import max_mult

    OUTDIR = _join(proj.path.cctm, '{}km'.format(dom.size), dom.name)
    NEW_START = 'TRUE' if new_start else 'FALSE'
    PROC = 'serial' if proc_type == 'serial' else 'mpi'
    START_DATE = '{:04d}-{:02d}-{:02d}'.format(year, month, start_day)
    END_DATE = '{:04d}-{:02d}-{:02d}'.format(year, month, end_day)

    if isinstance(ncores, str):
        if ncores == 'auto':
            ncores = psutil.cpu_count(logical=False)
    if isinstance(ncores, int):
        NPROW, NPCOL = max_mult(ncores)
    else:
        raise ValueError("ncores must be a positive integer or 'auto'")

    icon_type = icon_type.lower()
    if icon_type == 'auto':
        icon_type = 'regrid' if dom.__parent__ is not None else 'profile'

    bcon_type = bcon_type.lower()
    if bcon_type == 'auto':
        bcon_type = 'regrid' if dom.__parent__ is not None else 'profile'

    mn = calendar.month_name[month].lower()
    script = """
setenv compiler {}
echo 'Start Model Run At ' `date`

setenv CTM_DIAG_LVL 0

if ( ! $?compiler ) then
  setenv compiler gcc
endif
if ( ! $?compilerVrsn ) then
  setenv compilerVrsn Empty
endif

set CURDIR = $PWD
source /mnt/ssd2/APPS/CMAQ/config_cmaq.csh $compiler $compilerVrsn
cd $CURDIR

set month_name = {}
set proj_name = {}
set dom_name = {}
set dom_size2 = {}
set dom_size = ${{dom_size2}}km

set ICTYPE = {}
set BCTYPE = {}
set VRSN = v{}
set PROC = {}
set MECH = cb6r3_ae7_aq
set APPL = ${{proj_name}}_2015_${{dom_size}}

setenv RUNID  ${{VRSN}}_${{compilerString}}_${{APPL}}

set BLD  = ${{CMAQ_HOME}}/CCTM/scripts/BLD_CCTM_${{VRSN}}_${{compilerString}}
set EXEC = CCTM_${{VRSN}}.exe

if ( $CTM_DIAG_LVL != 0 ) set echo

setenv WORKDIR ${{CMAQ_HOME}}/CCTM/scripts
setenv OUTDIR  {}
setenv INPDIR  {}
setenv MCIPDIR {}
setenv LOGDIR  ${{OUTDIR}}/LOGS
setenv NMLpath ${{BLD}}

echo ""
echo "Current Dir: $CURDIR"
echo "Working Directory is $WORKDIR"
echo "Build Directory is $BLD"
echo "Output Directory is $OUTDIR"
echo "Log Directory is $LOGDIR"
echo "Executable Name is $EXEC"

setenv NEW_START {}
set START_DATE = "{}"
set END_DATE   = "{}"

set STTIME = 000000
set NSTEPS = 250000
set TSTEP  = 010000

if ( $PROC == serial ) then
  setenv NPCOL_NPROW "1 1"; set NPROCS   = 1
else
  # number of cores nxm
  @ NPCOL  =  {}; @ NPROW =  {}
  @ NPROCS = $NPCOL * $NPROW
  setenv NPCOL_NPROW "$NPCOL $NPROW";
endif

setenv EXECUTION_ID "CMAQ_CCTM_${{VRSN}}_`id -u -n`_`date -u +%Y%m%d_%H%M%S_%N`"

echo ""
echo "---CMAQ EXECUTION ID: $EXECUTION_ID ---"

echo "HELLOOO"

set CLOBBER_DATA = TRUE

if (! -e $LOGDIR ) then
  mkdir -p $LOGDIR
endif
setenv PRINT_PROC_TIME Y
setenv STDOUT T

setenv GRID_NAME ${{dom_size}}
setenv GRIDDESC ${{MCIPDIR}}/${{dom_size}}/${{dom_name}}/${{month_name}}/GRIDDESC

set NZ = 35
set NX = `grep -A 1 ${{GRID_NAME}} ${{GRIDDESC}} | tail -1 | sed 's/  */ /g' | cut -d' ' -f6`
set NY = `grep -A 1 ${{GRID_NAME}} ${{GRIDDESC}} | tail -1 | sed 's/  */ /g' | cut -d' ' -f7`
set NCELLS = `echo "${{NX}} * ${{NY}} * ${{NZ}}" | bc -l`

setenv CONC_SPCS "ALL"
setenv AVG_CONC_SPCS "ALL"
setenv ACONC_BLEV_ELEV " 1 1"
setenv AVG_FILE_ENDTIME N

setenv CTM_MAXSYNC 300
setenv CTM_MINSYNC  60
setenv SIGMA_SYNC_TOP 0.7
setenv CTM_ADV_CFL 0.95

setenv CTM_OCEAN_CHEM Y
setenv CTM_WB_DUST N
setenv CTM_WBDUST_BELD BELD3
setenv CTM_LTNG_NO N
setenv KZMIN Y
setenv CTM_MOSAIC N
setenv CTM_FST N
setenv PX_VERSION Y
setenv CLM_VERSION N
setenv NOAH_VERSION N
setenv CTM_ABFLUX N
setenv CTM_BIDI_FERT_NH3 T
setenv CTM_HGBIDI N
setenv CTM_SFC_HONO Y
setenv CTM_GRAV_SETL Y
setenv CTM_BIOGEMIS N

setenv VERTEXT N
setenv VERTEXT_COORD_PATH ${{WORKDIR}}/lonlat.csv

setenv IOAPI_LOG_WRITE F
setenv FL_ERR_STOP N
setenv PROMPTFLAG F
setenv IOAPI_OFFSET_64 YES
setenv IOAPI_CHECK_HEADERS N
setenv CTM_EMISCHK N
setenv EMISDIAG F
setenv EMISDIAG_SUM F

setenv EMIS_SYM_DATE N
setenv CTM_CKSUM Y
setenv CLD_DIAG N

setenv CTM_PHOTDIAG N
setenv NLAYS_PHOTDIAG "1"

setenv CTM_PMDIAG N
setenv CTM_APMDIAG Y
setenv APMDIAG_BLEV_ELEV "1 1"

setenv CTM_SSEMDIAG N
setenv CTM_DUSTEM_DIAG N
setenv CTM_DEPV_FILE N
setenv VDIFF_DIAG_FILE N
setenv LTNGDIAG N
setenv B3GTS_DIAG N
setenv CTM_WVEL Y

set ICpath    = {}
set BCpath    = {}
set EMISpath  = {}/${{dom_size}}/${{dom_name}}/${{month_name}}

set IN_LTpath = $INPDIR/lightning
set METpath   = ${{MCIPDIR}}/${{dom_size}}/${{dom_name}}/${{month_name}}

set OMIpath   = $BLD
set LUpath    = $INPDIR/land
set SZpath    = $INPDIR/land

set rtarray = ""

set TODAYG = ${{START_DATE}}
set TODAYJ = `date -ud "${{START_DATE}}" +%Y%j`
set START_DAY = ${{TODAYJ}}
set STOP_DAY = `date -ud "${{END_DATE}}" +%Y%j`
set NDAYS = 31

while ($TODAYJ <= $STOP_DAY )

  set NDAYS = `echo "${{NDAYS}} + 1" | bc -l`

  set YYYYMMDD = `date -ud "${{TODAYG}}" +%Y%m%d`
  set YYYYMM = `date -ud "${{TODAYG}}" +%Y%m`
  set YYMMDD = `date -ud "${{TODAYG}}" +%y%m%d`
  set YYYYJJJ = $TODAYJ

  set YESTERDAY = `date -ud "${{TODAYG}}-1days" +%Y%m%d`

  echo ""
  echo "Set up input and output files for Day ${{TODAYG}}."

  setenv CTM_APPL ${{RUNID}}_${{YYYYMMDD}}
  set APPL2 = ${{proj_name}}_${{dom_size}}_${{dom_name}}_${{YYYYMMDD}}

  if ( ! -d "$OUTDIR" ) mkdir -p $OUTDIR
  cp $BLD/CCTM_${{VRSN}}.cfg $OUTDIR/CCTM_${{CTM_APPL}}.cfg


  if ($NEW_START == true || $NEW_START == TRUE ) then
     setenv ICFILE ICON_${{VRSN}}_${{ICTYPE}}_${{APPL2}}
     setenv INIT_MEDC_1 notused
     setenv INITIAL_RUN Y
  else
     set ICpath = $OUTDIR
     setenv ICFILE CCTM_CONC_${{RUNID}}_${{YESTERDAY}}.nc
     setenv INIT_MEDC_1 $ICpath/CCTM_MEDIA_CONC_${{RUNID}}_${{YESTERDAY}}.nc
     setenv INITIAL_RUN N
  endif

  set BCFILE = BCON_${{VRSN}}_${{BCTYPE}}_${{APPL2}}

  set OMIfile = OMI_1979_to_2019.dat

  set OPTfile = PHOT_OPTICS.dat

  setenv GRID_BDY_2D $METpath/GRIDBDY2D_${{APPL2}}.nc
  setenv GRID_CRO_2D $METpath/GRIDCRO2D_${{APPL2}}.nc
  setenv GRID_CRO_3D $METpath/GRIDCRO3D_${{APPL2}}.nc
  setenv GRID_DOT_2D $METpath/GRIDDOT2D_${{APPL2}}.nc
  setenv MET_CRO_2D $METpath/METCRO2D_${{APPL2}}.nc
  setenv MET_CRO_3D $METpath/METCRO3D_${{APPL2}}.nc
  setenv MET_DOT_3D $METpath/METDOT3D_${{APPL2}}.nc
  setenv MET_BDY_3D $METpath/METBDY3D_${{APPL2}}.nc
  setenv LUFRAC_CRO $METpath/LUFRAC_CRO_${{APPL2}}.nc

  setenv EMISSCTRL_NML ${{BLD}}/EmissCtrl_${{MECH}}.nml

  setenv CMAQ_MASKS $SZpath/ocean_${{dom_size}}_${{proj_name}}_${{dom_name}}.nc

  setenv N_EMIS_GR 1
  set EMISfile  = ${{dom_name}}_${{dom_size2}}_CB6_aero7_${{YYYYMMDD}}.nc
  setenv GR_EMIS_001 ${{EMISpath}}/${{EMISfile}}
  setenv GR_EMIS_LAB_001 GRIDDED_EMIS
  setenv GR_EM_SYM_DATE_001 F

  setenv N_EMIS_PT 0

  if ( $CTM_LTNG_NO == 'Y' ) then
     setenv LTNGNO "InLine"

     setenv USE_NLDN  Y
     if ( $USE_NLDN == Y ) then
        setenv NLDN_STRIKES ${{IN_LTpath}}/NLDN.12US1.${{YYYYMMDD}}_bench.nc
     endif
     setenv LTNGPARMS_FILE ${{IN_LTpath}}/LTNG_AllParms_12US1_bench.nc
  endif

  if ( $CTM_BIOGEMIS == 'Y' ) then
     set IN_BEISpath = ${{INPDIR}}/land
     setenv GSPRO      $BLD/gspro_biogenics.txt
     setenv B3GRD      $IN_BEISpath/b3grd_bench.nc
     setenv BIOSW_YN   Y
     setenv BIOSEASON  $IN_BEISpath/bioseason.cmaq.2016_12US1_full_bench.ncf

     setenv SUMMER_YN  Y
     setenv PX_VERSION Y
     setenv SOILINP    $OUTDIR/CCTM_SOILOUT_${{RUNID}}_${{YESTERDAY}}.nc
  endif

  if ( $CTM_WB_DUST == 'Y' ) then
     setenv DUST_LU_1 $LUpath/beld3_12US1_459X299_output_a_bench.nc
     setenv DUST_LU_2 $LUpath/beld4_12US1_459X299_output_tot_bench.nc
  endif

  setenv OCEAN_1 $SZpath/ocean_${{dom_size}}_${{proj_name}}_${{dom_name}}.nc

  if ( $CTM_ABFLUX == 'Y' ) then
     setenv E2C_SOIL ${{LUpath}}/epic_festc1.4_20180516/2016_US1_soil_bench.nc
     setenv E2C_CHEM ${{LUpath}}/epic_festc1.4_20180516/2016_US1_time${{YYYYMMDD}}_bench.nc
     setenv E2C_CHEM_YEST ${{LUpath}}/epic_festc1.4_20180516/2016_US1_time${{YESTERDAY}}_bench.nc
     setenv E2C_LU ${{LUpath}}/beld4_12kmCONUS_2011nlcd_bench.nc
  endif

  setenv CTM_PROCAN N
  if ( $?CTM_PROCAN ) then
     if ( $CTM_PROCAN == 'Y' || $CTM_PROCAN == 'T' ) then
        setenv PACM_INFILE ${{NMLpath}}/pa_${{MECH}}.ctl
        setenv PACM_REPORT $OUTDIR/"PA_REPORT".${{YYYYMMDD}}
     endif
  endif

 setenv CTM_ISAM N
 if ( $?CTM_ISAM ) then
    if ( $CTM_ISAM == 'Y' || $CTM_ISAM == 'T' ) then
       setenv SA_IOLIST ${{WORKDIR}}/isam_control.txt
       setenv ISAM_BLEV_ELEV " 1 1"
       setenv AISAM_BLEV_ELEV " 1 1"

       if ($NEW_START == true || $NEW_START == TRUE ) then
          setenv ISAM_NEW_START Y
          setenv ISAM_PREVDAY
       else
          setenv ISAM_NEW_START N
          setenv ISAM_PREVDAY "$OUTDIR/CCTM_SA_CGRID_${{RUNID}}_${{YESTERDAY}}.nc"
       endif

       setenv SA_ACONC_1 "$OUTDIR/CCTM_SA_ACONC_${{CTM_APPL}}.nc -v"
       setenv SA_CONC_1  "$OUTDIR/CCTM_SA_CONC_${{CTM_APPL}}.nc -v"
       setenv SA_DD_1    "$OUTDIR/CCTM_SA_DRYDEP_${{CTM_APPL}}.nc -v"
       setenv SA_WD_1    "$OUTDIR/CCTM_SA_WETDEP_${{CTM_APPL}}.nc -v"
       setenv SA_CGRID_1 "$OUTDIR/CCTM_SA_CGRID_${{CTM_APPL}}.nc -v"

    endif
 endif


 setenv STM_SO4TRACK N
 if ( $?STM_SO4TRACK ) then
    if ( $STM_SO4TRACK == 'Y' || $STM_SO4TRACK == 'T' ) then
      setenv STM_ADJSO4 Y
    endif
 endif

  setenv S_CGRID         "$OUTDIR/CCTM_CGRID_${{CTM_APPL}}.nc"
  setenv CTM_CONC_1      "$OUTDIR/CCTM_CONC_${{CTM_APPL}}.nc -v"
  setenv A_CONC_1        "$OUTDIR/CCTM_ACONC_${{CTM_APPL}}.nc -v"
  setenv MEDIA_CONC      "$OUTDIR/CCTM_MEDIA_CONC_${{CTM_APPL}}.nc -v"
  setenv CTM_DRY_DEP_1   "$OUTDIR/CCTM_DRYDEP_${{CTM_APPL}}.nc -v"
  setenv CTM_DEPV_DIAG   "$OUTDIR/CCTM_DEPV_${{CTM_APPL}}.nc -v"
  setenv B3GTS_S         "$OUTDIR/CCTM_B3GTS_S_${{CTM_APPL}}.nc -v"
  setenv SOILOUT         "$OUTDIR/CCTM_SOILOUT_${{CTM_APPL}}.nc"
  setenv CTM_WET_DEP_1   "$OUTDIR/CCTM_WETDEP1_${{CTM_APPL}}.nc -v"
  setenv CTM_WET_DEP_2   "$OUTDIR/CCTM_WETDEP2_${{CTM_APPL}}.nc -v"
  setenv CTM_PMDIAG_1    "$OUTDIR/CCTM_PMDIAG_${{CTM_APPL}}.nc -v"
  setenv CTM_APMDIAG_1   "$OUTDIR/CCTM_APMDIAG_${{CTM_APPL}}.nc -v"
  setenv CTM_RJ_1        "$OUTDIR/CCTM_PHOTDIAG1_${{CTM_APPL}}.nc -v"
  setenv CTM_RJ_2        "$OUTDIR/CCTM_PHOTDIAG2_${{CTM_APPL}}.nc -v"
  setenv CTM_RJ_3        "$OUTDIR/CCTM_PHOTDIAG3_${{CTM_APPL}}.nc -v"
  setenv CTM_SSEMIS_1    "$OUTDIR/CCTM_SSEMIS_${{CTM_APPL}}.nc -v"
  setenv CTM_DUST_EMIS_1 "$OUTDIR/CCTM_DUSTEMIS_${{CTM_APPL}}.nc -v"
  setenv CTM_IPR_1       "$OUTDIR/CCTM_PA_1_${{CTM_APPL}}.nc -v"
  setenv CTM_IPR_2       "$OUTDIR/CCTM_PA_2_${{CTM_APPL}}.nc -v"
  setenv CTM_IPR_3       "$OUTDIR/CCTM_PA_3_${{CTM_APPL}}.nc -v"
  setenv CTM_IRR_1       "$OUTDIR/CCTM_IRR_1_${{CTM_APPL}}.nc -v"
  setenv CTM_IRR_2       "$OUTDIR/CCTM_IRR_2_${{CTM_APPL}}.nc -v"
  setenv CTM_IRR_3       "$OUTDIR/CCTM_IRR_3_${{CTM_APPL}}.nc -v"
  setenv CTM_DRY_DEP_MOS "$OUTDIR/CCTM_DDMOS_${{CTM_APPL}}.nc -v"
  setenv CTM_DRY_DEP_FST "$OUTDIR/CCTM_DDFST_${{CTM_APPL}}.nc -v"
  setenv CTM_DEPV_MOS    "$OUTDIR/CCTM_DEPVMOS_${{CTM_APPL}}.nc -v"
  setenv CTM_DEPV_FST    "$OUTDIR/CCTM_DEPVFST_${{CTM_APPL}}.nc -v"
  setenv CTM_VDIFF_DIAG  "$OUTDIR/CCTM_VDIFF_DIAG_${{CTM_APPL}}.nc -v"
  setenv CTM_VSED_DIAG   "$OUTDIR/CCTM_VSED_DIAG_${{CTM_APPL}}.nc -v"
  setenv CTM_LTNGDIAG_1  "$OUTDIR/CCTM_LTNGHRLY_${{CTM_APPL}}.nc -v"
  setenv CTM_LTNGDIAG_2  "$OUTDIR/CCTM_LTNGCOL_${{CTM_APPL}}.nc -v"
  setenv CTM_VEXT_1      "$OUTDIR/CCTM_VEXT_${{CTM_APPL}}.nc -v"


  setenv FLOOR_FILE ${{OUTDIR}}/FLOOR_${{CTM_APPL}}.txt

  ( ls CTM_LOG_???.${{CTM_APPL}} > $CURDIR/buff.txt ) >& /dev/null
  ( ls ${{LOGDIR}}/CTM_LOG_???.${{CTM_APPL}} >> $CURDIR/buff.txt ) >& /dev/null

  set log_test = `cat $CURDIR/buff.txt`; rm -f $CURDIR/buff.txt

  set OUT_FILES = ($FLOOR_FILE $S_CGRID $CTM_CONC_1 $A_CONC_1 \
             $MEDIA_CONC $CTM_DRY_DEP_1 $CTM_DEPV_DIAG $B3GTS_S \
             $SOILOUT $CTM_WET_DEP_1 $CTM_WET_DEP_2 $CTM_PMDIAG_1 \
             $CTM_APMDIAG_1 $CTM_RJ_1 $CTM_RJ_2 $CTM_RJ_3 $CTM_SSEMIS_1 \
             $CTM_DUST_EMIS_1 $CTM_IPR_1 $CTM_IPR_2 $CTM_IPR_3 $CTM_IRR_1 \
             $CTM_IRR_2 $CTM_IRR_3 $CTM_DRY_DEP_MOS $CTM_DRY_DEP_FST \
             $CTM_DEPV_MOS $CTM_DEPV_FST $CTM_VDIFF_DIAG $CTM_VSED_DIAG \
             $CTM_LTNGDIAG_1 $CTM_LTNGDIAG_2 $CTM_VEXT_1 )
  if ( $?CTM_ISAM ) then
     if ( $CTM_ISAM == 'Y' || $CTM_ISAM == 'T' ) then
        set OUT_FILES = (${{OUT_FILES}} ${{SA_ACONC_1}} ${{SA_CONC_1}} \
                         ${{SA_DD_1}} ${{SA_WD_1}} ${{SA_CGRID_1}} )
     endif
  endif

  set OUT_FILES = `echo $OUT_FILES | sed "s; -v;;g" | sed "s;MPI:;;g" `
  ( ls $OUT_FILES > $CURDIR/buff.txt ) >& /dev/null
  set out_test = `cat $CURDIR/buff.txt`; rm -f $CURDIR/buff.txt

  if ( $CLOBBER_DATA == true || $CLOBBER_DATA == TRUE  ) then
     echo
     echo "Existing Logs and Output Files for Day ${{TODAYG}} Will Be Deleted"

     foreach file ( ${{log_test}} )
        /bin/rm -f $file
     end

     foreach file ( ${{out_test}} )
        /bin/rm -f $file
     end
     /bin/rm -f ${{OUTDIR}}/CCTM_EMDIAG*${{RUNID}}_${{YYYYMMDD}}.nc

  else
     if ( "$log_test" != "" ) then
       echo "*** Logs exist - run ABORTED ***"
       echo "*** To overide, set CLOBBER_DATA = TRUE in run_cctm.csh ***"
       echo "*** and these files will be automatically deleted. ***"
       exit 1
     endif

     if ( "$out_test" != "" ) then
       echo "*** Output Files Exist - run will be ABORTED ***"
       foreach file ( $out_test )
          echo " cannot delete $file"
       end
       echo "*** To overide, set CLOBBER_DATA = TRUE in run_cctm.csh ***"
       echo "*** and these files will be automatically deleted. ***"
       exit 1
     endif
  endif

  setenv CTM_STDATE      $YYYYJJJ
  setenv CTM_STTIME      $STTIME
  setenv CTM_RUNLEN      $NSTEPS
  setenv CTM_TSTEP       $TSTEP
  setenv INIT_CONC_1 $ICpath/$ICFILE
  setenv BNDY_CONC_1 $BCpath/$BCFILE
  setenv OMI $OMIpath/$OMIfile
  setenv OPTICS_DATA $OMIpath/$OPTfile
  set TR_DVpath = $METpath
  set TR_DVfile = $MET_CRO_2D

  setenv gc_matrix_nml ${{NMLpath}}/GC_$MECH.nml
  setenv ae_matrix_nml ${{NMLpath}}/AE_$MECH.nml
  setenv nr_matrix_nml ${{NMLpath}}/NR_$MECH.nml
  setenv tr_matrix_nml ${{NMLpath}}/Species_Table_TR_0.nml

  setenv CSQY_DATA ${{NMLpath}}/CSQY_DATA_$MECH

  if (! (-e $CSQY_DATA ) ) then
     echo " $CSQY_DATA  not found "
     exit 1
  endif
  if (! (-e $OPTICS_DATA ) ) then
     echo " $OPTICS_DATA  not found "
     exit 1
  endif

  if ( $CTM_DIAG_LVL != 0 ) then
     ls -l $BLD/$EXEC
     size $BLD/$EXEC
     unlimit
     limit
  endif

  echo
  echo "CMAQ Processing of Day $YYYYMMDD Began at `date`"
  echo

  set MPI = =/usr/bin
  set MPIRUN = $MPI/mpirun
  ( /usr/bin/time -p mpirun -np $NPROCS --use-hwthread-cpus $BLD/$EXEC ) |& tee $CURDIR/buff_${{EXECUTION_ID}}.txt

  set rtarray = "${{rtarray}} `tail -3 buff_${{EXECUTION_ID}}.txt | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' | head -1` "
  rm -rf buff_${{EXECUTION_ID}}.txt

  echo
  echo "CMAQ Processing of Day $YYYYMMDD Finished at `date`"
  echo
  echo "\\\\\\=====\\\\\\=====\\\\\\=====\\\\\\=====///=====///=====///=====///"
  echo

  mv $CURDIR/CTM_LOG_???.${{CTM_APPL}} $LOGDIR
  if ( $CTM_DIAG_LVL != 0 ) then
    mv $CURDIR/CTM_DIAG_???.${{CTM_APPL}} $LOGDIR
  endif

  setenv NEW_START False

  set TODAYG = `date -ud "${{TODAYG}}+1days" +%Y-%m-%d`
  set TODAYJ = `date -ud "${{TODAYG}}" +%Y%j`

end

set RTMTOT = 0
foreach it ( `seq ${{NDAYS}}` )
    set rt = `echo ${{rtarray}} | cut -d' ' -f${{it}}`
    set RTMTOT = `echo "${{RTMTOT}} + ${{rt}}" | bc -l`
end

set RTMAVG = `echo "scale=2; ${{RTMTOT}} / ${{NDAYS}}" | bc -l`
set RTMTOT = `echo "scale=2; ${{RTMTOT}} / 1" | bc -l`

echo
echo "=================================="
echo "  ***** CMAQ TIMING REPORT *****"
echo "=================================="
echo "Start Day: ${{START_DATE}}"
echo "End Day:   ${{END_DATE}}"
echo "Number of Simulation Days: ${{NDAYS}}"
echo "Domain Name:               ${{GRID_NAME}}"
echo "Number of Grid Cells:      ${{NCELLS}}  (ROW x COL x LAY)"
echo "Number of Layers:          ${{NZ}}"
echo "Number of Processes:       ${{NPROCS}}"
echo "   All times are in seconds."
echo
echo "Num  Day        Wall Time"
set d = 0
set day = ${{START_DATE}}
foreach it ( `seq ${{NDAYS}}` )
    set d = `echo "${{d}} + 1"  | bc -l`
    set n = `printf "%02d" ${{d}}`

    set rt = `echo ${{rtarray}} | cut -d' ' -f${{it}}`

    echo "${{n}}   ${{day}}   ${{rt}}"

    set day = `date -ud "${{day}}+1days" +%Y-%m-%d`
end
echo "     Total Time = ${{RTMTOT}}"
echo "      Avg. Time = ${{RTMAVG}}"

exit""".format(proj.compiler, mn, proj.name, dom.name, dom.size, icon_type,
               bcon_type, proj.cmaq_ver, PROC, OUTDIR, proj.path.proj,
               proj.path.mcip, NEW_START, START_DATE, END_DATE, NPCOL, NPROW,
               proj.path.icon, proj.path.bcon, proj.path.emis)
    return script


def _parse_args_():
    from _helper_functions_ import _create_argparser_
    DESCRIPTION = 'icon script\n\n' + \
                  'Project: {}\n  Path: {}\nCMAQ\n  Path: {}\n  version: {}'
    DESCRIPTION = DESCRIPTION.format(proj.name, proj.path.proj,
                                     proj.path.cmaq_app, proj.cmaq_ver)
    EPILOG = 'Example of use:\n' + \
             ' %(prog)s -n 11 -y 2015 -m 2 -s 2\n'
    p = _create_argparser_(DESCRIPTION, EPILOG)
    p.add_argument('-n', '--domain', nargs='+', type=int,
                   default=proj.get_dom_ids(),
                   help="domain Id(s). Default is all domains in the project.")
    p.add_argument('-nc', '--ncores', type=int, default=-1,
                   help='Number of cores. Default is all cores.')
    p.add_argument('--ictype', default='auto',
                   choices=['auto', 'profile', 'regrid'],
                   help="default is 'auto'.")
    p.add_argument('--bctype', default='auto',
                   choices=['auto', 'profile', 'regrid'],
                   help="default is 'auto'.")
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

      if a.ncores < 1:
         a.ncores = 'auto'

    year, month, start_day, end_day = a.years, a.months, a.start_day, a.end_day
    day = [start_day] + [end_day + 1]
    doms = [proj.get_dom_by_id(i) for i in a.domain]
    ym = expandgrid(year, month)  # Year and months

    log_dir = _join(proj.path.logs, 'cctm')
    os.makedirs(log_dir, exist_ok=True)

    log.info('Creating cctm files')
    log.info('log dir: {}'.format(log_dir))

    # make sure output dir exist
    dir_out = proj.path.cctm
    os.makedirs(dir_out, exist_ok=True)
    log.info('Output Dir: {}'.format(dir_out))

    time_fmt = '[Time: {:.2f} secs]'
    flag = ExitHelper()
    cctm_timer_start = timer()
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

            script = get_script(y, m, start_day, end_day, dom, proc_type='mpi',
                                new_start=True, ncores=a.ncores,
                                icon_type=a.ictype, bcon_type=a.bctype)

            month_timer_start = timer()

            err = None
            try:
                print(script)
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

    el = month_timer_end - cctm_timer_start
    log.info('cctm completed ' + time_fmt.format(el))
