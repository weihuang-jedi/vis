#!/bin/bash

 datadir=/work2/noaa/gsienkf/weihuang/jedi/run/run_80.40t4n_156p/analysis/increment
 obsfile=/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/ioda_v2_data/sfc_ps_obs_2020011006.nc4

 output=0

 gridfile=regrid/fv3latlon.nc
 if [ ! -f ${gridfile} ]
 then
   cd regrid
  #interp.sh ${case_dir}/sfc-letkf/analysis/increment/
   interp.sh ${datadir}/
   cd ..
 fi

 python all-obs-jedi.py \
   --output=${output} \
   --casename=AllObs \
   --varname=T \
   --gridfile=${gridfile} \
   --obsfile=${obsfile}

