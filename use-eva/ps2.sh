#!/bin/bash

 set -x

 casedir=/work2/noaa/gsienkf/weihuang/jedi/case_study
 run_dir=run_80.40t1n_36p

 caselist=(surf     surf       surf)
 var_list=(sfc_ps   sfcship_ps sondes_ps)
 namelist=(surface_pressure surface_pressure surface_pressure)

#---------------------------------------------------------------------------
 for i in ${!var_list[@]}
 do
   echo "element $i is ${var_list[$i]}"
   obsfile=${casedir}/${caselist[$i]}/${run_dir}/obsout/${var_list[$i]}_obs_2020011006_0000.nc4
   varname=${namelist[$i]}
   plotdir=${var_list[$i]}

   sed -e "s?OBSFILE?${obsfile}?g" \
       -e "s?VARNAME?${varname}?g" \
       -e "s?PLOTDIR?${plotdir}?g" \
       obscoef.yaml.template > pscoef.yaml

   eva pscoef.yaml
 done

