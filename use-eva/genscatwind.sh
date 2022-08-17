#!/bin/bash

 set -x

 casedir=/work2/noaa/gsienkf/weihuang/jedi/case_study
 run_dir=run_80.40t1n_36p

 caselist=(scatwind)
 var_list=(uv)
 namelist=(eastward_wind,northward_wind)

#---------------------------------------------------------------------------
 for i in ${!var_list[@]}
 do
   echo "element $i is ${var_list[$i]}"
   varname=${namelist[$i]}
   for j in ${!caselist[@]}
   do
     plotdir=${caselist[$j]}
     obsfile=${casedir}/${caselist[$j]}/${run_dir}/obsout/${caselist[$j]}_obs_2020011006_0000.nc4

     if [ -f ${obsfile} ]
     then
       if [ ! -d ${casename}/eastward_wind ]
       then
         sed -e "s?OBSFILE?${obsfile}?g" \
             -e "s?VARNAME?${varname}?g" \
             -e "s?PLOTDIR?${plotdir}?g" \
             obscoef.yaml.template > obscoef.yaml

         eva obscoef.yaml
       fi
     fi
   done
 done

