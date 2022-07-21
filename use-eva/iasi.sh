#!/bin/bash

 set -x

 casedir=/work2/noaa/gsienkf/weihuang/jedi/case_study
 run_dir=run_80.40t1n_36p

 varlist=(metop-a metop-b)
 caselist=(iasi)
 namelist=(brightness_temperature)

#---------------------------------------------------------------------------
 for k in ${!varlist[@]}
 do
 for i in ${!namelist[@]}
 do
   echo "element $i is ${namelist[$i]}"
   varname=${namelist[$i]}
   for j in ${!caselist[@]}
   do
     plotdir=${caselist[$j]}_${varlist[$k]}
     obsfile=${casedir}/${caselist[$j]}/${run_dir}/obsout/${caselist[$j]}_${varlist[$k]}_obs_2020121500_m_0000.nc4

     if [ -f ${obsfile} ]
     then
      #if [ ! -d ${casename}/${varname} ]
      #then
         sed -e "s?OBSFILE?${obsfile}?g" \
             -e "s?VARNAME?${varname}?g" \
             -e "s?PLOTDIR?${plotdir}?g" \
             iasicoef.yaml.template > iasicoef.yaml

         eva iasicoef.yaml
      #fi
     fi
   done
 done
 done

