#!/bin/bash

 set -x

 casedir=/work2/noaa/gsienkf/weihuang/jedi/case_study
 run_dir=run_80.40t1n_36p

#caselist=(amsua_n19 iasi)
 caselist=(amsua)
 namelist=(brightness_temperature)

#---------------------------------------------------------------------------
 for i in ${!namelist[@]}
 do
   echo "element $i is ${namelist[$i]}"
   varname=${namelist[$i]}
   for j in ${!caselist[@]}
   do
     plotdir=${caselist[$j]}
    #obsfile=${casedir}/${caselist[$j]}/${run_dir}/obsout/${caselist[$j]}_n19_obs_2020011006_m_0000.nc4
     obsfile=${casedir}/${caselist[$j]}/${run_dir}/obsout/${caselist[$j]}_n19_obs_2020121500_m_0000.nc4

     if [ -f ${obsfile} ]
     then
      #if [ ! -d ${casename}/${varname} ]
      #then
         sed -e "s?OBSFILE?${obsfile}?g" \
             -e "s?VARNAME?${varname}?g" \
             -e "s?PLOTDIR?${plotdir}?g" \
             amsuacoef.yaml.template > amsuacoef.yaml

         eva amsuacoef.yaml
      #fi
     fi
   done
 done

