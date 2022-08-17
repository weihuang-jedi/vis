#!/bin/bash

 set -x

 casedir=/work2/noaa/gsienkf/weihuang/jedi/case_study

 varname=surface_pressure
 titlelist=(ship_observer ship_gsihofxbc sfc_observer sfc_gsihofxbc)
 obslist=(/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/obsout.after-observer/sfcship_ps_obs_2020011006_0000.nc4\
          /work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/obsout.gsihofxbc/sfcship_ps_obs_2020011006_0000.nc4\
          /work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/obsout.after_observer.sfc_ps/sfc_ps_obs_2020011006_0000.nc4\
          /work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/obsout.sfc_ps.gsihofxbc/sfc_ps_obs_2020011006_0000.nc4)

#---------------------------------------------------------------------------
 for i in ${!obslist[@]}
 do
   echo "element $i is ${obslist[$i]}"
   obsfile=${obslist[$i]}
   plotdir=${titlelist[$i]}
   casename=${titlelist[$i]}

   sed -e "s?OBSFILE?${obsfile}?g" \
       -e "s?VARNAME?${varname}?g" \
       -e "s?PLOTDIR?${plotdir}?g" \
       -e "s?CASENAME?${casename}?g" \
       obscoef.yaml.template > pscoef.yaml

   eva pscoef.yaml
 done

