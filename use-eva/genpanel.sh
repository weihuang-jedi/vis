#!/bin/bash

 set -x

 output=1

#caselist=(aircraft amsua iasi_metop-a iasi_metop-b sondes satwind scatwind vadwind windprof \
#          sfc-letkf_sfc_ps sfc-letkf_sfcship_ps sfc-letkf_sondes_ps \
#          sfc_ps sfcship_ps sondes_ps)
#namelist=(brightness_temperature air_temperature \
#          eastward_wind  northward_wind \
#          surface_pressure \
#          specific_humidity virtual_temperature)
#caselist=(sfc-letkf_sfcship_ps sfc-letkf_sfcship_ps-gsihofxbc)
#caselist=(sfc-letkf_sfc_ps-gsihofxbc sfc-letkf_sfc_ps-observer)
#namelist=(surface_pressure)
 caselist=(scatwind)
 namelist=(eastward_wind  northward_wind)

#---------------------------------------------------------------------------
 for i in ${!caselist[@]}
 do
   casename=${caselist[$i]}
   for j in ${!namelist[@]}
   do
     varname=${namelist[$j]}
     if [ -d ${casename}/${varname} ]
     then
       if [ ! -f panel-${casename}-${varname} ]
       then
         python make-panel.py --debug=1 \
           --output=${output} \
           --casename=${casename} \
           --varname=${varname} 
       fi
     fi
   done
 done

