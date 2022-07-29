#!/bin/bash

 set -x

 export PYTHONPATH=$PYTHONPATH:/work2/noaa/gsienkf/weihuang/ioda/ioda-bundle-bld/lib/python3.9/pyioda

#ioda-bundle build dir:
 export ioda_bld_dir=/work2/noaa/gsienkf/weihuang/jedi/ioda-bundle/build

#data dir
 data_dir=/work2/noaa/gsienkf/weihuang/jedi/case_study/jeff-GSI-data

 mkdir -p obs
 rm -f obs/*

 total_members=80
 n=80
 while [ $n -le $total_members ]
 do
   n=$(( $n + 1 ))
   if [ $n -lt 10 ]
   then
     member_str=mem00${n}
   elif [ $n -lt 100 ]
   then
     member_str=mem0${n}
   else
     member_str=mem${n}
   fi

   if [ $n -gt $total_members ]
   then
     member_str=ensmean
   fi

   obsfile=diag_conv_ps_ges.2020011006_${member_str}.nc4

   mkdir -p obs output

   cd obs
   ln -sf ${data_dir}/${obsfile} .
   cd ..

   rm -rf ${member_str}

  #Convert to ioda v2 data
   python ${ioda_bld_dir}/bin/proc_gsi_ncdiag.py -o output obs

   mv output ${member_str}

   rm obs/${obsfile}
 done

