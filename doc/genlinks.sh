#!/bin/bash

 set -x

 src_dir=/work2/noaa/gsienkf/weihuang/C96_psonly_delp/2020011006

 number_members=80
 n=1
 while [ $n -le $number_members ]
 do
   if [ $n -lt 10 ]
   then
     member_str=mem00${n}
   elif [ $n -lt 100 ]
   then
     member_str=mem0${n}
   else
     member_str=mem${n}
   fi

  #ln -sf ${src_dir}/${member_str}/INPUT ${member_str}
   cp coupler.res ${member_str}/.

   n=$(( $n + 1 ))
 done

