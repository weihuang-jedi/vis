#!/bin/bash

 set -x

 temp_dir=/work2/noaa/gsienkf/weihuang/jedi/run
 nodelist=(1 2 4 8)
 corelist=(36 72 144 288)

 for i in ${!nodelist[@]}
 do
   workdir=${temp_dir}/run_80.40t${nodelist[$i]}n_${corelist[$i]}p/stdoutNerr
   imgname=dist_${nodelist[$i]}n.png

   python tile-dist.py --output=1 \
	--workdir=${workdir} \
	--imagename=${imgname}
 done

