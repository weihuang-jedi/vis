#!/bin/bash

 set -x

#ioda-bundle build dir:
 export ioda_bld_dir=/work2/noaa/gsienkf/weihuang/jedi/ioda-bundle/build

#output dir.
 output_dir=/work2/noaa/gsienkf/weihuang/jedi/case_study/diag2iodav2/out

#input file fir.
 input_dir=/work2/noaa/gsienkf/weihuang/jedi/case_study/diag2iodav2/obs
 
#Convert to ioda_v1 data
 python ${ioda_bld_dir}/bin/proc_gsi_ncdiag.py -n 1 -o ${output_dir} ${input_dir}

 exit 0

#file name of the combined file.
 flnm=sondes_obs_2021010900.nc4

#Combine the files to a single ioda-v1 obs. file
 python ${ioda_bld_dir}/bin/combine_files.py \
  	-i sondes_q_obs_2021010900.nc4  sondes_tsen_obs_2021010900.nc4  sondes_tv_obs_2021010900.nc4 \
	sondes_uv_obs_2021010900.nc4 sondes_ps_obs_2021010900.nc4 \
  	-o $flnm

#ncrename -h -O -v datetime@MetaData,LaunchTime@MetaData ${flnm}

#Convert from ioda-v1 to ioda-v2.
 /work/noaa/gsienkf/weihuang/jedi/src/my.compile/bin/ioda-upgrade.x $flnm ioda_v2_${flnm}

