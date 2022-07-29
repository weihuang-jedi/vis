#!/bin/bash

 set -x

 module list

#ioda-bundle build dir:
#export ioda_bld_dir=/work2/noaa/gsienkf/weihuang/ioda/ioda-bundle-bld
 export ioda_bld_dir=/work2/noaa/gsienkf/weihuang/jedi/ioda-bundle/build

 export PYTHONPATH=$PYTHONPATH:${ioda_bld_dir}/lib/python3.9/pyioda

 mkdir -p output

#usage: proc_gsi_ncdiag.py [-h] [-o OBS_DIR] [-g GEOVALS_DIR] [-d OBSDIAG_DIR] [-b ADD_OBSBIAS] [-q ADD_QCVARS] [-r ADD_TESTREFS] input_dir

#Convert to ioda_v1 data
 python ${ioda_bld_dir}/bin/proc_gsi_ncdiag.py -o output obs

 exit 0

#file name of the combined file.
 input1=sondes_q_obs_2020011006.nc4
 input2=sondes_tsen_obs_2020011006.nc4
 input3=sondes_tv_obs_2020011006.nc4
 input4=sondes_uv_obs_2020011006.nc4
 flnm=sondes_obs_2020011006.nc4

#Combine the files to a single ioda-v1 obs. file
#python ${ioda_bld_dir}/bin/combine_files.py \
 python ${ioda_bld_dir}/bin/combine_obsspace.py \
 -i ${output_dir}/$input1 \
    ${output_dir}/$input2 \
    ${output_dir}/$input3 \
    ${output_dir}/$input4 \
 -o $flnm

#Convert from ioda-v1 to ioda-v2.
 /work/noaa/gsienkf/weihuang/jedi/src/my.compile/bin/ioda-upgrade.x $flnm ioda_v2_${flnm}

