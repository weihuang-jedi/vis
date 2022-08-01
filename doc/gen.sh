#!/bin/bash

#set -x

 taskspernode=40
 NUMMEM=80
 MYLAYOUT=1,20

#cpu_list=(   6    12    18    24    30    36     72     78)
#nodelist=(   1     1     1     1     1     1      2      2)
#lay_list=("1,1" "1,2" "1,3" "1,4" "1,5" "1,6" "1,12" "1,13")

#cpu_list=(  36)
#nodelist=(   1)
#lay_list=("2,3")

#cpu_list=(  78)
#nodelist=(   2)
#lay_list=(1,13)

#cpu_list=( 156)
#nodelist=(   4)
#lay_list=(2,13)

#cpu_list=( 312)
#nodelist=(   8)
#lay_list=(4,13)

 cpu_list=( 36   78  156  312)
 nodelist=(  1    2    4    8)
 lay_list=(3,2 1,13 2,13 4,13)

 templatedir=/work2/noaa/gsienkf/weihuang/jedi/case_study/templates
 topdir=/work2/noaa/gsienkf/weihuang/jedi/case_study
#caselist=(aircraft amsua iasi satwind scatwind sfcship sondes surf vadwind windprof)
 caselist=(amsua iasi)

#------------------------------------------------------------------------------
 n=0
 for j in ${!caselist[@]}
 do
   case=${caselist[$j]}
   echo "Case: $case"
   casedir=${topdir}/${case}

   cd ${casedir}
   ln -sf ../Data .
   ln -sf ../ioda_v2_data .

   for i in ${!cpu_list[@]}
   do
     echo "element $i is ${myArray[$i]}"
     totalcpus=${cpu_list[$i]}
     nodes=${nodelist[$i]}
     MYLAYOUT=${lay_list[$i]}

     workdir=${casedir}/run_${NUMMEM}.${taskspernode}t${nodes}n_${totalcpus}p
     mkdir -p ${workdir}
     cd ${workdir}

     sed -e "s?TASKSPERNODE?${taskspernode}?g" \
         -e "s?TOTALNODES?${nodes}?g" \
         -e "s?TOTALCPUS?${totalcpus}?g" \
         -e "s?WORKDIR?${workdir}?g" \
         -e "s?NUMMEM?${NUMMEM}?g" \
         -e "s?MYLAYOUT?${MYLAYOUT}?g" \
         ${templatedir}/slurm.template > run.slurm

     sed -e "s?LAYOUT?${MYLAYOUT}?" \
         -e "s?NUMBEROFMEMBERS?${NUMMEM}?" \
         ${templatedir}/getkf.yaml.template > getkf.yaml

     cat ${templatedir}/${case}.obs.yaml.template >> getkf.yaml

    #if [ $n -lt 1 ]
    #then
      #DEPEND=$(sbatch --parsable run.slurm)
      #echo "job_id: ${DEPEND}"
      #sbatch run.slurm
    #else
      #DEPEND=$(sbatch --dependency=afterany:${DEPEND} --parsable run.slurm)
       sbatch run.slurm
    #fi

     n=$((n+1))
     cd ${casedir}
   done
 done

