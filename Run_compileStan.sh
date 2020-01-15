#!/bin/bash

cwd=$(pwd)

declare -a array=(    "BVEP_centered"   "BVEP_noncentered"    )

arraylength=${#array[@]}

echo  ..........................................................................................


for (( i=1; i<${arraylength}+1; i++ ));
do

    echo $i " / " ${arraylength} " : " ${array[$i-1]}

    model=${array[$i-1]}


    log_file=report_convergence_CV_compilestan_${model}.txt
    
    echo "Compile Started!" >> ${log_file}


    echo  ......................................................
    curr_dir=$(pwd)
    echo $(pwd) >> ${log_file}
    cd /soft/stan/cmdstan-2.17.1 && make $curr_dir/${array[$i-1]} && cd $curr_dir
    echo  ......................................................


    echo "Compile Finished for" ${alg}_${model} >> ${log_file}

done

wait
echo "Compile done!" >> ${log_file}
echo  ..........................................................................................
