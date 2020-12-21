#!/bin/bash

cwd=$(pwd)


iter=50000
tol_rel_obj=0.0001

num_warmup=200
num_samples=200
delta=0.95
max_depth=10

declare -a data_array=("DatainputSLDecim84nodes_patient1"  "DatainputSLDecim84nodes_patient2" )




declare -a model_array=(    "BVEP_centered"    "BVEP_noncentered"  )

dataarraylength=${#data_array[@]}
modelarraylength=${#model_array[@]}

echo  ..........................................................................................
echo "Job Started!" >> ${log_file}
echo  ..........................................................................................

for (( i=1; i<${dataarraylength}+1; i++ ));
    do
          echo $i " / " ${dataarraylength} " : " ${data_array[$i-1]}

          data_name=${data_array[$i-1]}
          data_input_R=data_input_files/config2/Sourcelevel/${data_name}.R

          for (( ii=1; ii<${modelarraylength}+1; ii++ ));
          do

              echo $ii " / " ${modelarraylength} " : " ${model_array[$ii-1]}

              model=${model_array[$ii-1]}

              log_file=report_convergence_CV_${model}.txt
              
              echo $cwd >> ${log_file}   

              echo "...submitting hmc_"${model}"..." >> ${log_file}

              sbatch  Run_onejob.sh $model $data_name $data_input_R $num_warmup $num_samples $delta $max_depth $log_file &

              echo "...terminated hmc_"${model}"..." >> ${log_file}

          done

done
echo  ..........................................................................................
echo "Job Done!" >> ${log_file}
echo  ..........................................................................................
