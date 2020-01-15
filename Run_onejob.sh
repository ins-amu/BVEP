#!/bin/sh
#SBATCH --ntasks 4

cwd=$(pwd)
model=$1
data_name=$2
data_input_R=$3
num_warmup=$4
num_samples=$5
delta=$6
max_depth=$7
log_file=$8

mkdir -p  data_output_CV_hmc_${model}

echo "Running HMC job Started for"  ${model} >> ${log_file}

for i in 1 2 3 4; do
    echo "running HMC started for"  $i >> ${log_file}
    ./$model id=$((100*$i))\
    sample\
    save_warmup=0 num_warmup=$num_warmup num_samples=$num_samples\
    adapt delta=$delta \
    algorithm=hmc engine=nuts max_depth=$max_depth \
    data file=${data_input_R}\
    output diagnostic_file=data_output_CV_hmc_${model}/output_hmc_diagnostic_${model}_$i.R\
    file=data_output_CV_hmc_${model}/output_hmc_${model}_$i.csv refresh=1 \
    &> data_output_CV_hmc_${model}/output_hmc_${model}_$i.out &
    echo "running HMC finished for"  $i >> ${log_file}
done

wait

echo "Running HMC job Finished for"  ${model} >> ${log_file}
