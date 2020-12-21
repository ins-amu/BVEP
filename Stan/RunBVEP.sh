#!/bin/bash

# @author: meysamhashemi INS Marseille


cwd=$(pwd)

iter=50000
tol_rel_obj=0.0001

num_warmup=200
num_samples=200
delta=0.95
max_depth=10

Nchain=5
Max_opt_run=0
Max_advi_run=5
Max_hmc_run=5


declare -a data_array=("DatainputSLDecim84nodes_patient1"  "DatainputSLDecim84nodes_patient2")


dataarraylength=${#data_array[@]}

for (( i=1; i<${dataarraylength}+1; i++ ));
    do
          echo $i " / " ${dataarraylength} " : " ${data_array[$i-1]}

          data_name=${data_array[$i-1]}
          
          data_input_R=data_input_files/config2/Sourcelevel/${data_name}.R


          declare -a model_array=(   "BVEP_centered"     "BVEP_noncentered" )    

 
          modelarraylength=${#model_array[@]}



          for (( ii=1; ii<${modelarraylength}+1; ii++ ));
          do

              echo $ii " / " ${modelarraylength} " : " ${model_array[$ii-1]}

              model=${model_array[$ii-1]}

              log_file=report_convergence_benchmark_${model}_${data_name}.txt
              
              echo "Job Started!" >> ${log_file}


              echo  ......................................................
              echo creating models ...>> ${log_file}
              cd /Users/meysamhashemi/cmdstan-2.20.0 && make $cwd/${model_array[$i-1]} && cd $cwd
              echo  compiled models ... >> ${log_file}
              echo $(pwd) >> ${log_file}
              echo  ......................................................



              echo $cwd >> ${log_file}

              mkdir -p  data_output_opt_${model}_${data_name}

              echo "Running Started for" ${model}_${data_name} >> ${log_file}

              j=1

              while [ $j -le $Max_opt_run ] 
              do
                  
                                    echo "...running opt_"${model}_${data_name}_$j"..." >> ${log_file}

                                    ./$model id=$((100*$j + $ii))\
                                    optimize \
                                    data file=${data_input_R}  \
                                    output file=data_output_opt_${model}_${data_name}/output_opt_${model}_${data_name}_$j.R refresh=1 \
                                    &> data_output_opt_${model}_${data_name}/output_opt_${model}_${data_name}_$j.out
                                    

                                    echo "...checking convergence of opt_"${model}_${data_name}_$j"..." >> ${log_file}
                                    opt_convergence_check=$(python checking_converged_opt.py data_output_opt_${model}_${data_name}/output_opt_${model}_${data_name}_$j.out 2>&1 )

                                    echo $opt_convergence_check  >> ${log_file}

                                    if [ "$opt_convergence_check" == "opt converged" ]; then
                                        echo "converged output_opt_"${model}_${data_name}_$j >> ${log_file}
                                        ((j++))
                                    else
                                        echo "removed output_opt_"${model}_${data_name}_$j >> ${log_file}
                                        rm -rf data_output_opt_${model}_${data_name}/output_opt_${model}_${data_name}_$j.{out,R,csv}
                                    fi


                                    if [ $j -eq $Nchain ]; then
                                        echo "Opt converged for all chains." >> ${log_file}
                                        break
                                    fi  

                                    if [ $j -eq $Max_opt_run ]; then
                                        echo "Opt arrived max iterations." >> ${log_file}
                                        break
                                    fi  

              done

              wait
              echo  ..........................................................................................
               
               
              mkdir -p  data_output_advi_${model}_${data_name}
           
              jj=1

              while [ $jj -le $Max_advi_run ] 

              do
                  
                                    echo "...running advi_"${model}_${data_name}_$jj"..." >> ${log_file}

                                    ./$model id=$((100*$jj + $ii))\
                                    variational\
                                    iter=$iter tol_rel_obj=0.0001 \
                                    data file=${data_input_R}\
                                    output diagnostic_file=data_output_advi_${model}_${data_name}/output_advi_diagnostic_${model}_${data_name}_$jj.R\
                                    file=data_output_advi_${model}_${data_name}/output_advi_${model}_${data_name}_$jj.csv refresh=1 \
                                    &> data_output_advi_${model}_${data_name}/output_advi_${model}_${data_name}_$jj.out


                                    echo "...checking convergence of advi_"${model}_${data_name}_$jj"..." >> ${log_file}
                                    advi_convergence_check=$(python checking_converged_advi.py data_output_advi_${model}_${data_name}/output_advi_${model}_${data_name}_$jj.out 2>&1 )

                                    echo $advi_convergence_check >> ${log_file}

                                    if [ "$advi_convergence_check" == "advi converged" ]; then
                                        echo "converged output_advi_"${model}_${data_name}_$jj >> ${log_file}
                                        ((jj++))
                                    else
                                        echo "removed output_advi_"${model}_${data_name}_$jj >> ${log_file}
                                        rm -rf data_output_advi_${model}_${data_name}/output_advi_${model}_${data_name}_$jj.{out,R,csv}
                                    fi


                                    if [ $jj -eq $Nchain ]; then
                                        echo "ADVI converged for all chains." >> ${log_file}
                                        break
                                    fi  

                                    if [ $jj -eq $Max_advi_run ]; then
                                        echo "ADVI arrived max iterations." >> ${log_file}
                                        break
                                    fi  

              done

              wait
              echo  ..........................................................................................
               

              mkdir -p  data_output_hmc_${model}_${data_name}

              jjj=1

              while [ $jjj -le $Max_hmc_run ] 

              do
                  
                                    echo "...running hmc_"${model}_${data_name}_$jjj"..." >> ${log_file}

                                    ./$model id=$((100*$jjj + $ii))\
                                    sample\
                                    save_warmup=0 num_warmup=$num_warmup num_samples=$num_samples\
                                    adapt delta=$delta \
                                    algorithm=hmc engine=nuts max_depth=$max_depth \
                                    data file=${data_input_R}\
                                    output diagnostic_file=data_output_hmc_${model}_${data_name}/output_hmc_diagnostic_${model}_${data_name}_$jjj.R\
                                    file=data_output_hmc_${model}_${data_name}/output_hmc_${model}_${data_name}_$jjj.csv refresh=1 \
                                    &> data_output_hmc_${model}_${data_name}/output_hmc_${model}_${data_name}_$jjj.out


                                    echo "...checking convergence of hmc_"${model}_${data_name}_$jjj"..." >> ${log_file}
                                    hmc_convergence_check=$(python checking_converged_hmc.py data_output_hmc_${model}_${data_name}/output_hmc_${model}_${data_name}_$jjj.out 2>&1 )

                                    echo $hmc_convergence_check >> ${log_file}

                                    if [ "$hmc_convergence_check" == "hmc converged" ]; then
                                        echo "converged output_hmc_"${model}_${data_name}_$jjj >> ${log_file}
                                        ((jjj++))
                                    else
                                        echo "removed output_hmc_"${model}_${data_name}_$jjj >> ${log_file}
                                        rm -rf data_output_hmc_${model}_${data_name}/output_hmc_${model}_${data_name}_$jjj.{out,R,csv}
                                    fi


                                    if [ $jjj -eq $Nchain ]; then
                                        echo "HMC converged for all chains." >> ${log_file}
                                        break
                                    fi  

                                    if [ $jjj -eq $Max_hmc_run ]; then
                                        echo "HMC arrived max iterations." >> ${log_file}
                                        break
                                    fi  

              done

              wait
              echo  ..........................................................................................
            
              echo "Running Finished for"  ${model}_${data_name} >> ${log_file}

          done

done
wait
echo  ..........................................................................................
echo  ..........................................................................................
echo "Job done!" >> ${log_file}