#!/bin/bash

cwd=$(pwd)

alg="hmc"

num_warmup=2000
num_samples=2000
delta=0.99
max_depth=10

iter=50000
tol_rel_obj=0.000001


Nchain=5
Max_run=11
N_dataset=190


declare -a model_array=("Epileptor2D_Hypos")    

modelarraylength=${#model_array[@]}


for (( ii=1; ii<${modelarraylength}+1; ii++ ));
do

    echo $ii " / " ${modelarraylength} " : " ${model_array[$ii-1]}

    model=${model_array[$ii-1]}

    rm -rf report_convergence_CV*
    log_file=report_convergence_CV_${alg}_${model}.txt              
    echo "Job Started!" >> ${log_file}

    echo  ......................................................
    echo creating models ...>> ${log_file}
    cd /home/meysam/cmdstan && make $cwd/${model_array[$i-1]} && cd $cwd
    echo  compiled models ... >> ${log_file}
    echo $(pwd) >> ${log_file}
    echo  ......................................................

    mkdir -p  data_output_CV_${alg}_${model}

    echo "Running Started for" ${alg}_${model} >> ${log_file}


    for ((k = 1; k <= N_dataset; k++));
    do

                echo "running data_input2D"$k" started..." >> ${log_file}

                j=1

                while [ $j -le $Max_run ] 
                do   

                    if [ ${alg} == "opt" ]; then

                            echo "...running opt"$j" started..." >> ${log_file}

                            ./$model id=$((100*$j+$k))\
                                optimize\
                                data file=data_input2D/data_input2D$k.R\
                                output diagnostic_file=data_output_CV_${alg}_${model}/output_${alg}_diagnostic_${model}_$k-$j.R\
                                file=data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.csv refresh=1 \
                                &> data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out 
  

                            echo "...checking convergence of opt"$j"..." >> ${log_file}
                            opt_convergence_check=$(python checking_converged_${alg}.py data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out   2>&1 )
                            echo $opt_convergence_check >> ${log_file}

                            if [ "$opt_convergence_check" == "opt converged" ]; then
                                    echo "Opt converged for run"$j >> ${log_file}                                   
                                    ((j++))
                            else
                                    echo "...removing output_opt"$j"..."
                                    rm -rf data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.{out, R, csv}
                            fi
                            if [ $j -eq $Nchain ]; then
                                    echo "opt converged for all chains." >> ${log_file}
                                    break
                            fi  
                            if [ $j -eq $Max_run ]; then
                                    echo "opt arrived max iterations." >> ${log_file}
                                    break
                            fi  


                    elif [ ${alg} == "advi" ]; then
                            echo "...running ADVI"$j" started..." >> ${log_file}

                    
                            ./$model id=$((100*$j+$k))\
                                variational\
                                iter=$iter tol_rel_obj=0.000001 \
                                data file=data_input2D/data_input2D$k.R\
                                output diagnostic_file=data_output_CV_${alg}_${model}/output_${alg}_diagnostic_${model}_$k-$j.R\
                                file=data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.csv refresh=1 \
                                &> data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out

                            echo "...checking convergence of advi"$j"..." >> ${log_file}
                            advi_convergence_check=$(python checking_converged_${alg}.py data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out   2>&1 )
                            echo $advi_convergence_check >> ${log_file}

                            if [ "$advi_convergence_check" == "advi converged" ]; then
                                    echo "ADVI converged for run"$j >> ${log_file}                                   
                                    ((j++))
                            else
                                    echo "...removing output_advi"$j"..."
                                    rm -rf data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.{out, R, csv}
                            fi
                            if [ $j -eq $Nchain ]; then
                                    echo "ADVI converged for all chains." >> ${log_file}
                                    break
                            fi  
                            if [ $j -eq $Max_run ]; then
                                    echo "ADVI arrived max iterations." >> ${log_file}
                                    break
                            fi  


                    else
                            echo "...running HMC"$j" started..." >> ${log_file}

                            ./$model id=$((100*$j+$k))\
                                sample\
                                save_warmup=0 num_warmup=$num_warmup num_samples=$num_samples\
                                adapt delta=$delta \
                                algorithm=hmc engine=nuts max_depth=$max_depth \
                                data file=data_input2D/data_input2D$k.R\
                                output diagnostic_file=data_output_CV_${alg}_${model}/output_${alg}_diagnostic_${model}_$k-$j.R\
                                file=data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.csv refresh=1 \
                                &> data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out 
             

                            echo "...checking convergence of hmc"$j"..." >> ${log_file}
                            hmc_convergence_check=$(python checking_converged_${alg}.py data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.out   2>&1 )
                            echo $hmc_convergence_check >> ${log_file}

                            if [ "$hmc_convergence_check" == "hmc converged" ]; then
                                    echo "HMC converged for run"$j >> ${log_file}                                   
                                    ((j++))
                            else
                                    echo "...removing output_hmc"$j"..."
                                    rm -rf data_output_CV_${alg}_${model}/output_${alg}_${model}_$k-$j.{out, R, csv}
                            fi
                            if [ $j -eq $Nchain ]; then
                                    echo "HMC converged for all chains." >> ${log_file}
                                    break
                            fi  
                            if [ $j -eq $Max_run ]; then
                                    echo "HMC arrived max iterations." >> ${log_file}
                                    break
                            fi  


                    fi

                done

    done

    wait
    echo  ..................................................................................................

    echo "Running Finished for" ${alg}_${model} >> ${log_file}

done

wait
echo  ..........................................................................................
echo  ..........................................................................................
echo "Job done!" >> ${log_file}
