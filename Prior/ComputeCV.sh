#!/bin/bash
model=Epileptor2D_Hypos
alg=hmc

python ComputeWAIC_Hypos.py data_input2D    ${model}  data_output_CV_${alg}_${model}

python ComputePSIS_Hypos.py data_input2D  ${model}  data_output_CV_${alg}_${model}