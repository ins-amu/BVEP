#!/bin/bash

#SBATCH --nodelist=c[1,3-4]
#SBATCH --mem=120000

echo  ..........................................................................................

cwd=$(pwd)

model="hmc"

log_file=report_convergence_${model}.txt

echo  ..........................................................................................
echo "Job starts!" >> ${log_file}

python3 BVEP_ode_sbi_seeg_GrExp_patient1.py

echo "Job done!" >> ${log_file}
echo  ..........................................................................................
