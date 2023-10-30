#!/bin/bash

#SBATCH -t 04:00:00
#SBATCH -n 6
#SBATCH --mem-per-cpu=1G

module load gcc openmpi/4.1.1_ft3 python

echo "srun ./$1 $2 $3" 
srun ./$1 $2 $3
