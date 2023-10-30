#!/bin/bash

#SBATCH -t 39:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -N 1

module load gcc python


echo "srun ./$1 $2 $3" 
srun ./$1 $2 $3
