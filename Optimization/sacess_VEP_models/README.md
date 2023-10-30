# saCeSS global optimization library + VEP models  #

Here is a version of the Self-Adaptive Cooperative Enhanced Scatter Search (saCeSS) that integrates a set of calibration problems for Virtual Epileptic Patient (VEP) models. For further details about the implementation or installation of saCeSS, please refer to the original repository:

https://bitbucket.org/DavidPenas/sacess-library

The saCeSS code has been implemented using Fortran 90 and C. Parallelization has been achieved using MPI and OpenMP. It has been tested on Linux clusters running Rocky Linux release 8.4 and Ubuntu 22.04.1 LTS.

The saCeSS library facilitates the solution of Nonlinear Programming (NLP) and Mixed-Integer Nonlinear Programming (MINLP) problems. It also provides efficient local solvers for nonlinear parameter estimation problems associated with complex models (e.g., those described by differential equations).

The input files for VEP models are located in the "inputs" folder, with the following filenames:
1. Python_P1.xml -> Forward simulation of VEP model with weak coupling at source level.
2. Python_P2.xml -> Forward simulation of VEP model with weak coupling at sensor level.
3. Python_P3.xml -> Forward simulation of VEP model with strong coupling at source level.
4. Python_P4.xml -> Forward simulation of VEP model with strong coupling at sensor level.
5. Python_P5.xml -> Forward simulation of SDE VEP model with large tau (stiff equations) at source level.
6. Python_P6.xml -> Forward simulation of SDE VEP model with large tau (stiff equations) at sensor level.


# Reproducibility 

Scripts in the "reproducibility_scripts" folder can be used to reproduce the results on a cluster machine with a queue system. 


## Example PROBLEM 1
### Sequential scripts (compile SaCeSS without MPI)

sbatch ./send_1proc_time_short.sh  ./bin/paralleltestbed  inputs/Python_P1_SEQ.xml output/P1_output_1proc_run01


### Parallel scripts (compile SaCeSS using MPI)

### 6 processors
sbatch ./send_6proc_time_short.sh  ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_6proc_run01

sbatch ./send_6proc_time_short.sh  ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_6proc_run02

...

sbatch ./send_6proc_time_short.sh  ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_6proc_run10


### 12 processors
sbatch ./send_12proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_12proc_run01

sbatch ./send_12proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_12proc_run02

...

sbatch ./send_12proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_12proc_run10


### 24 processors
sbatch ./send_24proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_24proc_run01

sbatch ./send_24proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_24proc_run02

...

sbatch ./send_24proc_time_short.sh ./bin/paralleltestbed  inputs/Python_P1.xml output/P1_output_24proc_run10

## Example PROBLEM 4

### Sequential scripts (compile SaCeSS without MPI)
sbatch ./send_1proc_time_long.sh  ./bin/paralleltestbed  inputs/Python_P4_SEQ.xml output/P4_output_1proc_run01


### Parallel scripts (compile SaCeSS using MPI)

### 6 processors
sbatch ./send_6proc_time_long.sh  ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_6proc_run01

sbatch ./send_6proc_time_long.sh  ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_6proc_run02

...

sbatch ./send_6proc_time_long.sh  ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_6proc_run10


#### 12 processors
sbatch ./send_12proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_12proc_run01

sbatch ./send_12proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_12proc_run02

...

sbatch ./send_12proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_12proc_run10

#### 24 processors
sbatch ./send_24proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_24proc_run01

sbatch ./send_24proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_24proc_run02

...

sbatch ./send_24proc_time_long.sh ./bin/paralleltestbed  inputs/Python_P4.xml output/P4_output_24proc_run10


## REFERENCES ##

### SaCeSS reference: ###

Penas, D.R., P. Gonzalez, J.A. Egea, R. Doallo and J.R. Banga (2017) Parameter estimation in large-scale systems biology models: a parallel and self-adaptive cooperative strategy. BMC Bioinformatics 18:52.

