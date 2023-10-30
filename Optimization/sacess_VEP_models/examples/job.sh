# You need load the required modules, in this case the intel 
# compilers and intel mpi.

  module load intel/2.144
  module load impi/5.0.0.028
  
# It is very important putting number of the OMP variables 
# on the correspondent system variable. In this case the value
# have to be equal to 1 (only one thread per MPI processor).

  export OMP_NUM_THREADS=1

# This is the typical call to saCeSS solver with mpirun. Do not 
# change the name of variable $NSLOTS, because it loads the number
# of procesors specified in qsub call. $1 and $2 are the input 
# parameters of the script, corresponding with path of the XML
# and output files.

  mpirun -np  $NSLOTS $HOME/SACESS_TOOLS/bin/paralleltestbed  $1 $2
