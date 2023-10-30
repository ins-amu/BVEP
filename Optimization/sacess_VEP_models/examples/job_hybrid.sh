# You need load the required modules, in this case the intel 
# compilers and intel mpi.

  module load intel/2.144
  module load impi/5.0.0.028
  
# It is very important putting number of the OMP variables 
# on the correspondent system variable. In this case the value
# have to be equal to 1 (only one thread per MPI processor).

  export OMP_NUM_THREADS=3

# This is the typical call to saCeSS solver with mpirun for hybrid jobs. 
# Do not  use the variable $NSLOTS, because now the number of slots asked
# by qsub command (mpi+openmp) is different than the number of MPI processes.
# The parameter --map-by slot:pe=$OMP_NUM_THREADS means that for each mpi
# process, the job needs a number of slots equal to the $OMP_NUM_THREADS 
# value. We need to put this parameter for right allocation of the jobs in
# the resources, avoiding the oversubscription. 
# $1 and $2 are the input parameters of the script, corresponding
#  with path of the XML and output files.

  mpirun  --map-by slot:pe=$OMP_NUM_THREADS -np  11 $HOME/SACESS_TOOLS/bin/paralleltestbed  $1 $2


