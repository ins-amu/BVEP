BUILDDIR=.

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
export ROOT_DIR

LIBRARY=$(BUILDDIR)/lib
AMIGO_PATH=$(LIBRARY)/libAMIGO

PYTHON_INC=/usr/include/python3.10
#############################################################################
## PARALLEL GCC/GFORTRAN - OPENMP
#############################################################################
CC:=mpicc 
FC:=mpif90 
CLIBS+=  -lstdc++ -lpthread -lrt -lgfortran  -cpp -MMD -lm -ldl -lz 
CFLAGS+= -O3 -cpp -DGNU -fPIC -no-pie 
CPARALLEL+=-DOPENMP -DMPI2 -fopenmp -lmpi
FLIBS+=
FFLAGS+= -O3 -cpp -DGNU  -fallow-argument-mismatch -std=gnu 
FPARALLEL+=-DOPENMP -DMPI2 -DGNU -fopenmp

LIBS+= -L$(LIBRARY)/BLAS -lblas  -Dnullptr=0 -lpython3.10
INC+=  -I$(LIBRARY)/misqp/gnu  -I$(PYTHON_INC)
LIBS+= -L$(LIBRARY)/misqp/gnu
#MISQP= -L$(LIBRARY)/misqp/gnu -lmisqp -I$(LIBRARY)/misqp/gnu -lgfortran  $(LIBRARY)/misqp/gnu/libmisqp.so
MISQP=$(LIBRARY)/misqp/gnu/libmisqp.so
FLIBS+= $(MISQP)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/python/*.c)
INC+=-I$(BUILDDIR)/benchmarks/python
#############################################################################

#############################################################################
## PARALLEL ICC/IFORT - INTELMPI
############################################################################
#CC:=   mpiicc 
#FC:=   mpiifort 
#AR:=   xiar#
#LD:=   xild#
#CLIBS+=  -cxxlib -lrt -lhdf5  -limf -lifcore 
#CLIBS+=  -lm -MMD -lz 
#CFLAGS+= -O3 -ipo -xHost -DINTEL 
#CPARALLEL+= -openmp -pthread -mt_mpi  -DOPENMP -DMPI2 
#FLIBS+=  -cxxlib 
#FFLAGS+= -O3 -ipo -xHost -fpp -DEXPORT   
#FPARALLEL+= -openmp -DOPENMP -DMPI2 -DINTEL -mt_mpi -Dnullptr=0 -lpython3.9
#BLAS+= $(LIBRARY)/BLAS/libblas.a 


#LIBS+= -L$(LIBRARY)/BLAS -lblas
#INC+=  -I$(LIBRARY)/misqp/intel -I$(PYTHON_INC)
#LIBS+= -L$(LIBRARY)/misqp/intel
#MISQP=   $(LIBRARY)/misqp/intel/libmisqp.so
#############################################################################

#############################################################################
## PARALLEL ICC/IFORT - OPENMPI
#############################################################################
#OPENMPI:=/opt/MPI/openmpi-1.10.2/linux/intel_16.0.2.181
#CC:=   mpicc
#FC:=   mpif90
#AR:=   xiar#
#LD:=   xild#
#CLIBS+=  -cxxlib -lrt -lhdf5  -lxml2 -limf -lifcore
#CLIBS+=  -lm -MMD -lz
#CFLAGS+= -O3 -ipo -xHost -DINTEL
#CPARALLEL+= -openmp -pthread -lmpi -DOPENMP -DMPI2 
#FLIBS+=  -cxxlib
#FFLAGS+= -O3 -ipo -xHost -fpp -DEXPORT -DINTEL 
#FPARALLEL+= -openmp -DOPENMP -DMPI2 -lmpi 
#BLAS+= $(LIBRARY)/BLAS/libblas.a
#LIBS+= -L$(LIBRARY)/BLAS -lblas
#INC+=  -I$(LIBRARY)/misqp/intel
#LIBS+= -L$(LIBRARY)/misqp/intel
#MISQP=   $(LIBRARY)/misqp/intel/libmisqp.so
##############################################################################

#############################################################################
## SEQUENTIAL GCC/GFORTRAN
#############################################################################
#CC:= gcc 
#FC:= gfortran
#CLIBS+= -lstdc++ -lpthread -lrt -lgfortran  -cpp -MMD -lm -ldl -lz
#CFLAGS+= -O3 -cpp -DGNU  -fPIC -no-pie 
#FLIBS+=
#FFLAGS+= -O3 -cpp -DGNU   -fallow-argument-mismatch -std=gnu 
#INC+=  -I$(LIBRARY)/misqp/gnu
#LIBS+= -L$(LIBRARY)/misqp/gnu
#MISQP=   $(LIBRARY)/misqp/gnu/libmisqp.so

#LIBS+= -L$(LIBRARY)/BLAS -lblas  -Dnullptr=0 -lpython3.9
#INC+=  -I$(LIBRARY)/misqp/gnu  -I$(PYTHON_INC)

#SRC+=$(wildcard $(BUILDDIR)/benchmarks/python/*.c)
#INC+=-I$(BUILDDIR)/benchmarks/python
#############################################################################


#############################################################################
## SEQUENTIAL ICC/IFORT 
##############################################################################
#CC:=   icc
#FC:=   ifort
#AR:=   xiar#
#LD:=   xild#
#CLIBS+=  -cxxlib -lrt -lhdf5  -lxml2 -limf -lifcore  -lm -MMD -lz 
#CFLAGS+= -O3 -ipo -xHost -DINTEL
#FLIBS+=  -cxxlib 
#FFLAGS+= -O3 -ipo -xHost -fpp -DEXPORT -DINTEL 
#BLAS+= $(LIBRARY)/BLAS/libblas.a  
#LIBS+= -L$(LIBRARY)/BLAS -lblas
#INC+=  -I$(LIBRARY)/misqp/intel
#LIBS+= -L$(LIBRARY)/misqp/intel
#MISQP=   $(LIBRARY)/misqp/intel/libmisqp.so
#############################################################################

#############################################################################
## XML
#############################################################################
INC+= -I$(LIBRARY)/libxml2/include
LIBS+= -L$(LIBRARY)/libxml2/lib -lxml2 
#############################################################################

#############################################################################
## HDF5
#############################################################################
INC+=-I$(LIBRARY)/hdf5-1.8.12/include
LIBS+=-L$(LIBRARY)/hdf5-1.8.12/lib  -lhdf5 
#############################################################################

#############################################################################
## MATLAB
#############################################################################
#MATLAB_HOME = /state/partition1/apps/matlab/R2009_x64
#INC+=-I$(MATLAB_HOME)/extern/include
#LIBS+=  -L$(MATLAB_HOME)/bin/glnxa64 -L$(MATLAB_HOME)/extern/include -leng -lmat -lmex
#CFLAGS+= -DMATLAB
#FFLAGS+= -DMATLAB
#SRC+=$(wildcard $(BUILDDIR)/benchmarks/matlab/*.c)
#INC+=-I$(BUILDDIR)/benchmarks/matlab
#SRCFORTRAN_MATLAB=$(wildcard $(BUILDDIR)/benchmarks/matlab/matlabobjfunc.f90)
#############################################################################

#############################################################################
## R
#############################################################################
R_HOME = /usr/share/R
INC+=-I$(R_HOME)/include
LIBS+=-L$(R_HOME)/bin -L$(R_HOME)/lib -lR #-lRblas
#############################################################################

export AR
export CC
export FC
export CLIBS
export CFLAGS
export FLIBS
export FFLAGS
SRC+=$(wildcard $(BUILDDIR)/benchmarks/bbob/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/customized/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B1/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B2/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B3/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B4/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B5/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/others/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/others/3-step_pathway/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/others/circadian/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/others/Nfkb/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/MINLP/SSP/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/MINLP/HePG2/*.c)
SRC+=$(wildcard $(BUILDDIR)/benchmarks/systemsBiology/MINLP/HPN-DREAM/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/input_module/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/output/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/structure_data/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/method_module/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/method_module/DE/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/method_module/eSS/*.c)
SRC+=$(wildcard $(BUILDDIR)/src/*.c)

SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/scattersearchtypes.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/common_functions.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/funcevalinterface.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/qsort_mod.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/outputhdf5.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/dhc/dhc.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/misqp/misqp_interface.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/localsolver.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/scattersearchfunctions.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/localsolverinterfacec.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/parallelscattersearchfunctions.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/scattersearch.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/cess.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/sacess.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/acessdist.f90)
SRCFORTRAN+=$(wildcard $(BUILDDIR)/src/method_module_fortran/eSS/essm.f90)
SRCFORTRAN+=$(SRCFORTRAN_MATLAB)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_cvodes/*.c)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_amigo/*.c)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_de/*.c)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_fullC/*.c)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_julia/*.c)
SRCAMIGO+=$(wildcard $(AMIGO_PATH)/src/src_SRES/*.c)

INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B1/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/fly/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/util/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B5/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B4/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B3/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B2/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/others/circadian/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/others/Nfkb/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/others/3-step_pathway/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/MINLP/HePG2/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/MINLP/SSP/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/MINLP/HPN-DREAM/
INC+=-I$(BUILDDIR)/benchmarks/systemsBiology/
INC+=-I$(BUILDDIR)/benchmarks/customized/
INC+=-I$(BUILDDIR)/benchmarks/bbob/
INC+=-I$(AMIGO_PATH)/include/include_cvodes/
INC+=-I$(AMIGO_PATH)/include/include_amigo/
INC+=-I$(AMIGO_PATH)/include/include_de/
INC+=-I$(AMIGO_PATH)/include/include_SRES/
INC+=-I$(BUILDDIR)/include/
INC+=-I$(BUILDDIR)/include/input/
INC+=-I$(BUILDDIR)/include/method_module/
INC+=-I$(BUILDDIR)/include/method_module/eSS/
INC+=-I$(BUILDDIR)/include/method_module/DE/
INC+=-I$(BUILDDIR)/include/bbob/
INC+=-I$(BUILDDIR)/include/structure_data/
INC+=-I$(BUILDDIR)/include/output/
INC+=-I$(BUILDDIR)/include/error/
INC+=-I$(BUILDDIR)/include/method_module_fortran/
INC+=-I$(BUILDDIR)/src/method_module_fortran/eSS/
INC+=-I$(BUILDDIR)/src/method_module_fortran/eSS/localsolvers
INC+=-I$(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/dhc
INC+=-I$(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/misqp
INC+=-I$(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/nl2sol

INC+=  -I$(LIBRARY)/gsl-1.14 -I$(LIBRARY)/gsl-1.14/include  -I$(LIBRARY)/gsl-1.14/include/gsl
LIBS+= -L$(LIBRARY)/gsl-1.14/lib -lgsl
LIBS+= -lAMIGO -fPIC -DEXPORT
LIBS+= -L$(AMIGO_PATH)/lib

PROG := bin/paralleltestbed

OBJFILES := $(SRC:.c=.o)
DEPFILES := $(SRC:.c=.d)

OBJN2SOL=$(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/nl2sol/src_nl2sol/*.o
OBJB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/util/*.o
OBJB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/fly/*.o
OBJB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/*.o
DEPB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/util/*.d
DEPB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/fly/*.d
DEPB6FILES+=$(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6/*.d

OBJFORTRANFILES := $(SRCFORTRAN:.f90=.o)
DEPFORTRANFILES := $(SRCFORTRAN:.f90=.d)
MODFORTRANFILES := $(SRCFORTRAN:.f90=.mod)

OBJAMIGO := $(SRCAMIGO:.c=.o)
DEPAMIGO := $(SRCAMIGO:.c=.d)	


all: BLAS B6 N2SOL libAMIGO.a  $(OBJFORTRANFILES)  $(OBJCPPFILES) $(PROG)

BLAS:
	cd $(BUILDDIR)/lib/BLAS && $(MAKE) 

B6:
	cd $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6 && $(MAKE) deps
	cd $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6 && $(MAKE)

N2SOL:
	cd $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/nl2sol/src_nl2sol && $(MAKE) 	

libAMIGO.a :  $(OBJAMIGO)  
	$(AR) r $(AMIGO_PATH)/lib/libAMIGO.a $^ $(OBJN2SOLFILES)

%.o: %.c
	$(CC)  -c $< -o $@  $(LIBS) $(INC) $(CLIBS) $(CFLAGS) $(CPARALLEL)

%.o: %.f90
	$(FC) -c $< -o $@ $(LIBS) $(INC) $(FLIBS) $(FFLAGS) $(FPARALLEL)
	@mv $(BUILDDIR)/*.mod $(BUILDDIR)/include/method_module_fortran/
	
$(PROG) : $(OBJFILES)
	$(LINK.o) -o $(PROG) $(MISQP) $(INC) $(BLAS) $(OBJN2SOL) $(OBJB6FILES) $(OBJFILES) $(OBJFORTRANFILES) $(CFLAGS) $(LIBS) $(CLIBS) $(FLIBS) $(CPARALLEL) 

clean :
	rm -f $(PROG) $(OBJFILES) $(DEPFILES)
	rm -f $(OBJFORTRANFILES) $(DEPFORTRANFILES) $(MODFORTRANFILES)

veryclean :
	rm -f $(PROG) $(OBJFILES) $(DEPFILES)
	rm -f $(OBJFORTRANFILES) $(DEPFORTRANFILES) $(MODFORTRANFILES)
	rm -f $(AMIGO_PATH)/lib/libAMIGO.a $(OBJAMIGO) $(DEPAMIGO)
	rm -f $(BUILDDIR)/include/method_module_fortran/*.mod
	cd $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6 && $(MAKE) clean
	cd $(BUILDDIR)/benchmarks/systemsBiology/BioPredyn/B6 && $(MAKE) veryclean
	cd $(BUILDDIR)/src/method_module_fortran/eSS/localsolvers/nl2sol/src_nl2sol && $(MAKE) clean
	cd $(BUILDDIR)/lib/BLAS && make clean
