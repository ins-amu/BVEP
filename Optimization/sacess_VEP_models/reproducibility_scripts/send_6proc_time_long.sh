#!/bin/bash

#SBATCH -t 40:00:00
#SBATCH -n 6
#SBATCH --mem-per-cpu=2G

module load gcc openmpi/4.1.1_ft3 boost/1.68 gsl hdf5 python jdk

export CPPFLAGS="  -O3"
export CFLAGS=" -O3"

export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/ThirdParty-HSL/install/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/opt/libffi/include:$C_INCLUDE_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/Python-3.8.12/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/install/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/ceres-solver-2.0.0/build/install/lib64:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/build:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/ceres-solver-2.0.0/build/install:$CMAKE_PREFIX_PATH
export ParPE=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/build
export ParPE_INCLUDE_DIRS=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/include
export CMAKE_PREFIX_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/install:$CMAKE_PREFIX_PATH
export PKG_CONFIG_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/ceres-solver-2.0.0/build:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/mnt/netapp1/Optcesga_FT2_RHEL7/2020/gentoo/22072020/usr/bin/pkg-config:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/CBLAS/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/CBLAS/include:$C_INCLUDE_PATH
export HDF5_BASE=/opt/cesga/2020/software/Compiler/gcccore/system/hdf5/1.12.1
export PYTHONPATH=/home/csic/gim/dro/.local/lib/python3.9/site-packages
export PKG_CONFIG_PATH=/opt/cesga/2020/software/Compiler/gcccore/system/hdf5/1.12.1/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/ThirdParty-HSL/install/lib/pkgconfig:$PKG_CONFIG_PATH
export C_INCLUDE_PATH=/home/csic/gim/dro/.local/lib/python3.9/site-packages/amici/include:$C_INCLUDE_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/spdlog/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/install/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=/mnt/netapp1/Store_CSIC/home/csic/gim/dro/parPE/ThirdParty/Ipopt-releases-3.13.3/install/lib:$LD_RUN_PATH


echo "srun ./$1 $2 $3" 
srun ./$1 $2 $3
