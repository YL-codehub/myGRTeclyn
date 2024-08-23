#!/bin/bash
# To load the modules below call this file as an argument to `source`

# Check that we are running on an Icelake (either login or compute) node
if [[ $(hostname) != *"-q-"* ]]; then
	echo "These modules will only work on an Icelake (Rocky Linux 8) node." 1>&2
	echo "See https://docs.hpc.cam.ac.uk/hpc/user-guide/connecting.html"
	echo "Not loading anything..." 1>&2
	return 1
fi




module purge
module load rhel8/default-icl
module load intel-oneapi-compilers/2022.1.0/gcc/b6zld2mz
module load intel-oneapi-mpi/2021.6.0/intel/guxuvcpm
module load intel-oneapi-mkl/2022.1.0/intel/intel-oneapi-mpi/qqwrrcxw 
module load fftw/3.3.10/intel/intel-oneapi-mpi/uqgpuef7

# str1="/../../../Chombo/lib/" 
# newdir="$PWD$str1"
# export CHOMBO_HOME=$newdir
# Uncomment if using AHFinder
#module load petsc/3.17-icl

# Uncomment if using TwoPunctures
#module load gsl/2.7.1/gcc/cjb5f4ui