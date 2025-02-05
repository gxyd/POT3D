set -ex

MPICC=mpicc
FC=gfortran
FFLAGS="-O3 -march=native -lmpi"

cd src
${MPICC} -c mpi_wrapper.c
${FC} ${FFLAGS} -c mpi_c_bindings.f90
${FC} ${FFLAGS} -c mpi.f90
${FC} ${FFLAGS} -c psi_io.f90
${FC} -E ${FFLAGS} -cpp pot3d.F90 > pot3d_cpp.f90
${FC} -c ${FFLAGS} pot3d_cpp.f90 
${FC} ${FFLAGS} mpi_wrapper.o mpi_c_bindings.o mpi.o psi_io.o pot3d_cpp.o -o pot3d
cp pot3d ../bin/
