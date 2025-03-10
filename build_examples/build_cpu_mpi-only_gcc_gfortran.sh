set -ex

if [[ "$(uname)" == "Linux" ]]; then
  CC=${CONDA_PREFIX}/bin/gcc
else
  CC=${CONDA_PREFIX}/bin/clang
fi

FC=${CONDA_PREFIX}/bin/gfortran

cd src
${CC} -I$CONDA_PREFIX/include -c mpi_wrapper.c
${FC} -O3 -march=native -c mpi_c_bindings.f90
${FC} -O3 -march=native -c mpi.f90
${FC} -O3 -march=native -c psi_io.f90
${FC} -c -O3 -march=native -cpp pot3d.F90
${FC} -O3 -march=native -lmpi mpi_wrapper.o mpi_c_bindings.o mpi.o psi_io.o pot3d.o -o pot3d
cp pot3d ../bin/
