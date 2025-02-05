# build script which uses LFortran as it's compiler
# along with using our own MPI wrappers and call C MPI instead
# and doesn't use any HDF5
POT3D_HOME=$PWD
MPICC=${CONDA_PREFIX}/bin/mpicc
FC="$PWD/../src/bin/lfortran"
FFLAGS="-lmpi"

echo "==> Using Fortran compiler: " ${FC}
echo "=== Starting build process ==="

echo "Entering src directory"
cd src
${MPICC} -c mpi_wrapper.c
${FC} -c ${FFLAGS} mpi_c_bindings.f90
${FC} -c ${FFLAGS} mpi.f90
${FC} -c ${FFLAGS} psi_io.f90
${FC} -E ${FFLAGS} --cpp pot3d.F90 > pot3d_cpp.f90
${FC} -c ${FFLAGS} --implicit-interface pot3d_cpp.f90 --show-asr
