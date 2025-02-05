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
# even if the commands fail, they don't necessarily exit with non-zero status
if [ ! -f mpi_wrapper.o ]; then
    echo "Error: mpi_wrapper.o not created!" >&2
    exit 1
fi

${FC} -c ${FFLAGS} mpi_c_bindings.f90
if [ ! -f mpi_c_bindings.o ]; then
    echo "Error: mpi_c_bindings.o not created!" >&2
    exit 1
fi

${FC} -c ${FFLAGS} mpi.f90
if [ ! -f mpi.o ]; then
    echo "Error: mpi.o not created!" >&2
    exit 1
fi

${FC} -c ${FFLAGS} psi_io.f90
if [ ! -f psi_io.o ]; then
    echo "Error: psi_io.o not created!" >&2
    exit 1
fi

${FC} -E ${FFLAGS} --cpp pot3d.F90 > pot3d_cpp.f90
if [ ! -f pot3d_cpp.f90 ]; then
    echo "Error: pot3d_cpp.f90 not created!" >&2
    exit 1
fi

${FC} -c ${FFLAGS} --implicit-interface pot3d_cpp.f90 --show-asr > /dev/null 2>&1
rm -f *.mod *.o
