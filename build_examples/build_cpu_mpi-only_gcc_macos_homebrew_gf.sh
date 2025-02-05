# build script which uses GFortran as it's compiler
# along with using our own MPI wrappers and call C MPI instead
# and doesn't use any HDF5
POT3D_HOME=$PWD
MPICC=${CONDA_PREFIX}/bin/mpicc
FC=${CONDA_PREFIX}/bin/gfortran
FFLAGS="-O3 -march=native -lmpi"

echo "==> Using Fortran compiler: " ${FC}
echo "=== Starting build process ==="

echo "Entering src directory"
cd src
${MPICC} -c mpi_wrapper.c
${FC} -c ${FFLAGS} mpi_c_bindings.f90
${FC} -c ${FFLAGS} mpi.f90
${FC} -c ${FFLAGS} psi_io.f90
${FC} -E ${FFLAGS} -cpp pot3d.F90 > pot3d_cpp.f90
${FC} -c ${FFLAGS} pot3d_cpp.f90 
OBJS="mpi_wrapper.o mpi_c_bindings.o mpi.o psi_io.o pot3d_cpp.o"

${FC} ${FFLAGS} ${OBJS} -o pot3d

if [ ! -e pot3d ]; then
  echo "!!> ERROR!  pot3d executable not found.  Build most likely failed."
  echo "     Contents of src/build.err:"
  cat build.err
  echo ""
  exit 1
fi

echo "current working directory is: ${PWD}"

echo "POT3D executable created successfully"

echo "==> Copying pot3d executable to: ${POT3D_HOME}/bin/pot3d"
cp pot3d ${POT3D_HOME}/bin/
echo "==> Build complete!"
echo "Please add the following to your shell startup (e.g. .bashrc, .profile, etc.):"
echo "export PATH=${POT3D_HOME}/bin:\$PATH"
