## Steps to build POT3D on Linux/macOS with custom MPI wrappers and HDF5 dependency

Custom MPI module can be seen at dedicated repo [c_mpi](https://github.com/gxyd/c_mpi/)
Current Branch on POT3D without HDF5 dependencies and our latest MPI module can be seen at [Latest branch](https://github.com/gxyd/POT3D/tree/build_pot3d_custom_mpi_wrappers_and_without_hdf5)

### Clone POT3D with current branch

```console
> git clone git@github.com:gxyd/POT3D.git
> git checkout -b your_branch_name origin/build_pot3d_custom_mpi_wrappers_and_without_hdf5
```

## Create an environment

```console
> cd src
> conda env create -f ci/environment.yml   # this will create pot3d_env environment, which has GFortran and OpenMPI
> conda activate pot3d_env
```

### For compiling POT3D with GFortran

- For Ubuntu and macOS users both can use the same script
```console
> ./build_examples/build_cpu_mpi-only_gcc_gfortran.sh   # To run with GFortran
```

### For compiling POT3D with LFortran

**NOTE**: you might need to modify `FC` value in `./build_examples/build_cpu_mpi-only_gcc_lfortran.sh` script
```console
./build_examples/build_cpu_mpi-only_gcc_lfortran.sh     # To run with LFortran upto --show-asr level
```

### Run with a sample input

```console
./validate.sh
```

*The output should look like this*

```console
> ./build_examples/build_cpu_mpi-only_gcc_gfortran.sh
=== STARTING POT3D BUILD ===
==> Entering src directory...
==> Removing old Makefile...
==> Generating Makefile from Makefile.template...
==> Compiling code...
==> Copying pot3d executable to: /home/aditya-trivedi/thats_me/main/POT3D/bin/pot3d
==> Build complete!
    Please add the following to your shell startup (e.g. .bashrc, .profile, etc.):
    export PATH=/home/aditya-trivedi/thats_me/main/POT3D/bin:$PATH

> ./validate.sh 
Running POT3D with 1 MPI rank...
Done!
Wall clock time:                16.716739 seconds
 
Run has PASSED validation!
 
Running POT3D with 2 MPI ranks...
Done!
Wall clock time:                13.319068 seconds
 
Run has PASSED validation!

```

## (NOT important section) Steps to build POT3D on Ubuntu/macOS with HDF5 as dependency

**NOTE**: this section might need an update anyways

make sure we've `hdf5` installed as a dependency. You can do so via homebrew as:
```bash
$ brew install hdf5
```

it uses `FC=mpif90` in `build_examples/build_cpu_mpi-only_gcc_macos_homebrew_lf.sh`

```bash
$ git clone https://github.com/gxyd/POT3D
$ cd POT3D
$ git checkout -t origin/build_pot3d_unmodified
$ git checkout b547560acbb5e92de9f4220ce3a5eed5a79ac1b7
$ ./build_examples/build_cpu_mpi-only_gcc_macos_homebrew_lf.sh
$ export PATH=$PWD/bin:$PATH
```

## Running with a sample input

and then run:
```bash
$ ./validate.sh
```

the above command will produce an output like:
```console
Running POT3D with 1 MPI rank...
Done!
Wall clock time:                15.784941 seconds
 
Run has PASSED validation!
 
Running POT3D with 2 MPI ranks...
Done!
Wall clock time:                12.093250 seconds
 
Run has PASSED validation!
```
