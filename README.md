# CGTOOLS v0.4
##README

The aim of this project is to create a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.  The program currently calculates bonded parameters for CG models from atomistic simulations in GROMACS, i.e. bond lengths, angles, dihedrals and their force constants.
Bonded parameters are output into a GROMACS/MARTINI style forcefield file and a coarse-grained representation of the system created in GROMACS GRO or LAMMPS DATA format to allow initial simulations to be started with minimal hand tuning.

The completed tools will be compatible with a range of MD simulators and forcefield types, including support for dipoles (work in progress).  Current targets are GROMACS/MARTINI and LAMMPS/ELBA.

RAMSi stands for Rapid Analysis of Membrane Simulations and is a tool to aid in analysing simulations of biomembranes.  It is able to calculate membrane thickness, curvature and surface area per lipid.  RAMSi makes use of multiple CPU cores via OpenMP.
Input is GROMACS XTC and GRO files, along with a custom config file based which specifies options for the analysis. Membrane properties can be averaged over the complete simulation trajectory or output in batches.  An example config file lists the available options.

These programs are work-in-progress and are not yet extensively tested.

### How do I get set up? ###

Required to compile:

* CMake 2.8.4 or newer
* GCC or Clang compiler supporting the C++11 standard
* Boost C++ libraries with program\_options module
* Optional: Doxygen to build developer documentation

The CMake file has been tested only on Ubuntu 14.04 (GCC/Clang) but should work on similar systems.
I intend to make executables for common OSes available once the project is more complete.

To compile the program:

* Create a new directory 'build' within the main distribution directory
* From the build directory `cmake ..`
* `make` to build the executables
* `make check` to compile and run tests
* `make doc` for developer documentation - requires Doxygen

To use the programs:

CGTOOL
* Help text is available with `cgtool -h`
* Required inputs are a GROMACS GRO and XTC file and a config file
* The program should be called using `cgtool -c <cfg file> -x <xtc file> -g <gro file>`
* The config file specifies the mapping to be applied, an example is present in the test\_data directory

RAMSi
* Help text is available with `ramsi --help`
* The program should be called using `ramsi  -c <CFG file> -x <XTC file> -g <GRO file>`
* A config file is required which specifies the analysis options, in the format seen in the examples directory

### Testing ###
The Bitbucket repo is polled every 15 minutes by a Jenkins build server for unit and integration testing.  Builds are tested on Ubuntu Linux.
Currently only a few source files have complete unit tests, but this is being fixed.  Integration tests ensure that the program compiles successfully and produces correct output for a test dataset.

### Contact ###

* James Graham: <J.A.Graham@soton.ac.uk>
