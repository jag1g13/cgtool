# README #

The aim of this project is to create a tool to aid in parametrising coarse-grained molecular mechanics models.  The program currently extracts bond parameters from atomistic simulations in GROMACS and performs some initial dipole calculations.  The final aim is to be able to export force field data files directly.

This program is very much work-in-progress and is not yet extensively tested.

### How do I get set up? ###

Required to compile:

* Boost C++ libraries including program_options
* FLENS linear algebra headers (will be included in repo in the future)
* CMake 2.8.4 or newer
* GCC or Clang compiler supporting the C++11 standard
* Optional: Google Test to compile tests
* Optional: Doxygen to build developer documentation

The CMake file has been tested on Ubuntu 14.04.1 (GCC/Clang) and Mac OSX 10.6 (GCC) and should work on similar systems.
I intend to make executables for common OSes available once the project is more complete.

To compile the program:

* Create a new directory 'build' within the main distribution directory
* From the build directory 'cmake ..'
* 'make cgtool' to build the main executable
* 'make check' to compile and run tests
* 'make doc' for developer documentation - requires Doxygen

To use the program:

* Help text is available with 'cgtool --help'
* The program should be called using 'cgtool --dir <path to files>' or 'cgtool --cfg <cfg file> --xtc <xtc file> --itp <itp file>'
* A config file is required which specifies the mapping to be applied, in the format seen in the test_data directory

### Testing ###
The Bitbucket repo is polled every 15 minutes by a Jenkins build server for unit and integration testing.  Builds are tested on Ubuntu Linux and Mac OSX.
Currently only a few source files have complete unit tests, but this is being fixed.  Integration tests ensure that the program compiles successfully and produces correct output for several test datasets.

### Contact ###

* James Graham: <J.A.Graham@soton.ac.uk>