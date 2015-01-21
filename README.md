# README #

The aim of this project is to create a tool to aid in parametrising coarse-grained molecular mechanics models.  The program currently extracts bond parameters from atomistic simulations in GROMACS and performs some initial dipole calculations.  The final aim is to be able to export force field data files directly.

This program is very much work-in-progress and is not yet extensively tested.

### How do I get set up? ###

Required to compile:

* Boost C++ libraries
* Gromacs 5.0.x libraries
* CMake 2.8.4 or newer
* GCC or Clang compiler supporting the C++11 standard
* Google Test to compile tests

The CMake file has been tested on Ubuntu 14.04.1 (GCC/Clang) and Mac OSX 10.6 (GCC) and should work on similar systems.
I intend to make executables for common OSes available once the project is more complete.

To compile the program:

* Create a new directory 'build' within the main distribution directory
* From the build directory 'cmake ..'
* 'make'

To use the program:

* Usage text is available with 'cgtool -h' or 'cgtool --help'
* The program should be called using 'cgtool <path to xtc> <path to gro> <path to cfg> <path to top>'
* A config file is required which specifies the mapping to be applied, in the format seen in the test_data directory

### Contact ###

* James Graham: <J.A.Graham@soton.ac.uk>