# README #

RAMSi stands for Rapid Analysis of Membrane Simulations and is a tool to aid in analysing simulations of biomembranes.  It is able to calculate membrane thickness, curvature and surface area per lipid.  RAMSi makes use of multiple CPU cores via OpenMP.

Input is GROMACS XTC and GRO files, along with a custom config file based which specifies options for the analysis. Membrane properties can be averaged over the complete simulation trajectory or output in batches.

This program is work-in-progress and is not yet extensively tested.

### How do I get set up? ###

Required to compile:

* CMake 2.8.4 or newer
* GCC or Clang compiler supporting the C++11 standard
* Boost C++ libraries and program\_options
* Optional: Doxygen to build developer documentation

To compile the program:

* Create a new directory `build` within the main distribution directory
* From the build directory `cmake ..`
* `make` to build the main executable
* `make check` to compile and run tests
* `make doc` for developer documentation - requires Doxygen

To use the program:

* Help text is available with `ramsi --help`
* The program should be called using `ramsi  -c <CFG file> -x <XTC file> -g <GRO file>`
* A config file is required which specifies the analysis options, in the format seen in the examples directory

### Testing ###
The Bitbucket repo is polled every 15 minutes by a Jenkins build server for unit and integration testing.  Builds are tested on Ubuntu Linux.
Currently only a few source files have complete unit tests, but this is being fixed.  Integration tests ensure that the program compiles successfully and produces correct output for several test datasets.

### Contact ###

* James Graham: <J.A.Graham@soton.ac.uk>
