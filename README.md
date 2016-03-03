
**Please consider using the Python implementation PyCGTOOL at [https://github.com/jag1g13/pycgtool](https://github.com/jag1g13/pycgtool) as it has several advantages over this version.**

# CGTOOLS v0.4
##README

The aim of this project is to provide a tool to aid in parametrising coarse-grained (CG) molecular mechanics models.  CGTOOL generates coarse-grained models from atomistic trajectories using a user-provided mapping.  Equilibrium values and force constants of bonded terms are calculated by Boltzmann Inversion of histograms collected from the input trajectory allowing good replication of target properties.
Input is GROMACS XTC and GRO files, along with a custom config file based which specifies the mapping and other parameters.
The output is a GROMACS/MARTINI style forcefield file and a coarse-grained representation of the system created in GROMACS GRO format to allow initial simulations to be started with minimal hand tuning.


Experimental support is present for output in LAMMPS formats to be used with the ELBA forcefield, however this is not yet fully implemented.

These programs are work-in-progress.

Input is GROMACS XTC and GRO files, along with a custom config file based which specifies options for the analysis.  Membrane properties can be averaged over the complete simulation trajectory or output in batches.  An example config file lists the available options.

### How do I get set up? ###

Required to compile:

* CMake 2.8.4 or newer
* C++ compiler supporting the C++11 standard
* Boost C++ libraries (program\_options module recommended but optional)
* Optional: Doxygen to build developer documentation

To compile the program:

* Create a new directory 'build' within the main distribution directory
* From the build directory `cmake ..`
* `make` to build the executables
* `make check` to compile and run tests
* `make doc` for developer documentation - requires Doxygen

To use the programs:

CGTOOL
* Help text is available with `cgtool -h` or `cgtool --help`
* Required inputs are a GROMACS GRO and XTC file and a config file
* The program should be called using `cgtool -c <cfg file> -x <xtc file> -g <gro file>` (order not important)
* An optional GROMACS ITP file may be provided with the `-i <itp file>` option to allow calculation of charges
* The config file specifies the mapping to be applied, an example is present in the test\_data directory

RAMSi
* Help text is available with `ramsi -h` or `ramsi --help`
* The program should be called using `ramsi  -c <CFG file> -x <XTC file> -g <GRO file>` (order not important)
* A config file is required which specifies the analysis options, in the format seen in the examples directory

### Testing ###
The Bitbucket repo is polled every 15 minutes by a Jenkins build server for unit and integration testing.  Builds are tested automatically on Ubuntu and CentOS7 with occasional manual tests on OSX and RHEL6.
Currently only a few source files have complete unit tests.  Integration tests ensure that the program compiles successfully and produces correct output for a test dataset.

### Contact ###

* James Graham: <J.A.Graham@soton.ac.uk>