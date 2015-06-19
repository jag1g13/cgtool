#include <string>
#include <vector>

#include "common.h"

using std::string;
using std::vector;

int main(const int argc, const char *argv[]){
    const string version_string =
            "CGTOOL v0.3.247:e9eaf2ee95e9";

    const string help_header =
            "CGTOOL James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Performs mapping from atomistic to coarse-grained molecular dynamics\n"
            "trajectories and outputs a GROMACS ITP file containing the full mapping,\n"
            "equilibrium bond parameters and force constants.\n\n"
            "Requires GROMACS XTC and ITP files for the atomistic simulation and a\n"
            "configuration file as input.  The config file provides the mapping and\n"
            "bond parameters to be calculated as well as serveral other options.\n\n"
            "Usage:\n"
            "cgtool -c <cfg file> -x <xtc file> -i <itp file>\n";
    // Option syntax is <long flag> \t <comment> \t <flag type> [ \t <default value>]
    // Flag types are 0 - path, 1 - string, 2 - int, 3 - float, 4 - bool
    const string help_options =
            "--cfg\tCGTOOL config file\t0\n"
            "--xtc\tGROMACS XTC file\t0\n"
            "--itp\tGROMACS ITP file\t0\n"
            "--gro\tGROMACS GRO file\t0\n"
            "--fld\tGROMACS forcefield file\t0\n"
            "--dir\tDirectory containing all of the above\t0\n"
            "--frames\tNumber of frames to read\t2\t-1\n"
            "--csv\tOutput bond measurements to CSV\t4\t0";

    Common common;
    common.setHelpStrings(version_string, help_header, help_options);

    vector<string> req_files = {"cfg", "xtc", "gro"};
    vector<string> opt_files = {"itp", "fld"};

    common.collectInput(argc, argv, req_files, opt_files);
    int ok = common.run();
    printf("This segmentation fault doesn't break the program output\n");
    return ok;
}


