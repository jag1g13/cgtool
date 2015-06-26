#include <string>
#include <vector>

#include "common.h"

using std::string;
using std::vector;

int main(const int argc, const char *argv[]){
    const string version_string =
            "RAMSi v0.3.254:0446c2a51604";

    const string help_header =
            "James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Performs several analysis functions for biomembrane simulations in GROMACS.\n\n"
            "Requires GROMACS XTC and GRO files from the simulation and a configuration\n"
            "file as input.  The config file contains the output settings and a\n"
            "coarse-grain mapping if a AA->CG transformation is required.\n\n"
            "Usage:\n"
            "ramsi -c <CFG file> -x <XTC file> -g <GRO file>\n";
    // Option syntax is <long flag> \t <comment> \t <flag type> [ \t <default value>]
    // Flag types are 0 - path, 1 - string, 2 - int, 3 - float, 4 - bool
    const string help_options =
            "--cfg\tRAMSi config file\t0\n"
            "--xtc\tGROMACS XTC file\t0\n"
            "--gro\tGROMACS GRO file\t0";

    Common common;
    common.setHelpStrings(version_string, help_header, help_options);

    vector<string> req_files = {"cfg", "xtc", "gro"};

    common.collectInput(argc, argv, req_files);
    return common.run();
}


