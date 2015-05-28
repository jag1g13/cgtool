
#include <iostream>
#include <vector>

#include <sysexits.h>
#include <locale.h>

#include "frame.h"
#include "cg_map.h"
#include "itp_writer.h"
#include "parser.h"
#include "small_functions.h"
#include "file_io.h"
#include "field_map.h"
#include "cmd.h"
#include "rdf.h"

using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;


int main(const int argc, const char *argv[]){
    const double very_start = start_timer();
    double start = very_start;

    const string version_string =
            "CGTOOL v0.3.234:d59a67674392";

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
            "--csv\tOutput bond measurements to CSV\t4\t0\n"
            "--nomap\tDon't perform cg mapping\t4\t0\n"
            "--fcround\tRound force constants\t4\t0\n"
            "--field\tCalculate electric field\t4\t0";

    // Allow comma separators in numbers for printf
    setlocale(LC_ALL, "");

    // ##############################################################################
    // Input Collection
    // ##############################################################################

    // Get input files
    split_text_output(version_string, start);
    CMD cmd_parser(help_header, help_options, argc, argv);
    const string cfgname = cmd_parser.getFileArg("cfg");
    const string xtcname = cmd_parser.getFileArg("xtc");
    const string itpname = cmd_parser.getFileArg("itp");
    const string groname = cmd_parser.getFileArg("gro");
    const string fldname = cmd_parser.getFileArg("fld");

    cout << "CFG file: " << cfgname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "GRO file: " << groname << endl;
    if(itpname != "") cout << "ITP file: " << itpname << endl;
    if(fldname != "") cout << "FLD file: " << fldname << endl;
    // GRO file is not required but helps find residues
    if(!file_exists(xtcname) || !file_exists(cfgname) || !file_exists(groname)){
        cout << "Input file does not exist" << endl;
        exit(EX_NOINPUT);
    }

    // ##############################################################################
    // System Setup
    // ##############################################################################

    // Read number of frames from config, if number not found read them all
    Parser cfg_parser(cfgname);
    vector<string> tokens;
    int num_frames_max = -1;
    if(cfg_parser.getLineFromSection("frames", tokens, 1)) num_frames_max = stoi(tokens[0]);
    if(cmd_parser.getIntArg("frames") != 0) num_frames_max = cmd_parser.getIntArg("frames");

    bool do_map = !cmd_parser.getBoolArg("nomap") && cfg_parser.findSection("mapping");

    vector<Residue> residues;
    while(cfg_parser.getLineFromSection("residues", tokens, 1)){
        residues.emplace_back(Residue());
        Residue *res = &residues.back();
        res->resname = tokens[0];

        if(tokens.size() > 1) res->num_residues = stoi(tokens[1]);
        if(tokens.size() > 2){
            res->num_atoms = stoi(tokens[2]);
            res->calc_total();
            res->populated = true;
        }
        if(tokens.size() > 3) res->ref_atom_name = tokens[3];
        if(tokens.size() > 4) throw std::runtime_error("Old input file");
    }

    const int num_residues = residues.size();
    bool pop_so_far = residues[0].populated;
    residues[0].start = 0;
    for(int i=1; i<num_residues; i++){
        if(pop_so_far){
            residues[i].start = residues[i - 1].total_atoms + residues[i - 1].start;
            pop_so_far = residues[i].populated;
        }
    }

    // Open files and do setup
    split_text_output("Frame setup", start);
    Frame frame(xtcname, groname, residues);
    if(itpname != "") frame.initFromITP(itpname);
    if(fldname != "") frame.initFromFLD(fldname);
    for(Residue &res : residues) res.print();

    Frame cg_frame(frame);
    BondSet bond_set(cfgname, residues);
    CGMap mapping(residues);
    if(do_map){
        mapping.fromFile(cfgname);
        mapping.initFrame(frame, cg_frame);
        mapping.correctLJ();
        cg_frame.setupOutput();
    }

    bool do_rdf = cfg_parser.findSection("rdf");
    const int rdf_freq = cfg_parser.getIntKeyFromSection("rdf", "freq", 1);
    RDF rdf;
    if(do_rdf){
        const double cutoff = cfg_parser.getDoubleKeyFromSection("rdf", "cutoff", 2.);
        const int resolution = cfg_parser.getIntKeyFromSection("rdf", "resolution", 100.);
        rdf.init(residues, cutoff, resolution);
    }

    bool do_field = cfg_parser.findSection("field");
    const int electric_field_freq = cfg_parser.getIntKeyFromSection("field", "freq", 100);
    FieldMap field;
    if(do_field) field.init(100, 100, 100, mapping.numBeads_);

    // Read and process simulation frames
    split_text_output("Reading frames", start);
    start = start_timer();
    const int full_xtc_frames = get_xtc_num_frames(xtcname);
    printf("Total of %'6d frames in XTC\n", full_xtc_frames);
    if(num_frames_max < 0){
        printf("Reading all frames from XTC\n");
    }else{
        printf("Reading %'6d frames from XTC\n", num_frames_max);
    }

    // ##############################################################################
    // Main loop
    // ##############################################################################

    int i = 1;
    int progress_update_loc = 0;
    const int progress_update_freqs[6] = {1, 2, 5, 10, 100, 1000};
    double last_update = start_timer();
    // Keep reading frames until something goes wrong (run out of frames) or hit limit
    while(frame.readNext() && (num_frames_max < 0 || i < num_frames_max)){
        // Process each frame as we read it, frames are not retained
        #ifdef UPDATE_PROGRESS
        if(i % progress_update_freqs[progress_update_loc] == 0){
            // Set time between progress updates to nice number
            const double time_since_update = end_timer(last_update);
            if(time_since_update > 0.5f && progress_update_loc > 0){
                progress_update_loc--;
            }else if(time_since_update < 0.01f && progress_update_loc < 6){
                progress_update_loc++;
            }

            const double time = end_timer(start);
            const double fps = i / time;

            double t_remain = (num_frames_max - i) / fps;
            if(num_frames_max < 0) t_remain = (full_xtc_frames - i) / fps;
            printf("Read %'9d frames @ %'d FPS %6.1fs remaining\r", i, static_cast<int>(fps), t_remain);
            std::flush(cout);

            last_update = start_timer();
        }
        #endif

        // Calculate bonds and store in BondStructs
        if(do_map){
            mapping.apply(frame, cg_frame);
            cg_frame.writeToXtc();
            bond_set.calcBondsInternal(cg_frame);
        }else{
            bond_set.calcBondsInternal(frame);
        }

        // Calculate electric field/dipoles
        if(do_field && i % electric_field_freq == 0){
            field.calculate(frame, cg_frame, mapping);
        }

        if(do_rdf && i % rdf_freq == 0) rdf.calculateRDF(frame);

        i++;
    }

    // ##############################################################################
    // Post processing / Averaging
    // ##############################################################################

    // Print some data at the end
    cout << string(80, ' ') << "\r";
    printf("Read %'9d frames", i);
    const double time = end_timer(start);
    const double fps = i / time;
    printf(" @ %'d FPS", static_cast<int>(fps));
    if(num_frames_max == -1){
        // Bitrate (in MiBps) of XTC input - only meaningful if we read whole file
        const double bitrate = file_size(xtcname) / (time * 1024 * 1024);
        printf("%6.1f MBps", bitrate);
    }
    printf("\n");

    // Post processing
    split_text_output("Post processing", start);
    if(do_map){
        cg_frame.printGRO();
        bond_set.BoltzmannInversion();

        const FileFormat file_format = FileFormat::GROMACS;
        const FieldFormat field_format = FieldFormat::MARTINI;

        cout << "Printing results to ITP" << endl;
        //TODO put format choice in config file or command line option
        ITPWriter itp(residues, file_format, field_format);
        if(cg_frame.atomHas_.lj) itp.printAtomTypes(mapping);
        itp.printAtoms(mapping);
        itp.printBonds(bond_set, cmd_parser.getBoolArg("fcround"));
    }else{
        bond_set.calcAvgs();
    }

    // Write out all frame bond lengths/angles/dihedrals to file
    // This bit is slow - IO limited
    if(cmd_parser.getBoolArg("csv")) bond_set.writeCSV();

    // Calculate RDF
    if(do_rdf) rdf.normalize();

    // Print something so I can check results by eye
    for(int j=0; j<6 && j<bond_set.bonds_.size(); j++){
        printf("%8.4f", bond_set.bonds_[j].avg_);
    }
    if(bond_set.bonds_.size() > 6) printf("  ...");
    cout << endl;

    // Final timer
    split_text_output("Finished", very_start);
    return EX_OK;
}


