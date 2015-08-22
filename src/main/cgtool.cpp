#include "cgtool.h"

#include <string>
#include <vector>

#include <sysexits.h>

#include "itp_writer.h"

using std::string;
using std::vector;

int main(const int argc, const char *argv[]){
    const string version_string =
            "CGTOOL v0.3"
            #include "revision_number.inc"
    ;

    const string help_header =
            "CGTOOL James Graham <J.A.Graham@soton.ac.uk> University of Southampton\n\n"
            "Performs mapping from atomistic to coarse-grained molecular dynamics\n"
            "trajectories and outputs a GROMACS ITP file containing the full mapping,\n"
            "equilibrium bond parameters and force constants.\n\n"
            "Requires GROMACS XTC and ITP files for the atomistic simulation and a\n"
            "configuration file as input.  The config file provides the mapping and\n"
            "bond parameters to be calculated as well as serveral other options.\n\n"
            "Usage:\n"
            "cgtool -c <cfg file> -x <xtc file> -g <gro file>\n";
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

    const string compile_info =
            #include "compile_info.inc"
    ;

    Cgtool cgtool;
    cgtool.setHelpStrings(version_string, help_header, help_options, compile_info);

    vector<string> req_files = {"cfg", "xtc", "gro"};
    vector<string> opt_files = {"itp", "fld"};

    cgtool.collectInput(argc, argv, req_files, opt_files);
    return cgtool.run();
}


void Cgtool::readConfig(){
    Parser cfg_parser(inputFiles_["cfg"].name);

    if(cfg_parser.findSection("membrane")){
        printf("CGTOOL no longer performs membrane analysis - use RAMSi\n");
        exit(EX_USAGE);
    }

    settings_["map"]["on"] =
            cfg_parser.findSection("mapping");

    settings_["bonds"]["on"] =
            cfg_parser.findSection("length") || cfg_parser.findSection("angle") ||
            cfg_parser.findSection("dihedral");

    settings_["csv"]["on"] =
            cfg_parser.findSection("csv");
    settings_["csv"]["molecules"] =
            cfg_parser.getIntKeyFromSection("csv", "molecules", 10000);

    settings_["rdf"]["on"] =
            cfg_parser.findSection("rdf");
    settings_["rdf"]["freq"] =
            cfg_parser.getIntKeyFromSection("rdf", "freq", 1);
    settings_["rdf"]["cutoff"] =
            cfg_parser.getIntKeyFromSection("rdf", "cutoff", 2);
    settings_["rdf"]["resolution"] =
            cfg_parser.getIntKeyFromSection("rdf", "resolution", 100);

    settings_["field"]["on"] =
            cfg_parser.findSection("field");
    settings_["field"]["freq"] =
            cfg_parser.getIntKeyFromSection("field", "freq", 100);
    settings_["field"]["resolution"] =
            cfg_parser.getIntKeyFromSection("field", "resolution", 100);
    settings_["field"]["export"] =
            cfg_parser.getIntKeyFromSection("field", "export", 100);

    string file_format = cfg_parser.getStringKeyFromSection("output", "program", "GROMACS");
    boost::to_upper(file_format);
    outProgram_ = getFileFormat.at(file_format);
    string field_format = cfg_parser.getStringKeyFromSection("output", "field", "MARTINI");
    boost::to_upper(field_format);
    outField_ = getFieldFormat.at(field_format);

    string potential = cfg_parser.getStringKeyFromSection("output", "bond", "HARMONIC");
    potentialTypes_[0] = getPotential.at(potential);
    potential = cfg_parser.getStringKeyFromSection("output", "angle", "COSSQUARED");
    potentialTypes_[1] = getPotential.at(potential);
    potential = cfg_parser.getStringKeyFromSection("output", "dihedral", "HARMONIC");
    potentialTypes_[2] = getPotential.at(potential);
}


void Cgtool::setupObjects(){
    // Open files and do setup
    frame_ = new Frame(inputFiles_["xtc"].name, inputFiles_["gro"].name, &residues_);
    if(inputFiles_["itp"].exists) frame_->initFromITP(inputFiles_["itp"].name);
    if(inputFiles_["fld"].exists) frame_->initFromFLD(inputFiles_["fld"].name);
    for(Residue &res : residues_) res.print();

    if(settings_["map"]["on"]){
        cgMap_ = new CGMap(&residues_, &cgResidues_);
        cgMap_->fromFile(inputFiles_["cfg"].name);
        cgFrame_ = new Frame(*frame_, &cgResidues_);
        cgMap_->initFrame(*frame_, *cgFrame_);
        cgFrame_->setupOutput();
        if(settings_["bonds"]["on"])
            bondSet_ = new BondSet(inputFiles_["cfg"].name, &cgResidues_, potentialTypes_);
    }else{
        // If not mapping make both frames point to the same thing
        cgFrame_ = frame_;
        if(settings_["bonds"]["on"])
            bondSet_ = new BondSet(inputFiles_["cfg"].name, &residues_, potentialTypes_);
    }

    if(settings_["rdf"]["on"])
        rdf_ = new RDF(&residues_, settings_["rdf"]["cutoff"]/100.,
                       settings_["rdf"]["resolution"]);

    if(settings_["field"]["on"]){
        if(settings_["map"]["on"]){
            field_ = new FieldMap(settings_["field"]["resolution"], &residues_, &cgResidues_);
        }else{
            printf("ERROR: Option 'field' requires 'mapping'\n");
            exit(EX_USAGE);
        }
    }
}

void Cgtool::mainLoop(){
    // Calculate bonds and store in BondStructs
    if(settings_["map"]["on"]){
        cgMap_->apply(*frame_, *cgFrame_);
        cgFrame_->writeToXtc();
        if(settings_["bonds"]["on"]) bondSet_->calcBondsInternal(*cgFrame_);

        // Calculate electric field/dipoles
        if(settings_["field"]["on"]){
            if(currFrame_ % settings_["field"]["freq"] == 0){
                field_->calculate(*frame_, *cgFrame_, *cgMap_);
            }
            if(currFrame_ % settings_["field"]["export"] == 0){
                field_->printFieldsToFile();
            }
        }
    }else{
        if(settings_["bonds"]["on"]) bondSet_->calcBondsInternal(*frame_);
    }

    if(settings_["rdf"]["on"] && currFrame_ % settings_["rdf"]["freq"] == 0){
        rdf_->calculateRDF(*frame_);
    }
}

void Cgtool::postProcess(){
    if(settings_["bonds"]["on"]){
        bondSet_->BoltzmannInversion();

        printf("Printing results to ITP\n");
        ITPWriter itp(&residues_, outProgram_, outField_);
        if(cgFrame_->atomHas_.lj) itp.printAtomTypes(*cgMap_);

        if(settings_["map"]["on"]){
            itp.printAtoms(*cgMap_);
        } else{
            bondSet_->calcAvgs();
        }

        itp.printBonds(*bondSet_);

        // Write out all frame bond lengths/angles/dihedrals to file
        // This bit is slow - IO limited
        if(settings_["csv"]["on"])
            bondSet_->writeCSV(settings_["csv"]["molecules"]);

        // Print something so to check results by eye
        for(int j=0; j<6 && j<bondSet_->bonds_.size(); j++){
            printf("%8.4f", bondSet_->bonds_[j].avg_);
        }
        if(bondSet_->bonds_.size() > 6) printf("  ...");
        printf("\n");
    }

    if(settings_["map"]["on"]) cgFrame_->printGRO();
    if(settings_["rdf"]["on"]) rdf_->normalize();
}

Cgtool::~Cgtool(){
    if(bondSet_) delete bondSet_;
    if(rdf_) delete rdf_;
    if(field_) delete field_;
}
