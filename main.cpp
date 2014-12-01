#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <ctime>
#include <map>

//#include <unistd.h>
#include <rpc/rpc.h>
#include <sys/stat.h>

#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO

#include <gromacs/fileio/xtcio.h>

#endif

#include "frame.h"
#include "cg_map.h"

#define DEBUG true

/*things from std that get used a lot*/
using std::ifstream;
using std::string;
using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::clock_t;

/*prototype functions*/
//TODO convert int functions to bool


vector<float> calc_bond_lens(Frame *, vector<float>);

float calc_avg(vector<float>);

void setup_cg_map(Frame *, Frame *);

void cg_map(Frame *, Frame *);

void split_text_output(const char *, clock_t);

bool file_exists(const char *);

//extern int read_next_xtc(t_fileio *fio,
//        int natoms, int *step, real *time,
//        matrix box, rvec *x, real *prec, gmx_bool *bOK);


int main(int argc, char *argv[]){
    int ok;
    float time, prec;
    rvec *x;
    matrix box;
    t_fileio *xtc;
    vector<float> bond_lens;
    char mode[1] = {'r'};

    /* Where does the user want us to look for GRO and XTC files? */
    clock_t start = std::clock();
    split_text_output("Identifying files", start);
    char groname[40], xtcname[40], mapname[40];
    int natoms, step;
    if(argc < 2){
        cout << "Using current directory" << endl;
        strcpy(groname, "npt.gro");
        strcpy(xtcname, "md.xtc");
        strcpy(mapname, "sacc.map");
    } else if(argc == 2){
        cout << "Using directory provided" << endl;
        strcpy(groname, argv[1]);
        strcat(groname, "/npt.gro");
        strcpy(xtcname, argv[1]);
        strcat(xtcname, "/md.xtc");
        strcpy(mapname, argv[1]);
        strcat(mapname, "/sacc.map");
    } else if(argc == 4){
        cout << "Using filenames provided" << endl;
        strcpy(groname, argv[1]);
        strcpy(xtcname, argv[2]);
        strcpy(mapname, argv[3]);
    } else{
        cout << "Wrong number of arguments given" << endl;
        throw std::runtime_error("Wrong number of arguments");
    }
    if(!file_exists(groname) || !file_exists(xtcname) || !file_exists(mapname)){
        cout << "Input file does not exist" << endl;
        throw std::runtime_error("File doesn't exist");
    }
    cout << "GRO file: " << groname << endl;
    cout << "XTC file: " << xtcname << endl;
    cout << "MAP file: " << mapname << endl;

    /* Open files and do setup */
    split_text_output("Frame setup", start);
    Frame frame = Frame(0, 0, "");
    xtc = open_xtc(xtcname, mode);
    ok = frame.setupFrame(groname, xtc, &natoms, &step, &time, box, &x, &prec);
    Frame cg_frame = Frame(&frame);
    CGMap mapping(mapname);
    mapping.initFrame(&frame, &cg_frame);

    /* Keep reading frames until something goes wrong (run out of frames) */
    split_text_output("Reading frames", start);
    start = std::clock();
    int i = 0;
    while(frame.readNext(xtc, &natoms, &step, &time, box, &x, &prec)){
        /* Process each frame as we read it, frames are not retained */
        //cg_map(&frame, &cg_frame);
        bond_lens = calc_bond_lens(&frame, bond_lens);
        if(i % 1000 == 0){
            cout << "Read " << i << " frames\r";
            std::flush(cout);
        }
        //usleep(1000);
        i++;
    }
    cout << "Read " << i << " frames" << endl;

    /* close remaining files */
    close_xtc(xtc);

    /* Post processing */
    split_text_output("Post processing", start);
    cout << "Avg bond[0] " << calc_avg(bond_lens) << endl;
    return !ok;
}



/*def map_cg_solvent_within_loop(curr_frame, frame, cg_frame=0):
    """
    perform CG mapping using cg_map list of lists
            with current cg_map does a simple heavy atom mapping

    will be CM or GC depending on how 'Atom.mass' was set previously
    if mapping is changed in cg_map to include other atoms

            should remove the setup code into its own function (or the main xtc setup)
    """
    global cg_atom_nums
    if curr_frame == 0:
        cg_frame = Frame(curr_frame, cg_atom_nums)
    cg_frame.num = curr_frame
    for i, site in enumerate(cg_sites):
        coords = np.zeros(3)
        tot_mass = 0.
        charge = 0.
        for atom in cg_map[i]:
            mass = frame.atoms[sugar_atom_nums[atom]].mass
            tot_mass = tot_mass + mass
            coords = coords + mass*frame.atoms[sugar_atom_nums[atom]].loc
            charge = charge + frame.atoms[sugar_atom_nums[atom]].charge
        coords /= tot_mass  # number of atoms cancels out
        if curr_frame == 0:
            cg_frame.atoms.append(Atom(site, coords, charge))
        else:
            cg_frame.atoms[i] = Atom(site, coords, charge)
        if curr_frame == 0:
            cg_atom_nums[site] = i
    j = len(cg_sites)
    for atom in frame.atoms:
        if atom.atom_type == "OW":
            if curr_frame == 0:
                cg_frame.atoms.append(Atom("OW", atom.loc, 0.0))
            else:
                cg_frame.atoms[j] = Atom("OW", atom.loc, 0.0)
            j += 1
    return cg_frame*/

void setup_cg_map(Frame *frame, Frame *cg_frame){
    /**
    * \brief Setup a Frame to hold the coarse grained trajectory
    */
    throw std::logic_error("Not implemented");
    int num_cg_atoms = 0;
    for(int i = 0; i < frame->num_atoms_; i++){

    }
    cg_frame->allocateAtoms(num_cg_atoms);
}

void cg_map(Frame *frame, Frame *cg_frame){
    throw std::logic_error("Not implemented");
}

vector<float> calc_bond_lens(Frame *frame, vector<float> bond_lens){
    /**
    * \brief Calculate all relevant bond lengths in Frame
    */
    bond_lens.push_back(frame->bondLength(0, 1));
    //cout << bond_lens.back() << endl;
    return bond_lens;
}

void read_cg_mapping(vector<char *> bead_names,
        std::map<string, vector<string>> bead_map){
    throw std::logic_error("Not implemented");
}

bool get_bond_measures(string dir){
    /**
    * \brief Read in all required bond length, angle and dihedral measurements from file
    */
    throw std::logic_error("Not implemented");
    string filename = dir + "bonds.conf";
    return false;
}

float calc_avg(vector<float> vec){
    /**
    * \brief Calculate the average of a Vector<float>
    */
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    float mean = sum / vec.size();
    return mean;
}

void split_text_output(const char *name, clock_t start){
    clock_t now = std::clock();
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << (float) (now - start) / CLOCKS_PER_SEC << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

bool file_exists(const char *name){
    struct stat buffer;
    return (stat(name, &buffer) == 0);
}