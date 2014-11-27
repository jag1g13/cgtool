#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO
#include <gromacs/fileio/xtcio.h>
#endif

#ifndef ATOMSTRUCT
#define ATOMSTRUCT
/**
* \brief Struct to hold atom data
*/
struct Atom{
    int atom_num;
    char atom_type[3];
    float coords[3];
    float charge;
    float mass;
    //Atom **neighbours; // list of pointers to neighbours
};
#endif



/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> an contains member functions to operate on this
*/
class Frame{
//TODO perhaps move xtc read functions into Frame class
public:
    int step, num_atoms;        // step is frame_num, wanted to be consistent with GMX
    std::vector<Atom> atoms;    // list of atoms
    //Atom* atoms;
    float time, prec;
    matrix box;
    std::string name;
    
    Frame(int, int, std::string);

    //TODO convert int functions to bool
    int allocate_atoms(int);

    bool write_to_xtc(t_fileio*);
    
    float bond_length(int, int);

    float bond_angle(int, int, int, int);
};
