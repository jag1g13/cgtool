#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO
#include <gromacs/fileio/xtcio.h>
#endif

#ifndef FRAME_H
#define FRAME_H

#include <vector>

using std::vector;

/**
* \brief Struct to hold atom data
*/
struct Atom{
    /** A serial number; no longer needed */
    int atom_num;
    /** A three character atom type specifier; e.g. "OH1"; Should replace with a string */
    char atom_type[3];
    /** Atomic coordinates in x, y, z */
    float coords[3];
    /** Atomic charge from the force field */
    float charge;
    /** Atomic mass */
    float mass;
    /** Vector of pointers to Atom listing all atoms bonded to this one */
    vector<Atom*> neighbours;
};


/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
//TODO perhaps move xtc read functions into Frame class
public:
    /** The simulation step corresponding to this frame */
    int step;
    /** The number of atoms stored in this frame */
    int num_atoms;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms;
    /** The simulation time of this frame, in picoseconds */
    float time;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec;
    /** Size of the simulation box */
    matrix box;
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name;

    /**
    * \brief Create Frame passing frame number, number of atoms to store and the frame name
    *
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocate_atoms()
    */
    Frame(int, int, std::string);

    //TODO convert int functions to bool
    /**
    * \brief Allocate space for a number of atoms
    *
    * Used if the number of atoms isn't known at time of creation
    */
    int allocate_atoms(int);

    /**
    * \brief Write Frame to XTC output file
    */
    bool write_to_xtc(t_fileio*);

    /**
    * \brief Calculate distance between two atoms
    */
    float bond_length(int, int);

    /**
    * \brief Calculate angle between vectors a->b and c->d
    *
    * To be used for bond angles (b=c) and dihedrals (b=/=c)
    */
    float bond_angle(int, int, int, int);
};

#endif
