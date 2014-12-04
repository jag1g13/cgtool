#ifndef FRAME_H_
#define FRAME_H_

#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO
#include <gromacs/fileio/xtcio.h>
#endif

#include <vector>
#include <string>
#include <map>

#include "bond_struct.h"

using std::vector;
using std::string;

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
    vector<Atom *> neighbours;
};

/**
* \brief Struct to hold CG bead data; inherits from Atom
*/
struct CGBead : Atom{
    /** Vector of atoms that this CG bead represents */
    vector<string> sub_atoms;
};


/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
public:
    /** The number of this Frame, starting at 0 */
    int num_;
    /** The simulation step corresponding to this frame */
    int step_;
    /** The number of atoms stored in this frame */
    int num_atoms_;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** The simulation time of this frame, in picoseconds */
    float time_;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec_;
    /** Size of the simulation box */
    matrix box_;
    /** Holds atomic coordinates for GROMACS */
    rvec *x_;
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_;
    /** Dictionary mapping atom numbers to atom names */
    std::map<int, string> num_to_name_;
    /** Dictionary mapping atom names to numbers */
    std::map<string, int> name_to_num_;

    /**
    * \brief Create Frame passing frame number, number of atoms to store and the frame name
    *
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocateAtoms()
    */
    Frame(int, int, std::string);

    /**
    * \brief Create Frame by copying name and step from another Frame
    *
    * Intended for creating a CG Frame from an atomistic one
    */
    Frame(const Frame*);

    /**
    * \brief Sets up Frame from XTC and GRO files
    *
    * Reads in first frame of XTC and allocates atoms
    */
    bool setupFrame(const char *groname, const char *topname, t_fileio *xtc);

    /**
    * \brief Read next frame from the open XTC file
    */
    bool readNext(t_fileio *xtc);

    /**
    * \brief Allocate space for a number of atoms
    *
    * Used if the number of atoms isn't known at time of creation
    */
    int allocateAtoms(int);

    /**
    * \brief Write Frame to XTC output file
    */
    bool writeToXtc(t_fileio *);

    /**
    * \brief Calculate distance between two atoms
    */
    float bondLength(int, int);

    /**
    * \brief Calculate distance between two atoms in a BondStruct object
    *
    * Wrapper around float bondLength(int, int)
    */
    float bondLength(BondStruct *bond);

    /**
    * \brief Calculate angle between vectors a->b and c->d
    *
    * To be used for bond angles (b=c) and dihedrals (b=/=c)
    */
    float bondAngle(int, int, int, int);


    /**
    * \brief Calculate angle or dihedral between atoms in a BondStruct object
    *
    * Wrapper around float bondAngle(int, int, int, int)
    */
    float bondAngle(BondStruct *bond);
};

#endif
