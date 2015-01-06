#ifndef FRAME_H_
#define FRAME_H_

#ifndef INCLUDE_GMXFIO
#define INCLUDE_GMXFIO
#include <gromacs/fileio/xtcio.h>
#endif

#include <vector>
#include <string>
#include <map>

#include <string.h>

#include "bond_struct.h"

/**
* \brief Struct to hold atom data
*/
struct Atom{
    /** A serial number; no longer needed */
    int atom_num;
    /** Residue number and name in the GRO file */
    //char resid[10];
    std::string resid;
    /** A three character atom type specifier; e.g. "OH1"; Should replace with a string */
    char atom_type[4];
    /** Atomtype as a string.  I don't want to be dealing with *char */
    std::string atom_type_string;
    /** Atomic coordinates in x, y, z */
    float coords[3];
    /** Atomic charge from the force field */
    float charge;
    /** Atomic mass */
    float mass;
    /** Vector of pointers to Atom listing all atoms bonded to this one */
    std::vector<Atom *> neighbours;
    /** Create an Atom and set its number */
    Atom(int num){atom_num = num;};
    /** Create a blank Atom instance */
    Atom(){};
};

/**
* \brief Struct to hold CG bead data; inherits from Atom
*/
struct CGBead : Atom{
    /** Vector of atoms that this CG bead represents */
    std::vector<std::string> sub_atoms;
};

/**
* \brief A residue from the GRO file, contains pointers to atoms
*/
struct Residue{
    /** The name of this Residue */
    std::string res_name;
    //char res_name[10];
    /** Atoms contained within this residue */
    std::vector<int> atoms;
    /** Atoms contained within this residue */
    //std::vector<std::string> atom_names;
    std::vector<char *> atom_names;
    /** Constructor to set res_name */
    //Residue(const char* tmp){strcpy(res_name, tmp);};
    Residue(const std::string tmp){res_name = tmp;};
    /** Blank constructor */
    Residue();
};


/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
public:
    /** Has the Frame been properly setup yet */
    bool isSetup_ = false;
    /** The number of this Frame, starting at 0 */
    int num_;
    /** The simulation step corresponding to this frame */
    int step_;
    //TODO refactor these names to match the new style
    /** The number of atoms stored in this frame */
    int num_atoms_;
    /** The number of atoms stored in this frame that we find interesting */
    int numAtomsTrack_;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** Vector of Residues; Each Residue contains pointers to atoms */
    std::vector<Residue> residues_;
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
    std::map<int, std::string> num_to_name_;
    /** Dictionary mapping atom names to numbers */
    std::map<std::string, int> name_to_num_;

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

    /** Print info for all atoms up to n.  Default print all. */
    void printAtoms(const int n=-1);

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
