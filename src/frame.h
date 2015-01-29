#ifndef FRAME_H_
#define FRAME_H_

#include <gromacs/fileio/xtcio.h>

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
    std::string resname;
    /** Atomtype as a string.  I don't want to be dealing with *char */
    std::string atom_type;
    /** Atomic coordinates in x, y, z */
    float coords[3];
    /** Atomic charge from the force field */
    float charge = 0.f;
    /** Atomic mass */
    float mass = 0.f;
    /** \brief Create an Atom, set its number and zero its coordinates */
    Atom(int num){atom_num = num; coords[0] = 0.f; coords[1] = 0.f; coords[2] = 0.f;};
    /** Create a blank Atom instance */
    Atom(){coords[0] = 0.f; coords[1] = 0.f; coords[2] = 0.f;};
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
    std::vector<std::string> atom_names;
    /** The number of atoms in the residue */
    int num_atoms;
    /** Constructor to set res_name */
    Residue(const std::string tmp){res_name = tmp;};
    /** Blank constructor */
    Residue(){};
};


/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
public:
    /** Has the Frame been properly setup yet? */
    bool isSetup_ = false;
    /** Has the XTC output been setup yet? */
    bool outputSetup_ = false;
    /** GROMACS xtc file to import frames */
    t_fileio *xtcInput_ = nullptr;
    /** GROMACS xtc file to export frames */
    t_fileio *xtcOutput_ = nullptr;
    /** The number of this Frame, starting at 0 */
    int num_ = 0;
    /** The simulation step corresponding to this frame */
    int step_ = 0;
    /** The number of atoms stored in this frame */
    int numAtoms_ = 0;
    /** The number of atoms stored in this frame that we find interesting */
    int numAtomsTrack_ = 0;
    /** The number of residues stored in this frame that we find interesting */
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** Vector of Residues; Each Residue contains pointers to atoms */
    std::vector<Residue> residues_;
    /** The simulation time of this frame, in picoseconds */
    float time_ = 0.f;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec_ = 0.f;
    /** Size of the simulation box */
    float box_[3][3];
    /** Holds atomic coordinates for GROMACS */
    rvec *x_ = NULL;
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_;
    /** Dictionary mapping atom numbers to atom names */
    std::map<int, std::string> numToName_;
    /** Dictionary mapping atom names to numbers */
    std::map<std::string, int> nameToNum_;
    /** Is frame invalid for some reason - molecule lies on pbc */
    bool invalid_ = false;

    /**
    * \brief Create Frame passing frame number, number of atoms to store and the frame name
    *
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocateAtoms()
    */
    Frame(int num, int natoms, std::string name);

    /**
    * \brief Create Frame passing config files.
    *
    * Replaces calls to the function Frame::setupFrame()
    */
    Frame(const std::string topname, const std::string xtcname);

    /**
    * \brief Create Frame by copying data from another Frame
    *
    * Intended for creating a CG Frame from an atomistic one.  Atoms are not copied.
    */
    Frame(const Frame &frame);

    /**
    * \brief Destructor to free memory allocated by GROMACS functions
    */
    ~Frame();

//    /**
//    * \brief Move constructor
//    */
//    Frame(Frame&& frame);

//    /** \brief Assignment operator
//    * Doesn't copy atoms.
//    */
//    Frame &operator=(const Frame &frame);

    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    * \throws runtime_error if Frame has already been setup
    *
    * GROMACS read_first_xtc() gets data from the XTC file about the system.
    * This function uses this data to create a Frame object to process this data
    */
    bool setupFrame(const std::string topname, t_fileio *xtc);

    /**
    * \brief Read next frame from the open XTC file
    */
    bool readNext();

    /**
    * \brief Allocate space for a number of atoms
    *
    * Used if the number of atoms isn't known at time of creation
    */
    int allocateAtoms(int);

    /**
    * \brief Prepare to write XTC output.
    * \throws std::runtime_error if memory cannot be allocated for atom array
    * \throws std::runtime_error if output TOP file cannot be opened
    * Allocate necessary atom array and create TOP file.
    */
    void setupOutput(const std::string xtcnameout, const std::string topnameout);

    /**
    * \brief Write Frame to XTC output file
    */
    bool writeToXtc();

    /** \brief Recentre simulation box on an atom
    * Avoids problems where a residue is split by the periodic boundary,
    * causing bond lengths to be calculated incorrectly
    */
    void recentreBox(const int atom_num);

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
    float bondAngle(int a, int b, int c, int d);


    /**
    * \brief Calculate angle or dihedral between atoms in a BondStruct object
    *
    * Wrapper around float bondAngle(int, int, int, int)
    */
    float bondAngle(BondStruct *bond);
};

#endif
