#ifndef FRAME_H_
#define FRAME_H_

#include "xdrfile.h"

#include <vector>
#include <string>
#include <map>

#include "bond_struct.h"
#include "residue.h"

/**
* \brief Struct to hold atom data
*/
struct Atom{
    /** Atomtype as a string.  I don't want to be dealing with *char */
    std::string atom_type = "";
    /** Atomic coordinates in x, y, z */
    double coords[3] = {0., 0., 0.};
    /** Atomic charge from the force field */
    double charge = 0.;
    /** Atomic mass */
    double mass = 0.;
    /** Residue number */
    int resnum = 0;
    /** Create a blank Atom instance */
    Atom(){};
};

enum class BoxType{CUBIC, TRICLINIC};

/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
protected:
    /** Has the XTC output been setup yet? */
    bool outputSetup_ = false;
    /** GROMACS xtc file to import frames */
    XDRFILE *xtcInput_ = nullptr;
    /** GROMACS xtc file to export frames */
    XDRFILE *xtcOutput_ = nullptr;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec_ = 0.f;
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_ = "";
    /** What box shape do we have?  Currently must be cubic */
    BoxType boxType_ = BoxType::CUBIC;

    void createAtoms(int natoms=-1);
    bool initFromXTC(const std::string &xtcname);
    void initFromGRO(const std::string &groname, std::vector<Residue> &residues);
    void copyCoordsIntoAtoms(int natoms=-1);

    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    * \throws logic_error if Frame has already been setup
    *
    * Uses libxdrfile to get number of atoms and allocate storage.
    * This function uses this data to create a Frame object to process this data.
    */
    void initFromITP(const std::string &topname);

    /** \brief Recentre simulation box on an atom
    * Avoids problems where a residue is split by the periodic boundary,
    * causing bond lengths to be calculated incorrectly */
    void recentreBox(const int atom_num);

public:
    /** Has the Frame been properly setup yet? */
    bool isSetup_ = false;
    /** The number of atoms stored in this frame that we find interesting */
    int numAtomsTrack_ = 0;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** The number of atoms stored in this frame */
    int numAtoms_ = 0;
    /** The simulation time of this frame, in picoseconds */
    float time_ = 0.f;
    /** The number of this Frame, starting at 0 */
    int num_ = 0;
    /** The simulation step corresponding to this frame */
    int step_ = 0;
    /** Map mapping atom names to numbers for each residue */
    std::map<std::string, int> nameToNum_;
    /** Size of the simulation box */
    float box_[3][3];
    std::vector<Residue> residues_;
    /** Holds atomic coordinates for GROMACS */
    rvec *x_ = nullptr;


    /** \brief Create Frame passing frame number, number of atoms to store and the frame name
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocateAtoms() */
    Frame(const int num, const int natoms, const std::string name);

    /** \brief Create Frame passing config files.
    * Replaces calls to the function Frame::setupFrame() */
    Frame(const std::string &itpname, const std::string &xtcname,
          const std::string &groname, std::vector<Residue> &residues,
          const bool do_itp=true);

    /** \brief Create Frame by copying data from another Frame
    * Intended for creating a CG Frame from an atomistic one.  Atoms are not copied. */
    Frame(const Frame &frame);

    /** \brief Destructor to free memory allocated by XDR functions */
    ~Frame();


    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    bool readNext();

    /**
    * \brief Prepare to write XTC output.
    * \throws std::runtime_error if memory cannot be allocated for atom array
    * \throws std::runtime_error if output TOP file cannot be opened
    * Allocate necessary atom array and create TOP file.
    */
    void setupOutput(std::string xtcnameout="", std::string top="");

    /** \brief Write Frame to XTC output file */
    bool writeToXtc();

    /** \brief Print info for all atoms up to n.  Default print all. */
    void printAtoms(int natoms=-1);

    /** \brief Print all atoms up to n to GRO file.  Default print all. */
    void printGRO(std::string filename="", int natoms=-1);

    /**
    * \brief Calculate distance between two atoms in a BondStruct object
    * Wrapper around float bondLength(int, int)
    */
    double bondLength(BondStruct &bond, const int offset=0);

    /**
    * \brief Calculate angle or dihedral between atoms in a BondStruct object
    * Wrapper around float bondAngle(int, int, int, int)
    */
    double bondAngle(BondStruct &bond, const int offset=0);

};

#endif
