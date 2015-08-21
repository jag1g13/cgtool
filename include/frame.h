#ifndef FRAME_H_
#define FRAME_H_

#include "xdrfile.h"

#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

#include "residue.h"
//#include "trj_output.h"
class TrjOutput;

/** \brief Contains data from a line of a GRO file */
struct GROLine{
    int resnum = 0;
    std::string resname = "";
    std::string atomname = "";
    int atomnum = 0;
    float coords[3] = {0.f, 0.f, 0.f};
    float velocity[3] = {0.f, 0.f, 0.f};

    /** \brief Populate from string - line of GRO file */
    int populate(const std::string &line){
        if(line.size() < 41) return 0;

        // Adjust for 7.2f formatting
        int float_len = 8;
        if(line.size() == 41) float_len = 7;

        resnum = std::stoi(line.substr(0, 5));
        resname = line.substr(5, 5);
        boost::trim(resname);
        atomname = line.substr(10, 5);
        boost::trim(atomname);
        atomnum = std::stoi(line.substr(15, 5));
        coords[0] = std::stof(line.substr(20, float_len));
        coords[1] = std::stof(line.substr(20 + float_len, float_len));
        coords[2] = std::stof(line.substr(20 + 2*float_len, float_len));

        if(line.size() >= 68){
            velocity[0] = std::stof(line.substr(44, float_len));
            velocity[1] = std::stof(line.substr(52, float_len));
            velocity[2] = std::stof(line.substr(60, float_len));
            return 2;
        }
        return 1;
    };

    /** \brief Print the parsed line for debugging */
    void print(){
        printf("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
               resnum, resname.c_str(), atomname.c_str(),
               atomnum, coords[0], coords[1], coords[2]);
    }
};

/** \brief Struct to keep track of which data have been loaded into atoms */
struct AtomsHave{
    /** \brief Atoms have been created */
    bool created = false;
    /** \brief Atoms have been assigned a type */
    bool atom_type = false;
    /** \brief Atoms have been assigned a name */
    bool atom_name = false;
    /** \brief Atoms have been assigned coordinates */
    bool coords = false;
    /** \brief Atoms have been assigned a charge */
    bool charge = false;
    /** \brief Atoms have been assigned a mass */
    bool mass = false;
    /** \brief Atoms have been assigned a residue number */
    bool resnum = false;
    /** \brief Atoms have been assigned Lennard Jones parameters */
    bool lj = false;
};

/** \brief Struct to hold atom data */
struct Atom{
    /** Atomtype as a string */
    std::string atom_type = "";
    /** Atomname as a string */
    std::string atom_name = "";
    /** Atomic coordinates in x, y, z */
    double coords[3] = {0., 0., 0.};
    /** Atomic charge from the force field */
    double charge = 0.;
    /** Atomic mass */
    double mass = 0.;
    /** Residue number */
    int resnum = 0;
    /** Lennard-Jones C6 parameter */
    double c06 = 0.;
    /** Lennard-Jones C12 parameter */
    double c12 = 0.;
    /** Create a blank Atom instance */
    Atom(){};

    Atom(const Atom &other){};
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
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_ = "";
    /** What box shape do we have?  Currently must be cubic */
    BoxType boxType_ = BoxType::CUBIC;

    /** \brief Output writer */
    TrjOutput *trjOut_ = nullptr;

    void createAtoms(int natoms=-1);

    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    * \throws logic_error if Frame has already been setup
    *
    * Uses libxdrfile to get number of atoms and allocate storage.
    * This function uses this data to create a Frame object to process this data.
    */
    bool initFromXTC(const std::string &xtcname);
    bool initFromGRO(const std::string &groname);
    void copyCoordsIntoAtoms(int natoms=-1);

    /** \brief Perform wraparound to put all atoms in box.
     * Equivalent to GROMACS trjconv -pbc atom */
    void pbcAtom(int natoms=-1);

public:
    /** Has the Frame been properly setup yet? */
    bool isSetup_ = false;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** The number of atoms stored in this frame */
    int numAtoms_ = 0;
    /** The simulation time of this frame, in picoseconds */
    float time_ = 0.f;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec_ = 0.f;
    /** The number of this Frame, starting at 0 */
    int num_ = 0;
    /** The simulation step corresponding to this frame */
    int step_ = 0;
    /** Size of the simulation box */
    float box_[3][3];
    /** Holds atomic coordinates for GROMACS */
    rvec *x_ = nullptr;
    /** Which data have been loaded into atoms? */
    AtomsHave atomHas_;
    std::vector<Residue> *residues_;

    /** \brief Create Frame passing config files.
    * Replaces calls to the function Frame::setupFrame() */
    Frame(const std::string &xtcname, const std::string &groname,
          std::vector<Residue> *residues);

    /** \brief Create Frame by copying data from another Frame
    * Intended for creating a CG Frame from an atomistic one.  Atoms are not copied. */
    Frame(const Frame &frame, std::vector<Residue> *residues=nullptr);

    /** \brief Destructor to free memory allocated by XDR functions */
    ~Frame();


    /**
    * \brief Read a frame from the XTC file into an existing Frame object
    *
    * Reads a frame into a pre-setup Frame object.
    * The same Frame object should be used for each frame to save time in allocation.
    */
    bool readNext();

    void initFromITP(const std::string &topname);
    void initFromFLD(const std::string &fldname);

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
    void printAtoms(int natoms=-1) const;

    /** \brief Print all atoms up to n to GRO file.  Default print all. */
    void printGRO(std::string filename="", int natoms=-1) const;

    /** \brief Print box vectors */
    void printBox() const;
};

#endif
