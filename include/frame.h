#ifndef FRAME_H_
#define FRAME_H_

#include <vector>
#include <string>
#include <array>

#include "residue.h"

class TrjOutput;
class TrjInput;


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
    /** Atomic coordinates in x, y, z */
    std::array<double, 3> coords = {{0., 0., 0.}};
    /** Atom dipole components in x, y, z and magnitude */
    std::array<double, 4> dipole = {{0., 0., 0., 0.}};

    /** Atomic charge from the force field */
    double charge = 0.;
    /** Atomic mass */
    double mass = 0.;

    /** Atomtype as a string */
    std::string atom_type = "\0\0\0\0\0";
    /** Atomname as a string */
    std::string atom_name = "\0\0\0\0\0";

    /** Lennard-Jones C6 parameter */
    double c06 = 0.;
    /** Lennard-Jones C12 parameter */
    double c12 = 0.;

    /** Residue number */
    int resnum = 0;
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
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_ = "";
    /** What box shape do we have?  Currently should be cubic */
    BoxType boxType_ = BoxType::CUBIC;

    /** \brief Input readers */
    TrjInput *trjIn_ = nullptr;

    bool initFromGRO(const std::string &groname);

    /** \brief Perform wraparound to put all atoms in box.
     * Equivalent to GROMACS trjconv -pbc atom */
    void pbcAtom(int natoms=-1);

    /** \brief Write the minimal GROMACS TOP file */
    void writeTOP(const std::string &filename);

public:
    /** Has the Frame been properly setup yet? */
    bool isSetup_ = false;
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
    /** Size of the simulation box */
    float box_[3][3];
    std::array<double, 3> boxDiag_;
    /** Which data have been loaded into atoms? */
    AtomsHave atomHas_;
    std::vector<Residue> &residues_;

    /** \brief Create Frame passing config files.
    * Replaces calls to the function Frame::setupFrame() */
    Frame(const std::string &xtcname, const std::string &groname,
          std::vector<Residue> &residues);

    Frame(const Frame &frame) : name_(frame.name_), boxType_(frame.boxType_),
                                time_(frame.time_), num_(frame.num_), step_(frame.step_),
                                residues_(frame.residues_){}

    Frame(const Frame &frame, std::vector<Residue> &residues) :
                                name_(frame.name_), boxType_(frame.boxType_),
                                time_(frame.time_), num_(frame.num_), step_(frame.step_),
                                residues_(frame.residues_){}

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

    /** \brief Output frame to trajectory file. */
    bool outputTrajectoryFrame(TrjOutput &output);

    /** \brief Print info for all atoms up to n.  Default print all. */
    void printAtoms(int natoms=-1) const;

    /** \brief Print box vectors */
    void printBox() const;
};

#endif
