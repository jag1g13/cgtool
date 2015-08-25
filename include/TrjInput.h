//
// Created by james on 24/08/15.
//

#ifndef CGTOOL_TRJINPUT_H
#define CGTOOL_TRJINPUT_H

#include <string>

#include "frame.h"
#include "residue.h"

struct TrjHas{
    bool num_atoms = false;
    bool atom_num = false;
    bool atom_name = false;
    bool atom_coord = false;
    bool atom_charge = false;
    bool atom_lj = false;
    bool atom_mass = false;

    bool num_residues = false;
    bool residue_num = false;
    bool residue_name = false;
    bool residue_num_atoms = false;
};

class TrjInput{
protected:
    /** \brief Number of atoms present. */
    int natoms_;

    /** \brief Open and prepare input file. */
    virtual int openFile(const std::string &filename) = 0;
    /** \brief Close input file. */
    virtual int closeFile() = 0;

public:
    /** \brief Read a frame from input file.  Pure virtual function. */
    virtual int readFrame(Frame &frame) = 0;

    /** \brief Empty destructor to be overwritten. */
    virtual ~TrjInput(){};

    /** \brief Read in residues from input file.  Does not have to be supported. */
    virtual void readResidues(std::vector<Residue> &residues){
        throw std::logic_error("Input file does not support reading residues");
    };

    int getNumAtoms() const{
        return natoms_;
    }

    friend class Frame;
};

#endif //CGTOOL_TRJINPUT_H
