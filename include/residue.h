//
// Created by james on 07/05/15.
//

#ifndef CGTOOL_RESIDUE_H
#define CGTOOL_RESIDUE_H

#include <string>
#include <map>

class Residue{
public:
    /** \brief Name of this residue */
    std::string resname = "";
    /** \brief Name of reference atom - given by user */
    std::string ref_atom_name = "";
    /** \brief Number of reference atom */
    int ref_atom = -1;
    /** \brief Number of atoms in this residue */
    int num_atoms = -1;
    /** \brief Number of this residue present in simulation */
    int num_residues = -1;
    /** \brief Total number of atoms across all residues of this type */
    int total_atoms = -1;
    /** \brief Atom number in GRO where residue starts */
    int start = -1;
    /** \brief Atom number in GRO where residue ends */
    int end = -1;
    /** \brief Has data been read into this residue */
    bool populated = false;
    /** Map mapping atom names to numbers for each residue */
    std::map<std::string, int> name_to_num;

    /** \brief Blank constructor */
    Residue(){};
    /** \brief Blank destructor */
    ~Residue(){};

    /** \brief Copy constructor */
    Residue& operator=(const Residue other);

    /** \brief Calcaulate the total number of atoms */
    void calc_total();

    /** \brief Initialise num_atoms and num_residues to 0 */
    void init();

    /** \brief Set num_atoms */
    void set_num_atoms(const int val);

    /** \brief Set num_residues */
    void set_num_residues(const int val);

    /** \brief Set start */
    void set_start(const int val);

    /** \brief Set and check the residue name */
    void set_resname(const std::string &val);

    /** \brief Set total_atoms */
    void set_total_atoms(const int val);

    /** \brief Print information about this residue */
    void print(const bool extra=false) const;
};

#endif //CGTOOL_RESIDUE_H
