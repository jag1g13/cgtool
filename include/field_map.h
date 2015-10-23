#ifndef FIELDMAP_H_
#define FIELDMAP_H_

#include <vector>
#include <array>

#include "light_array.h"
#include "frame.h"
#include "cg_map.h"

namespace constants{
    const double ENM2DEBYE = 2.0819434;
}

/**
* \brief Class to hold and operate on electric field maps
*
* Calculates an electric field from point charges in a molecule,
* then fits a set of dipoles to this
*/
class FieldMap{
private:
    /** \brief Residues of atomistic representation */
    const std::vector<Residue> &aaResidues_;
    /** \brief Residues of coarse-grained representation */
    const std::vector<Residue> &cgResidues_;
    /** \brief Number of atoms in a single aa residue */
    const int aaNumAtoms_ = 0;
    /** \brief Number of atoms in a single cg residue */
    const int cgNumAtoms_ = 0;

    /** Dimensions of the field grids (3 ints) */
    int gridDims_[3];
    /** \brief Border to leave around molecule in field grid
    * Also the radius of selection for the CHELPG style grid */
    const double border_ = 1.;    // 1nm

    /** Coordinates of each grid point */
    std::vector<std::array<double, 3>> gridContracted_;
    int numGridPoints_;
    /** Frame number */
    int frameNum_ = 0;

    std::vector<double> fieldMonopoleContracted_;
    std::vector<double> fieldDipoleContracted_;
    /** Dipole of each atom, coords, vector, magnitude */
    LightArray<double> dipoles_;
    double totalDipole_[6];
    double sumDipole_[6];

    /** Create a CHELPG style grid using only points in a shell around the molecule */
    void setupGridContracted(const Frame &frame);

    /** Calculate the electric field from point charges using the CHELPG style grid */
    void calcFieldMonopolesContracted(const Frame &frame);

    /** Calculate the electric field from point dipoles using the CHELPG style grid */
    void calcFieldDipolesContracted(const Frame &frame);

    /** \brief Calculate dipoles on beads directly from frame
    * Modifies charges in atomistic Frame */
    void calcDipolesDirect(const CGMap &cgmap, const Frame &cg_frame, const Frame &aa_frame);

    /** Calculate the total dipole of the system/residue */
    void calcTotalDipole(const Frame &aa_frame);

    /** Calculate the sum of calculated dipoles */
    void calcSumDipole();

    /** Print the dipole array */
    void printDipoles();

    /** Print the electric fields */
    void printFields();


public:
    /** Constructor for FieldMap to perform setup - uses init() */
    FieldMap(const int grid, const std::vector<Residue> &aa_res, const std::vector<Residue> &cg_res) :
            aaResidues_(aa_res), cgResidues_(cg_res),
            aaNumAtoms_(aa_res[0].num_atoms), cgNumAtoms_(cg_res[0].num_atoms){
        gridDims_[0] = grid; gridDims_[1] = grid; gridDims_[2] = grid;
        dipoles_.alloc(cgNumAtoms_, 6);
    }

    /** \brief Run all electric field calculations */
    void calculate(const Frame &aa_frame, const Frame &cg_frame, const CGMap &cgmap);

    /** Print the electric field to file */
    void printFieldsToFile();
};

#endif