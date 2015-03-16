#ifndef FIELDMAP_H_
#define FIELDMAP_H_

#include <vector>

#include "array.h"
#include "frame.h"
#include "cg_map.h"

namespace constants{
    const float ENM2DEBYE = 2.0819434f;
}

//TODO tidy up this monster - partially done
/**
* \brief Class to hold and operate on electric field maps
*
* Calculates an electric field from point charges in a molecule,
* then fits a set of dipoles to this
*/
class FieldMap{
private:
    /** Dimensions of the field grids (3 ints) */
    int gridDims_[3];
    /** Array to hold electric field calculated from monopoles */
    Array fieldMonopole_;
    /** Array to hold electric field calculated from dipoles */
    Array fieldDipole_;
    /** \brief Border to leave around molecule in field grid
    * Also the radius of selection for the CHELPG style grid */
    double border_ = 2.f;    // 2nm
    /** Array to hold grid bounds; needs to be reset each frame (or often) */
    Array gridBounds_;
    /** Coordinates of each grid point */
    Array gridCoords_;
    /** Centre of grid */
    double gridCentre_[3];
    Array gridContracted_;
    int numGridPoints_;
    std::vector<double> fieldMonopoleContracted_;
    std::vector<double> fieldDipoleContracted_;
    int numDipoles_;
    /** Dipole of each atom, coords, vector, magnitude */
    Array dipoles_;
    Array totalDipole_;
    Array sumDipoles_;
    /** Frame number */
    int frameNum_ = 0;
    /** Number of atoms in a single aa residue */
    int aaNumAtoms_ = 0;
    /** Number of atoms in a single cg residue */
    int cgNumAtoms_ = 0;

    /** Determine grid bounds from a Frame object and do setup each time */
    void setupGrid(const Frame &frame);

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

    /** Print the electric fields to files */
    void printFieldsToFile();

public:
    /** Constructor for a blank instance of an electric field map */
    FieldMap();

    /** Constructor for FieldMap to perform setup */
    FieldMap(const int a, const int b, const int c, const int natoms=0);

    /** \brief Run all electric field calculations */
    void calculate(const Frame &aa_frame, const Frame &cg_frame, const CGMap &cgmap);

};

/** \brief Calculate the square of the distance between two points */
inline double distSqr(const double c1[3], const double c2[3]);
/** 3d vector dot product */
inline double dot(const double A[3], const double B[3]);
/** 3d vector magnitude */
inline double abs(const double vec[3]);
/** 3d polar coordinate conversion */
void polar(const double cart[3], double polar[3]);

#endif