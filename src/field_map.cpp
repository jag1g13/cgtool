#include "field_map.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "small_functions.h"

using std::string;
using std::min;
using std::max;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;

FieldMap::FieldMap(){
}

FieldMap::FieldMap(const int a, const int b, const int c, const int ndipoles){
    init(a, b, c, ndipoles);
}

FieldMap::FieldMap(const int res, const vector<Residue> *aa_residues, const vector<Residue> *cg_residues){
    aaResidues_ = aa_residues;
    cgResidues_ = cg_residues;
    init(res, res, res, (*cgResidues_)[0].num_atoms);
}

void FieldMap::init(const int a, const int b, const int c, const int ndipoles){
    gridDims_[0] = a; gridDims_[1] = b; gridDims_[2] = c;

    gridBounds_.init(3, 2, 1, false);
    gridCoords_.init(3, max(a, max(b, c)), 1, false);
    gridContracted_.init(a*b*c, 4, 1, true);

    fieldMonopole_.init(a, b, c, false);
    fieldDipole_.init(a, b, c, false);

    numDipoles_ = ndipoles;
    dipoles_.init(ndipoles, 6, 1, false);
    totalDipole_.init(6, 1, 1, false);
    sumDipoles_.init(6, 1, 1, false);
}

void FieldMap::calculate(const Frame &aa_frame, const Frame &cg_frame, const CGMap &cgmap){
    frameNum_ = aa_frame.num_;
    if(frameNum_ != cg_frame.num_) throw std::logic_error("Frame numbers do not match");
    // Calculate only the first molecule
    aaNumAtoms_ = (*aaResidues_)[0].num_atoms;
    cgNumAtoms_ = (*cgResidues_)[0].num_atoms;

    calcDipolesDirect(cgmap, cg_frame, aa_frame);
    setupGrid(aa_frame);
    setupGridContracted(aa_frame);
    calcFieldMonopolesContracted(aa_frame);
    calcFieldDipolesContracted(cg_frame);
    calcTotalDipole(aa_frame);
    calcSumDipole();
    printDipoles();

    StatsBox sb = vector_stats(fieldMonopoleContracted_, fieldDipoleContracted_);
    cout << "\tRMS: " << sb.rmsd << "\tRRMS: " << sb.nrmsd << endl;
}

void FieldMap::setupGrid(const Frame &frame){
    // Create min and max initial values
    gridBounds_(0, 0) = frame.atoms_[0].coords[0]; gridBounds_(0, 1) = frame.atoms_[0].coords[0];
    gridBounds_(1, 0) = frame.atoms_[0].coords[1]; gridBounds_(1, 1) = frame.atoms_[0].coords[1];
    gridBounds_(2, 0) = frame.atoms_[0].coords[2]; gridBounds_(2, 1) = frame.atoms_[0].coords[2];

    for(int i=0; i < aaNumAtoms_; i++){
        // Find bounding box of molecule
        const Atom *atom = &(frame.atoms_[i]);
        gridBounds_(0, 0) = min(gridBounds_(0, 0), atom->coords[0]);
        gridBounds_(0, 1) = max(gridBounds_(0, 1), atom->coords[0]);
        gridBounds_(1, 0) = min(gridBounds_(1, 0), atom->coords[1]);
        gridBounds_(1, 1) = max(gridBounds_(1, 1), atom->coords[1]);
        gridBounds_(2, 0) = min(gridBounds_(2, 0), atom->coords[2]);
        gridBounds_(2, 1) = max(gridBounds_(2, 1), atom->coords[2]);
    }

    gridBounds_(0, 0) -= border_; gridBounds_(0, 1) += border_;
    gridBounds_(1, 0) -= border_; gridBounds_(1, 1) += border_;
    gridBounds_(2, 0) -= border_; gridBounds_(2, 1) += border_;

    gridCentre_[0] = (gridBounds_(0, 1) - gridBounds_(0, 0)) / 2.;
    gridCentre_[1] = (gridBounds_(1, 1) - gridBounds_(1, 0)) / 2.;
    gridCentre_[2] = (gridBounds_(2, 1) - gridBounds_(2, 0)) / 2.;

    gridCoords_.linspace(0, gridDims_[0], gridBounds_(0, 0), gridBounds_(0, 1));
    gridCoords_.linspace(1, gridDims_[1], gridBounds_(1, 0), gridBounds_(1, 1));
    gridCoords_.linspace(2, gridDims_[2], gridBounds_(2, 0), gridBounds_(2, 1));
}

void FieldMap::setupGridContracted(const Frame &frame){
    double radmin2, dist2;
    // Use an rmax equal to border_ around molecule
    const double rmax2 = border_ * border_;
    // Reject if within 1A of an atom (inside atomic radius)
    const double vrad2 = 0.1 * 0.1;
    bool accepted;
    int accepted_count = 0, close_count = 0, far_count = 0;
    double coords[3];

    gridContracted_.resetAppendedRows();

//    #pragma omp parallel for private(coords, radmin2, dist2, accepted) reduction(+: accepted_count, far_count, close_count)
    for(int i=0; i < gridDims_[0]; i++){
        coords[0] = gridCoords_(0, i);
        for(int j=0; j < gridDims_[1]; j++){
            coords[1] = gridCoords_(1, j);
            for(int k=0; k < gridDims_[2]; k++){
                coords[2] = gridCoords_(2, k);
                radmin2 = 500.;
                accepted = true;

                for(int ii=0; ii < aaNumAtoms_; ii++){
                    dist2 = distSqr(frame.atoms_[ii].coords, coords);
                    radmin2 = min(radmin2, dist2);
                    if(dist2 < vrad2){
                        accepted = false;
                        ++close_count;
                        break;
                    }
                }

                if(accepted && radmin2 > rmax2){
                    ++far_count;
                    continue;
                }

                if(accepted){
//                    #pragma omp critical
                    {
                        gridContracted_.append(coords);
                        ++accepted_count;
                    }
                }
            }
        }
    }
    numGridPoints_ = gridContracted_.getAppendedRows();
    fieldMonopoleContracted_.resize(numGridPoints_);
    fieldDipoleContracted_.resize(numGridPoints_);
}

void FieldMap::calcFieldMonopolesContracted(const Frame &frame){
    double coords[3];
//#pragma omp parallel for private(coords)
    for(int i=0; i < numGridPoints_; i++){
        fieldMonopoleContracted_[i] = 0.;
        for(int j=0; j<aaNumAtoms_; j++){
            coords[0] = gridContracted_(i, 0);
            coords[1] = gridContracted_(i, 1);
            coords[2] = gridContracted_(i, 2);
            fieldMonopoleContracted_[i] += frame.atoms_[j].charge /
                    distSqr(frame.atoms_[j].coords, coords);
        }
    }
}


void FieldMap::printFields(){
    cout << "ELECTRIC FIELDS" << endl;
    for(int i=0; i<numGridPoints_; i++){
        cout << fieldMonopoleContracted_[i] << "\t" << fieldDipoleContracted_[i] << endl;
    }
}

void FieldMap::calcFieldDipolesContracted(const Frame &frame){
    double vec_a[3], vec_b[3];
    double abs_a;
    double coords[3];
//#pragma omp parallel for private(coords, vec_a, vec_b, abs_a)
    for(int i=0; i < numGridPoints_; i++) {
        fieldDipoleContracted_[i] = 0.;
        // For charge on the cg bead
        for(int j=0; j < cgNumAtoms_; j++){
//            fieldDipoleContracted_[i] += dipoles_(j, 5) / (abs_a*abs_a);

            for(int k=0; k<3; k++) {
                coords[k] = gridContracted_(i, k);
                vec_a[k] = coords[k] - frame.atoms_[j].coords[k];
                vec_b[k] = dipoles_(j, k);
            }

            abs_a = abs(vec_a);
            double cos_dip_angle = dot(vec_a, vec_b) / (dipoles_(j, 5) * abs_a);

            // If it's NaN, there's no dipole - don't add anything to the field
            if(cos_dip_angle != cos_dip_angle) cos_dip_angle = 0.;

            fieldDipoleContracted_[i] += dipoles_(j, 5) * cos_dip_angle / (abs_a*abs_a);
            // Do I need to include the field from the charge on the bead?
//            fieldDipoleContracted_[i] += frame.atoms_[j].charge /
//                    distSqr(frame.atoms_[j].coords, coords);
        }
    }
}


/**
* This uses a hack to allow direct calculation: charges within a bead are
* rescaled to make each bead neutral.  This means that dipoles are no longer
* dependent on the frame of reference.
* The electric field from these dipoles should be compared against the field
* from atomic point charges to determine validity.  They may need to be rescaled.
* I don't see a better way to do this.
*/
void FieldMap::calcDipolesDirect(const CGMap &cgmap, const Frame &cg_frame, const Frame &aa_frame){
    dipoles_.zero();
    for(int i=0; i<cgNumAtoms_; i++){
        // For each bead in the CG frame
        const BeadMap &bead_type = cgmap.mapping_[i];
        const Atom &cg_atom = cg_frame.atoms_[i];

        for(const int &j : bead_type.atom_nums){
            // For each atom inside the bead
            double charge = aa_frame.atoms_[j].charge;
            // Rescale charges so bead charge is zero
            // This is how GMX_DIPOLE does it
            charge -= cg_atom.charge * aa_frame.atoms_[j].mass / bead_type.mass;
            dipoles_(i, 0) += aa_frame.atoms_[j].coords[0] * charge;
            dipoles_(i, 1) += aa_frame.atoms_[j].coords[1] * charge;
            dipoles_(i, 2) += aa_frame.atoms_[j].coords[2] * charge;
        }
        double charge = cg_frame.atoms_[i].charge;
        dipoles_(i, 0) -= cg_frame.atoms_[i].coords[0] * charge;
        dipoles_(i, 1) -= cg_frame.atoms_[i].coords[1] * charge;
        dipoles_(i, 2) -= cg_frame.atoms_[i].coords[2] * charge;

        // Calculate magnitude here
        dipoles_(i, 5) = sqrt(dipoles_(i, 0)*dipoles_(i, 0) +
                              dipoles_(i, 1)*dipoles_(i, 1) +
                              dipoles_(i, 2)*dipoles_(i, 2));
    }
}

void FieldMap::calcTotalDipole(const Frame &aa_frame){
    totalDipole_.zero();

    for(int i=0; i < aaNumAtoms_; i++){
        double charge = aa_frame.atoms_[i].charge;
        totalDipole_(0) += aa_frame.atoms_[i].coords[0] * charge;
        totalDipole_(1) += aa_frame.atoms_[i].coords[1] * charge;
        totalDipole_(2) += aa_frame.atoms_[i].coords[2] * charge;
    }

    totalDipole_(5) = sqrt(totalDipole_(0)*totalDipole_(0) +
                           totalDipole_(1)*totalDipole_(1) +
                           totalDipole_(2)*totalDipole_(2));
}

void FieldMap::calcSumDipole(){
    sumDipoles_.zero();
    for(int i=0; i < numDipoles_; i++){
        sumDipoles_(0) += dipoles_(i, 0);
        sumDipoles_(1) += dipoles_(i, 1);
        sumDipoles_(2) += dipoles_(i, 2);
    }

    // Calc magnitude
    sumDipoles_(5) = sqrt(sumDipoles_(0)*sumDipoles_(0) +
                          sumDipoles_(1)*sumDipoles_(1) +
                          sumDipoles_(2)*sumDipoles_(2));
}

void FieldMap::printDipoles(){
    cout << "Total molecular dipole and sum of bead dipoles" << endl;
    totalDipole_.print(8, 4, constants::ENM2DEBYE);
    sumDipoles_.print(8, 4, constants::ENM2DEBYE);
}

void FieldMap::printFieldsToFile(){
    const string aa_file = "field" + std::to_string(frameNum_) + "_aa.csv";
    const string cg_file = "field" + std::to_string(frameNum_) + "_cg.csv";
    FILE *faa = std::fopen(aa_file.c_str(), "w");
    if(faa == NULL) throw std::runtime_error("Could not open file to write electric field");
    FILE *fcg = std::fopen(cg_file.c_str(), "w");
    if(fcg == NULL) throw std::runtime_error("Could not open file to write electric field");

    for(int i=0; i<numGridPoints_; i++){
        fprintf(faa, "%12.5f%12.5f%12.5f%12.5f\n",
                gridContracted_(i, 0), gridContracted_(i, 1),
                gridContracted_(i, 2), fieldMonopoleContracted_[i]);
        fprintf(fcg, "%12.5f%12.5f%12.5f%12.5f\n",
                gridContracted_(i, 0), gridContracted_(i, 1),
                gridContracted_(i, 2), fieldDipoleContracted_[i]);
    }

    std::fclose(faa);
    std::fclose(fcg);
    cout << aa_file << endl << cg_file << endl;
}

