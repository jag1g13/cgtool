#include "field_map.h"

#include <cmath>
#include <iostream>

using std::min;
using std::max;
using std::vector;
using std::cout;
using std::endl;

FieldMap::FieldMap(){
}

FieldMap::FieldMap(const int a, const int b, const int c, const int ndipoles){
    gridDims_.reserve(3);
    gridDims_[0] = a; gridDims_[1] = b; gridDims_[2] = c;
    gridCentre_.reserve(3);
    fieldMonopole_.init(a, b, c, false);
    fieldDipole_.init(a, b, c, false);
    gridBounds_.init(3, 2, 1, false);
    gridCoords_.init(3, max(a, max(b, c)), 1, false);
    numDipoles_ = ndipoles;
    //TODO init this later
    dipoles_.init(ndipoles, 6, 1, false);
    //TODO try out Boost arrays, are they faster/better?
    gridContracted_.init(a*b*c, 4, 1, true);
    totalDipole_.init(6, 1, 1, false);
    sumDipoles_.init(6, 1, 1, false);
}

void FieldMap::setupGrid(const Frame &frame){
    // create min and max initial values
    gridBounds_(0, 0) = 1e6f; gridBounds_(0, 1) = -1e6f;
    gridBounds_(1, 0) = 1e6f; gridBounds_(1, 1) = -1e6f;
    gridBounds_(2, 0) = 1e6f; gridBounds_(2, 1) = -1e6f;
//    for(int ii : frame->residues_[0].atoms){
    for(int i=0; i < frame.numAtomsTrack_; i++){
        // find bounding box of molecule
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

    gridCentre_[0] = (gridBounds_(0, 1) - gridBounds_(0, 0)) / 2.f;
    gridCentre_[1] = (gridBounds_(1, 1) - gridBounds_(1, 0)) / 2.f;
    gridCentre_[2] = (gridBounds_(2, 1) - gridBounds_(2, 0)) / 2.f;

    gridCoords_.linspace(0, gridDims_[0], gridBounds_(0, 0), gridBounds_(0, 1));
    gridCoords_.linspace(1, gridDims_[1], gridBounds_(1, 0), gridBounds_(1, 1));
    gridCoords_.linspace(2, gridDims_[2], gridBounds_(2, 0), gridBounds_(2, 1));
}

void FieldMap::setupGridContracted(const Frame &frame){
    float radmin, dist;
    float rmax = border_;    // use an rmax equal to border_ around molecule
    float vrad = 0.1f;       // reject if within 1A of an atom (inside atomic radius)
    bool accepted;
    int accepted_count = 0, close_count = 0, far_count = 0;
    vector<float> coords(3);
    gridContracted_.appendedRows_ = 0;

#pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        coords[0] = gridCoords_(0, i);
        for(int j=0; j < gridDims_[1]; j++){
            coords[1] = gridCoords_(1, j);
            for(int k=0; k < gridDims_[2]; k++){
                coords[2] = gridCoords_(2, k);
                radmin = 500.f;
                accepted = true;
                for(int ii=0; ii < frame.numAtomsTrack_; ii++){
                    dist = float(sqrt(distSqr(frame.atoms_[ii].coords, coords[0], coords[1], coords[2])));
                    radmin = min(radmin, dist);
                    if(dist < vrad){
                        accepted = false;
                        ++close_count;
                        break;
                    }
                }
                if(accepted && radmin > rmax){
                    ++far_count;
                    continue;
                }
                if(accepted){
                    gridContracted_.append(coords);
                    ++accepted_count;
                }
            }
        }
    }
    numGridPoints_ = gridContracted_.appendedRows_;
    fieldMonopoleContracted_.resize(numGridPoints_);
    fieldDipoleContracted_.resize(numGridPoints_);
}

void FieldMap::calcFieldMonopolesContracted(const Frame &frame){
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        fieldMonopoleContracted_[i] = 0.f;
        for(Atom atom : frame.atoms_) {
            fieldMonopoleContracted_[i] += atom.charge /
                    distSqr(atom.coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));
        }
    }
}

inline float dot(const float *A, const float* B){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

inline float abs(const float* vec){
    return float(sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]));
}

void FieldMap::printFields(){
    cout << "ELECTRIC FIELDS" << endl;
    for(int i=0; i<numGridPoints_; i++){
        cout << fieldMonopoleContracted_[i] << "\t" << fieldDipoleContracted_[i] << endl;
    }
}

//TODO dipole field - to check that they're right - pretty much done, needs testing
void FieldMap::calcFieldDipolesContracted(const Frame &frame){
    float vec_a[3], vec_b[3];
    float abs_a;
//    cout << numGridPoints_ << endl;
//    cout << frame->numAtomsTrack_ << endl;
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        fieldDipoleContracted_[i] = 0.f;
        // for charge on the cg bead
//        fieldDipoleContracted_[i] += dipoles_(j, 5) / (abs_a*abs_a);
        for(int j=0; j < frame.numAtomsTrack_; j++){
//            cout << "j=" << j << endl;
            for(int k=0; k<3; k++) {
                vec_a[k] = gridContracted_(i, k) - frame.atoms_[j].coords[k];
                vec_b[k] = dipoles_(j, k);
            }
            abs_a = abs(vec_a);
            float cos_dip_angle = dot(vec_a, vec_b) / (dipoles_(j, 5) * abs_a);
            // If it's NaN, there's no dipole - don't add anything to the field
            if(cos_dip_angle != cos_dip_angle) cos_dip_angle = 0.f;
            fieldDipoleContracted_[i] += dipoles_(j, 5) * cos_dip_angle / (abs_a*abs_a);
            // Do I need to include the field from the charge on the bead?
//            fieldDipoleContracted_[i] += frame->atoms_[j].charge /
//                    distSqr(frame->atoms_[j].coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));
        }
    }

    #pragma omp master
    {
//        StatsBox sb = vector_stats(&fieldMonopoleContracted_, &fieldDipoleContracted_);
//        cout << "\tRMS: " << sb.rms << "\tRRMS: " << sb.rrms << endl;
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
void FieldMap::calcDipolesDirect(const CGMap &cgmap, const Frame &cg_frame, Frame &aa_frame){
    dipoles_.zero();
    for(int i=0; i<cgmap.num_beads; i++){
        // for each bead in the CG frame
        const BeadMap &bead_type = cgmap.mapping_[i];
        const Atom &cg_atom = cg_frame.atoms_[i];

        for(const int &j : bead_type.atom_nums){
            // for each atom inside the bead
            float charge = aa_frame.atoms_[j].charge;
            // rescale charges so bead charge is zero
//            charge -= cg_atom.charge / bead_type.num_atoms;
            // this is how GMX_DIPOLE does it
            charge -= cg_atom.charge * aa_frame.atoms_[j].mass / bead_type.mass;
            dipoles_(i, 0) += aa_frame.atoms_[j].coords[0] * charge;
            dipoles_(i, 1) += aa_frame.atoms_[j].coords[1] * charge;
            dipoles_(i, 2) += aa_frame.atoms_[j].coords[2] * charge;
        }
        float charge = cg_frame.atoms_[i].charge;
        dipoles_(i, 0) -= cg_frame.atoms_[i].coords[0] * charge;
        dipoles_(i, 1) -= cg_frame.atoms_[i].coords[1] * charge;
        dipoles_(i, 2) -= cg_frame.atoms_[i].coords[2] * charge;

        // calculate magnitude here
        dipoles_(i, 5) = float(sqrt(dipoles_(i, 0)*dipoles_(i, 0) +
                              dipoles_(i, 1)*dipoles_(i, 1) +
                              dipoles_(i, 2)*dipoles_(i, 2)));
    }
    printDipoles();
}

void FieldMap::calcDipolesFit(const CGMap &cgmap, const Frame &cg_frame, const Frame &aa_frame){
    dipoles_.zero();
    // keep track of which dipoles we know
    vector<int> dipoles_calculated;
    // calculate dipoles on all uncharged beads
    for(int i = 0; i < cgmap.num_beads; i++){
        const BeadMap &bead_type = cgmap.mapping_[i];
        // skip beads that are charged
        if(bead_type.charge > 0.01f || bead_type.charge < -0.01f) continue;

        const Atom &cg_atom = cg_frame.atoms_[i];
        // for each bead in the CG frame
        for(const int &j : bead_type.atom_nums){
            // for each atom inside the bead
            float charge = aa_frame.atoms_[j].charge;
            dipoles_(i, 0) += aa_frame.atoms_[j].coords[0] * charge;
            dipoles_(i, 1) += aa_frame.atoms_[j].coords[1] * charge;
            dipoles_(i, 2) += aa_frame.atoms_[j].coords[2] * charge;
        }
        dipoles_calculated.push_back(i);
    }
    int num_remaining = cgmap.num_beads - int(dipoles_calculated.size());

    // calculate residual molecular dipole
    calcTotalDipole(aa_frame);
    calcSumDipole(dipoles_calculated);
    // keep residual dipole in place of totalDipole_
    totalDipole_ -= sumDipoles_;

    // divide residual between remaining beads
    for(int i = 0; i < cgmap.num_beads; i++){
        const BeadMap &bead_type = cgmap.mapping_[i];
        // skip beads that are NOT charged
        if(bead_type.charge < 0.01f && bead_type.charge > -0.01f) continue;

        for(int j = 0; j < 3; j++){
            dipoles_(i, j) = totalDipole_(j) / num_remaining;
        }
    }

    // calculate magnitudes
    for(int i = 0; i < cgmap.num_beads; i++){
        dipoles_(i, 5) = float(sqrt(dipoles_(i, 0) * dipoles_(i, 0) +
                dipoles_(i, 1) * dipoles_(i, 1) +
                dipoles_(i, 2) * dipoles_(i, 2)));
    }
//    printDipoles();
}

//TODO get fit dipoles to match total dipole
void FieldMap::calcTotalDipole(const Frame &aa_frame, int num_atoms){
    if(num_atoms == 0) num_atoms = aa_frame.numAtomsTrack_;
    totalDipole_.zero();
    for(int i=0; i < num_atoms; i++){
        float charge = aa_frame.atoms_[i].charge;
        totalDipole_(0) += aa_frame.atoms_[i].coords[0] * charge;
        totalDipole_(1) += aa_frame.atoms_[i].coords[1] * charge;
        totalDipole_(2) += aa_frame.atoms_[i].coords[2] * charge;
    }
    totalDipole_(5) = float(sqrt(totalDipole_(0)*totalDipole_(0) +
            totalDipole_(1)*totalDipole_(1) +
            totalDipole_(2)*totalDipole_(2)));

    cout << "Total molecular dipole and sum of bead dipoles" << endl;
    totalDipole_.print(8, 4, constants::ENM2DEBYE);

    sumDipoles_.zero();
    for(int i=0; i < numDipoles_; i++){
        for(int j=0; j<6; j++){
            sumDipoles_(j) += dipoles_(i, j);
        }
    }
    sumDipoles_.print(8, 4, constants::ENM2DEBYE);
}

void FieldMap::calcSumDipole(const vector<int> nums){
    sumDipoles_.zero();
    for(const int &i : nums){
        sumDipoles_(0) += dipoles_(i, 0);
        sumDipoles_(1) += dipoles_(i, 1);
        sumDipoles_(2) += dipoles_(i, 2);
    }
    sumDipoles_(5) = float(sqrt(sumDipoles_(0)*sumDipoles_(0) +
            sumDipoles_(1)*sumDipoles_(1) +
            sumDipoles_(2)*sumDipoles_(2)));

    sumDipoles_.print(8, 4, constants::ENM2DEBYE);
}

void FieldMap::printDipoles(){
    cout << "  Dipx    Dipy    Dipz    Polt    Polp    Polm" << endl;
    dipoles_.print(8, 4, constants::ENM2DEBYE);
}

//TODO move this outside the class - it doesn't need to be here
float FieldMap::distSqr(const float *coords, const float x, const float y, const float z) {
    float tmpx = coords[0] - x;
    float tmpy = coords[1] - y;
    float tmpz = coords[2] - z;
    return tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
}

void polar(const float *cart, float *polar){

}

