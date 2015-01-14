#include "field_map.h"

//#include <algorithm>
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
//    cout << "Field Monopole" << endl;
    fieldMonopole_.init(a, b, c, false);
//    cout << "Field Dipole" << endl;
    fieldDipole_.init(a, b, c, false);
//    cout << "Grid bounds" << endl;
    gridBounds_.init(3, 2, 1, false);
//    cout << "Coords" << endl;
    gridCoords_.init(3, max(a, max(b, c)), 1, false);
//    cout << "Dipoles" << endl;
    numDipoles_ = ndipoles;
    //TODO init this later
    dipoles_.init(ndipoles, 6, 1, false);
//    cout << "Grid contracted" << endl;
    //TODO try out Boost arrays, are they faster/better?
    gridContracted_.init(a*b*c, 4, 1, true);
    totalDipole_.init(6, 1, 1, false);
    sumDipoles_.init(6, 1, 1, false);
}

void FieldMap::setupGrid(const Frame *frame){
    /* create min and max initial values */
    gridBounds_(0, 0) = 1e6; gridBounds_(0, 1) = -1e6;
    gridBounds_(1, 0) = 1e6; gridBounds_(1, 1) = -1e6;
    gridBounds_(2, 0) = 1e6; gridBounds_(2, 1) = -1e6;
//    for(auto atom : frame->atoms_){
    for(int ii : frame->residues_[0].atoms){
        /* find bounding box of molecule */
        const Atom *atom = &(frame->atoms_[ii]);
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
    //cout << "Grid centre at: " << gridCentre_[0] << "," << gridCentre_[1] << "," << gridCentre_[2] << endl;

    gridCoords_.linspace(0, gridDims_[0], gridBounds_(0, 0), gridBounds_(0, 1));
    gridCoords_.linspace(1, gridDims_[1], gridBounds_(1, 0), gridBounds_(1, 1));
    gridCoords_.linspace(2, gridDims_[2], gridBounds_(2, 0), gridBounds_(2, 1));
//    for(int i=0; i<3; i++){
//        /* for x, y, z do linspace of grid coordinates */
//        gridCoords_.linspace(i, gridBounds_(i, 0), gridBounds_(i, 1));
//    }
//    cout << "Grid bounds" << endl;
//    cout << gridBounds_(0, 0) << " " << gridBounds_(0, 1) << endl;
//    cout << gridBounds_(1, 0) << " " << gridBounds_(1, 1) << endl;
//    cout << gridBounds_(2, 0) << " " << gridBounds_(2, 1) << endl;
//    cout << "Grid coords" << endl;
//    for(int i=0; i<gridDims_[0]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(0, i) << " ";
//    }
//    cout << endl;
//    for(int i=0; i<gridDims_[1]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(1, i) << " ";
//    }
//    cout << endl;
//    for(int i=0; i<gridDims_[2]; i++){
//        if(i%10 == 0) cout << endl;
//        cout << gridCoords_(2, i) << " ";
//    }
//    cout << endl;
//    cout << "Grid setup done" << endl;
}

void FieldMap::setupGridContracted(const Frame *frame){
    /*RADMIN=50.0D0                                                     CHE01860        still inside x,y,z, loops
    187       DO 100 I=1,NATOMS                                                 CHE01870        for each atom
    188       VRAD = RADII(I)                                                   CHE01880
    189       DIST = (P1 - C(1,I))**2 + (P2 - C(2,I))**2 + (P3 - C(3,I))**2     CHE01890        calculate distance of grid point to atom
    190       DIST = DSQRT(DIST)                                                CHE01900
    191       IF (DIST .LT. VRAD) GOTO 210                                      CHE01910        if grid point is inside atom - skip
    192       IF (DIST .LT. RADMIN) RADMIN = DIST                               CHE01920        keep track of closest grid point
    193   100 CONTINUE                                                          CHE01930        do next atom
        194       IF (RADMIN .GT. RMAX) GOTO 210                                    CHE01940        if grid point is too far away from all atoms - skip
    195 C                                                                       CHE01950
    196 C        STORE POINTS (IN ATOMIC UNITS)                                 CHE01960
    197 C                                                                       CHE01970
    198       IPOINT = IPOINT + 1                                               CHE01980
    199       P(1,IPOINT) = P1                                                  CHE01990        store grid point if accepted
        200       P(2,IPOINT) = P2                                                  CHE02000
    201       P(3,IPOINT) = P3                                                  CHE02010
    202       IF (IPO(2) .EQ. 1)                                                CHE02020
    203      $ WRITE(IOUT,*) 'POINT ',IPOINT,' X,Y,Z ',P1,P2,P3                 CHE02030
    204   210 CONTINUE                                                          CHE02040
    205   200 CONTINUE
    */
    float radmin = 500.f;
    float dist = 0.f;
    float rmax = border_;    // use an rmax equal to border_ around molecule
    float vrad = 0.1f;       // reject if within 1A of an atom (inside atomic radius)
//    float x, y, z;
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
                //for(Atom atom : frame->atoms_){
//                for(int ii : frame->residues_[0].atoms){
//                cout << i << j << k << endl;
                for(int ii=0; ii < frame->residues_[0].atoms.size(); ii++){
                    dist = sqrt(distSqr(frame->atoms_[ii].coords, coords[0], coords[1], coords[2]));
//                    cout << ii << endl;
//                    if(dist < border_) cout << dist << endl;
                    radmin = min(radmin, dist);
                    if(dist < vrad){
                        accepted = false;
//                        cout << "too close" << endl;
                        ++close_count;
                        break;
                    }
//                    if(dist < radmin) radmin = dist;
                }
                if(accepted && radmin > rmax){
//                    cout << "too far" << endl;
//                    cout << radmin << endl;
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
//    cout << "a " << fieldMonopoleContracted_.size() << "\tb " << fieldDipoleContracted_.size() << endl;
//    cout << "Acc\tFar\tClose" << endl;
//    cout << accepted_count << "\t" << far_count << "\t" << close_count << endl;
}

void FieldMap::calcFieldMonopoles(const Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // inveps = 8.9875517873681e9
    float x, y, z;
    #pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        x = gridCoords_(0, i);
        for(int j=0; j < gridDims_[1]; j++){
            y = gridCoords_(1, j);
            for(int k=0; k < gridDims_[2]; k++){
                z = gridCoords_(2, k);
                fieldMonopole_(i, j, k) = 0.;
                for(Atom atom : frame->atoms_){
                    fieldMonopole_(i, j, k) += atom.charge /
                            distSqr(atom.coords, x, y, z);
                }
            }
        }
    }
}

void FieldMap::calcFieldMonopolesContracted(const Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // inveps = 8.9875517873681e9
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        fieldMonopoleContracted_[i] = 0.f;
        for(Atom atom : frame->atoms_) {
            fieldMonopoleContracted_[i] += atom.charge /
                    distSqr(atom.coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));
        }
    }
}

inline float dot(const float *A, const float* B){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

inline float abs(const float* vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void FieldMap::printFields(){
    cout << "ELECTRIC FIELDS" << endl;
    for(int i=0; i<numGridPoints_; i++){
        cout << fieldMonopoleContracted_[i] << "\t" << fieldDipoleContracted_[i] << endl;
    }
}

//TODO dipole field - to check that they're right - pretty much done, needs testing
void FieldMap::calcFieldDipolesContracted(const Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // float inveps = 8.9875517873681e9
    float vec_a[3], vec_b[3];
    float abs_a;
//    cout << numGridPoints_ << endl;
//    cout << frame->numAtomsTrack_ << endl;
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        fieldDipoleContracted_[i] = 0.f;
//        cout << "i=" << i << endl;
        // for charge on the cg bead
//        fieldDipoleContracted_[i] += dipoles_(j, 5) / (abs_a*abs_a);
        for(int j=0; j < frame->numAtomsTrack_; j++){
//            cout << "j=" << j << endl;
            // each 'atom' in frame, actually cg bead
            for(int k=0; k<3; k++) {
                // vector from point dipole to grid point
                vec_a[k] = gridContracted_(i, k) - frame->atoms_[j].coords[k];
                // copy of dipole vector
                vec_b[k] = dipoles_(j, k);
            }
            abs_a = abs(vec_a);
//            cout << dipoles_(j, 5) << "\t" << abs_a << endl;
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
        StatsBox sb = vector_stats(&fieldMonopoleContracted_, &fieldDipoleContracted_);
        cout << "\tRMS: " << sb.rms << "\tRRMS: " << sb.rrms << endl;
    }
}

void FieldMap::calcFieldDipoles(const Frame *frame) {
    throw std::runtime_error("Not implemented");
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    //float inveps = 8.9875517873681e9;
    float x, y, z;
    float vec_a[3], vec_b[3];
    float cos_dip_angle;
    #pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        x = gridCoords_(i, 0);
        for(int j=0; j < gridDims_[1]; j++){
            y = gridCoords_(j, 1);
            for(int k=0; k < gridDims_[2]; k++){
                z = gridCoords_(k, 2);
                fieldDipole_(i, j, k) = 0.f;
                int ii = 0;
                for(Atom atom : frame->atoms_){
//                    vec_a[0] =
//                    cos_dip_angle = dot()
                    fieldDipole_(i, j, k) += dipoles_(ii, 5) * cos_dip_angle /
                            distSqr(atom.coords, x, y, z);
                    ii++;
                    //self.grid_dipole[i][j][k] += dipole[6] * cos(dip_angle)
                }
            }
        }
    }
}

/*
def calc_dipoles(cg_frame, frame, out_file, outfile_sum, export=True,
                 cg_internal_bonds=cg_internal_bonds,
                 sugar_atom_nums=sugar_atom_nums, adjacent=adjacent):
    """
    dipole of a charged fragment is dependent on where you measure it from
    so now includes flag to make beads neutral
    should modify this to include OW dipoles
    """
    charge_redist = True   # redistribute charge, make all beads neutral
    dipoles = []
    frame_dipoles = np.zeros((len(cg_sites), 3))
    for i, site in enumerate(cg_sites):
        num_atoms = float(len(cg_internal_map[site]))
        dipole = np.zeros(3)
        dipole_sum = np.zeros(3)
        cg_atom = cg_frame.atoms[cg_atom_nums[site]]
        for atom_name in cg_internal_map[site]:
            atom = frame.atoms[sugar_atom_nums[atom_name]]
            charge = atom.charge
            if charge_redist:
                charge -= cg_atom.charge / num_atoms
            dipole += atom.loc * charge    # first calc from origin
        if not charge_redist:
            dipole -= cg_atom.loc * cg_atom.charge  # then recentre it
        # sum the dipoles, check that they equal the total molecular dipole
        dipole_sum += dipole
        norm, bisec = cg_frame.norm_bisec(cg_atom_nums[adjacent[site][0]], i,
                                          cg_atom_nums[adjacent[site][1]])
        frame_dipoles[i] += polar_coords(dipole, norm, bisec)
        if export:
            np.savetxt(out_file, frame_dipoles, delimiter=",")
        dipoles.append(frame_dipoles)
    if export:
        np.savetxt(outfile_sum, dipole_sum, delimiter=",")
    return dipoles
*/

/**
* This uses a hack to allow direct calculation: charges within a bead are
* rescaled to make each bead neutral.  This means that dipoles are no longer
* dependent on the frame of reference.
* The electric field from these dipoles should be compared against the field
* from atomic point charges to determine validity.  They may need to be rescaled.
* I don't see a better way to do this.
*/
void FieldMap::calcDipolesDirect(const CGMap *cgmap, const Frame *cg_frame, Frame *aa_frame){
//    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    dipoles_.zero();
//    for(BeadMap &bead : cgmap->mapping_){
    for(int i=0; i<cgmap->num_beads; i++){
//        cout << "doing bead #" << i << endl;
        const BeadMap &bead_type = cgmap->mapping_[i];
        const Atom &cg_atom = cg_frame->atoms_[i];
        // for each bead in the CG frame
        for(const int &j : bead_type.atom_nums){
//            cout << "atom #" << j << endl;
            // for each atom inside the bead
            // rescale charges so bead charge is zero
            float charge = aa_frame->atoms_[j].charge;
//            charge -= cg_atom.charge / bead_type.num_atoms;
            dipoles_(i, 0) += aa_frame->atoms_[j].coords[0] * charge;
            dipoles_(i, 1) += aa_frame->atoms_[j].coords[1] * charge;
            dipoles_(i, 2) += aa_frame->atoms_[j].coords[2] * charge;
        }
        float charge = cg_frame->atoms_[i].charge;
//        cout << charge << endl;
        dipoles_(i, 0) -= cg_frame->atoms_[i].coords[0] * charge;
        dipoles_(i, 1) -= cg_frame->atoms_[i].coords[1] * charge;
        dipoles_(i, 2) -= cg_frame->atoms_[i].coords[2] * charge;
        // calculate magnitude here
        dipoles_(i, 5) = float(sqrt(dipoles_(i, 0)*dipoles_(i, 0) +
                              dipoles_(i, 1)*dipoles_(i, 1) +
                              dipoles_(i, 2)*dipoles_(i, 2)));
        //TODO Dipoles are about 10x smaller than they should be - why?
//        dipoles_(i, 0) *= 10.f;
//        dipoles_(i, 1) *= 10.f;
//        dipoles_(i, 2) *= 10.f;
//        dipoles_(i, 5) *= 10.f;
    }
    printDipoles();
}

void FieldMap::calcDipolesFit(const CGMap *cgmap, const Frame *cg_frame, const Frame *aa_frame){
    dipoles_.zero();
    // keep track of which dipoles we know
    vector<int> dipoles_calculated;
    // calculate dipoles on all uncharged beads
    for(int i = 0; i < cgmap->num_beads; i++){
        const BeadMap &bead_type = cgmap->mapping_[i];
        // skip beads that are charged
        if(bead_type.charge > 0.01f || bead_type.charge < -0.01f) continue;

        const Atom &cg_atom = cg_frame->atoms_[i];
        // for each bead in the CG frame
        for(const int &j : bead_type.atom_nums){
            // for each atom inside the bead
            float charge = aa_frame->atoms_[j].charge;
            dipoles_(i, 0) += aa_frame->atoms_[j].coords[0] * charge;
            dipoles_(i, 1) += aa_frame->atoms_[j].coords[1] * charge;
            dipoles_(i, 2) += aa_frame->atoms_[j].coords[2] * charge;
        }
        dipoles_calculated.push_back(i);
    }
    int num_remaining = cgmap->num_beads - int(dipoles_calculated.size());

    // calculate residual molecular dipole
    calcTotalDipole(aa_frame);
    calcSumDipole(dipoles_calculated);
    // keep residual dipole in place
    totalDipole_ -= sumDipoles_;

    // divide residual between remaining beads
    for(int i = 0; i < cgmap->num_beads; i++){
        const BeadMap &bead_type = cgmap->mapping_[i];
        // skip beads that are NOT charged
        if(bead_type.charge < 0.01f && bead_type.charge > -0.01f) continue;

        for(int j = 0; j < 3; j++){
            dipoles_(i, j) = totalDipole_(j) / num_remaining;
        }
    }

    // calculate magnitudes
    for(int i = 0; i < cgmap->num_beads; i++){
        dipoles_(i, 5) = float(sqrt(dipoles_(i, 0) * dipoles_(i, 0) +
                dipoles_(i, 1) * dipoles_(i, 1) +
                dipoles_(i, 2) * dipoles_(i, 2)));
    }
    printDipoles();
}

//TODO get fit dipoles to match total dipole
void FieldMap::calcTotalDipole(const Frame *aa_frame, int num_atoms){
    if(num_atoms == 0) num_atoms = aa_frame->numAtomsTrack_;
    totalDipole_.zero();
    for(int i=0; i < num_atoms; i++){
    // only add dipoles we're confident of
//    for(int i=0; i < 9; i++){
        float charge = aa_frame->atoms_[i].charge;
        totalDipole_(0) += aa_frame->atoms_[i].coords[0] * charge;
        totalDipole_(1) += aa_frame->atoms_[i].coords[1] * charge;
        totalDipole_(2) += aa_frame->atoms_[i].coords[2] * charge;
    }
    totalDipole_(5) = float(sqrt(totalDipole_(0)*totalDipole_(0) +
            totalDipole_(1)*totalDipole_(1) +
            totalDipole_(2)*totalDipole_(2)));

    cout << "Total molecular dipole" << endl;
    cout << "Sum of bead dipoles" << endl;
    totalDipole_.print(8, 4, constants::ENM2DEBYE);
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
    cout << "Dipx\tDipy\tDipz\tPolt\tPolp\tPolm" << endl;
    dipoles_.print(8, 4, constants::ENM2DEBYE);
}

//TODO move this outside the class - it doesn't need to be here
/**
* Originally used cmath pow - this version is much faster.
*/
float FieldMap::distSqr(const float *coords, const float x, const float y, const float z) {
//    return (coords[0] - x)*(coords[0] - x) +
//            (coords[1] - y)*(coords[1] - y) +
//            (coords[2] - z)*(coords[2] - z);
    float tmpx = coords[0] - x;
    float tmpy = coords[1] - y;
    float tmpz = coords[2] - z;
    return tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
}

void polar(const float *cart, float *polar){

}

