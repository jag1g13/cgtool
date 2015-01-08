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
    dipoles_.init(ndipoles, 6, 1, false);
//    cout << "Grid contracted" << endl;
    gridContracted_.init(a*b*c, 4, 1, true);
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
    fieldMonopoleContracted_.reserve(numGridPoints_);
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
        for(Atom atom : frame->atoms_) {
            fieldMonopoleContracted_[i] += atom.charge /
                    distSqr(atom.coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));
        }
    }
}

//TODO dipole field - to check that they're right
void FieldMap::calcFieldDipolesContracted(const Frame *frame){
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    // inveps = 8.9875517873681e9
#pragma omp parallel for
    for(int i=0; i < numGridPoints_; i++) {
        for(int j=0; j < frame->numAtomsTrack_; j++){
//            fieldMonopoleContracted_[i] += atom.charge /
//                    distSqr(atom.coords, gridContracted_(i, 0), gridContracted_(i, 1), gridContracted_(i, 2));

        }
    }
}

void FieldMap::calcFieldDipoles(const Frame *frame) {
    float inveps = 1. / (4 * M_PI * 8.854187817e-12);
    //float inveps = 8.9875517873681e9;
    float x, y, z;
    float cos_dip_angle = 0;
    #pragma omp parallel for
    for(int i=0; i < gridDims_[0]; i++){
        x = gridCoords_(i, 0);
        for(int j=0; j < gridDims_[1]; j++){
            y = gridCoords_(j, 1);
            for(int k=0; k < gridDims_[2]; k++){
                z = gridCoords_(k, 2);
                fieldDipole_(i, j, k) = 0.;
                int ii = 0;
                //TODO dot product: cos(theta) = A.B / |A||B|
                for(Atom atom : frame->atoms_){
                    
                    cos_dip_angle = dot()
                    fieldDipole_(i, j, k) += dipoles_(ii, 5) * cos_dip_angle /
                            distSqr(atom.coords, x, y, z);
                    ii++;
                    //self.grid_dipole[i][j][k] += dipole[6] * cos(dip_angle)
                }
            }
        }
    }
}

inline float dot(const float *A, const float* B){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
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
        dipoles_(i, 5) = sqrt(dipoles_(i, 0)*dipoles_(i, 0) +
                              dipoles_(i, 1)*dipoles_(i, 1) +
                              dipoles_(i, 2)*dipoles_(i, 2));
    }
}

void FieldMap::printDipoles(){
    cout << "Dipx\tDipy\tDipz\tPolt\tPolp\tPolm" << endl;
    dipoles_.print();
}

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

