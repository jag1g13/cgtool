//
// Created by james on 22/05/15.
//

#ifndef CGTOOL_RDF_H
#define CGTOOL_RDF_H

#include <vector>

#include "residue.h"
#include "frame.h"
#include "histogram.h"
#include "light_array.h"

class RDF{
protected:
    const std::vector<Residue> &residues_;
    double cutoff_ = 2.;
    int resolution_ = 100;

    int grid_ = 200;

    int frames_ = 0;
    double density_ = 0.;

    Histogram histogram_;
    LightArray<double> rdf_;

public:
    RDF(const std::vector<Residue> &residues, const double cutoff, const int resolution) :
        residues_(residues), cutoff_(cutoff), resolution_(resolution){
        grid_ = static_cast<int>(cutoff_ * resolution_);
        histogram_.init(grid_);
        rdf_.alloc(grid_);
    };

    void calculateRDF(const Frame &frame);
    void normalize();
};

#endif //CGTOOL_RDF_H
