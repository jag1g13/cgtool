//
// Created by james on 22/05/15.
//

#ifndef CGTOOL_RDF_H
#define CGTOOL_RDF_H

#include <vector>

#include "residue.h"
#include "frame.h"
#include "histogram.h"
#include "array.h"

class RDF{
protected:
    std::vector<Residue> residues_;
    Histogram histogram_;
    Array rdf_;
    int resolution_ = 100;
    double cutoff_ = 2.;
    int grid_ = 200;
    int frames_ = 0;
    double density_ = 0.;

public:
    RDF(){};
    RDF(const std::vector<Residue> &residues);
    RDF(const std::vector<Residue> &residues, const double cutoff, const int resolution);

    void init(const std::vector<Residue> &residues,
              const double cutoff, const int resolution);

    void calculateRDF(const Frame &frame);
    void normalize();
};

#endif //CGTOOL_RDF_H
