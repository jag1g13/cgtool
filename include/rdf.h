//
// Created by james on 22/05/15.
//

#ifndef CGTOOL_RDF_H
#define CGTOOL_RDF_H

#include <vector>

#include "residue.h"
#include "frame.h"

class RDF{
protected:
    std::vector<Residue> residues_;

public:
    RDF(const std::vector<Residue> &residues);

    void calculateRDF(const Frame &frame);
};

#endif //CGTOOL_RDF_H
