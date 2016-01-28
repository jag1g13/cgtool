//
// Created by james on 13/08/15.
//

#ifndef CGTOOL_CGTOOL_H
#define CGTOOL_CGTOOL_H

#include "common.h"

#include "bondset.h"
#include "rdf.h"

class Cgtool : public Common{
protected:
    // Output file formats
    FileFormat outProgram_;
    FieldFormat outField_;
    /** \brief Types of the three bonded potentials */
    PotentialType potentialTypes_[3];

    BondSet  *bondSet_ = nullptr;
    RDF      *rdf_ = nullptr;

    TrjOutput *trjOutput_ = nullptr;

    // Protected functions
    /** \brief Read config file and determine which functions should be performed */
    void readConfig();

    /** \brief Construct objects which are required to perform requested functions */
    void setupObjects();

    /** \brief Function executed within the main loop - performs most significant work*/
    void mainLoop();

    /** \brief Perform final calculations and end program */
    void postProcess();

public:
    Cgtool(){};
    virtual ~Cgtool();
};


#endif //CGTOOL_CGTOOL_H
