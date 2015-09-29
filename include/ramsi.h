//
// Created by james on 10/08/15.
//

#ifndef RAMSI_RAMSI_H
#define RAMSI_RAMSI_H

#include "common.h"

#include "membrane.h"

class Ramsi : public Common{
protected:
    Membrane *membrane_ = nullptr;
    std::vector<double> thickness_;

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
    Ramsi(){};
    virtual ~Ramsi();
};

#endif //RAMSI_RAMSI_H
