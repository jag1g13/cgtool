//
// Created by james on 26/05/15.
//

#ifndef CGTOOL_HISTOGRAM_H
#define CGTOOL_HISTOGRAM_H

class Histogram{
protected:
    int size_;
    int* array_ = nullptr;
    bool fast_ = false;
    bool allocated_ = false;
    double min_ = 0.;
    double max_ = 0.;
    double step_ = 0.;

public:
    // ##############################################################################
    // Con/Destructors
    // ##############################################################################
    Histogram(const int size, const bool fast=false);
    Histogram(){};
    ~Histogram();

    // ##############################################################################
    // Setup / Tear down
    // ##############################################################################
    void init(const int size, const bool fast=false);
    void zero(const int set=0);
    void free();

    // ##############################################################################
    // Add values
    // ##############################################################################
    void increment(int loc);
    void decrement(int loc);
    int at (int loc) const;

    void scale(const double mult);

    // ##############################################################################
    // Printing
    // ##############################################################################
    void print(const int width=8) const;
    void printGraph(const int scale=10) const;
};

#endif //CGTOOL_HISTOGRAM_H
