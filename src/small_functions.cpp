//
// Created by james on 20/03/15.
//
#include "small_functions.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

#include <sys/stat.h>

using std::string;
using std::cout;
using std::endl;

bool file_exists(const std::string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

bool backup_old_file(const std::string name){
    if(!file_exists(name)) return true;

    string newName = "#" + name + "#";
    int i = 1;
    while(file_exists((newName+std::to_string(i)).c_str())) i++;
    int check = std::rename(name.c_str(), (newName+std::to_string(i)).c_str());
    if(check == 0) return true;
    throw std::runtime_error("File could not be backed up");
}

void split_text_output(const string name, const clock_t start, const int num_threads){
    clock_t now = std::clock();
    // If time has passed, how much?  Ignore small times
    if(double(now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << double(now - start) / (CLOCKS_PER_SEC * num_threads) << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

double dot(const double A[3], const double B[3]){
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

double abs(const double vec[3]){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double distSqr(const double c1[3], const double c2[3]){
    return  (c1[0]-c2[0])*(c1[0]-c2[0]) +
            (c1[1]-c2[1])*(c1[1]-c2[1]) +
            (c1[2]-c2[2])*(c1[2]-c2[2]);
}

void polar(const double cart[3], double polar[3]){
    throw std::logic_error("Not implemented yet");
}

