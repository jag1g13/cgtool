//
// Created by james on 20/03/15.
//
#include "small_functions.h"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cassert>

#include <sys/stat.h>

using std::string;
using std::cout;
using std::endl;
using std::vector;

bool file_exists(const string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

long file_size(const string filename)
{
    struct stat buffer;
    int rc = stat(filename.c_str(), &buffer);
    return rc == 0 ? buffer.st_size : -1;
}

float time_since(const clock_t &since, const int num_threads){
    return float(std::clock() - since) / (CLOCKS_PER_SEC * num_threads);
}

bool backup_old_file(const string name){
    if(!file_exists(name)) return true;

    string newName = "#" + name + "#";
    int i = 1;
    while(file_exists((newName+std::to_string(i)).c_str())) i++;
    int check = std::rename(name.c_str(), (newName+std::to_string(i)).c_str());
    if(check == 0) return true;
    throw std::runtime_error("File could not be backed up");
}

void split_text_output(const string name, const clock_t start, const int num_threads){
    // If time has passed, how much?  Ignore small times
    float time = time_since(start, num_threads);
    cout << endl;
    if(time > 0.1f){
        cout << "--------------------" << endl;
        cout << time << " seconds" << endl;
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
    return (c1[0]-c2[0])*(c1[0]-c2[0]) +
           (c1[1]-c2[1])*(c1[1]-c2[1]) +
           (c1[2]-c2[2])*(c1[2]-c2[2]);
}

double distSqrPlane(const double c1[3], const double c2[3]){
    return (c1[0]-c2[0])*(c1[0]-c2[0]) +
           (c1[1]-c2[1])*(c1[1]-c2[1]);
}

void polar(const double cart[3], double polar[3]){
    throw std::logic_error("Not implemented yet");
}

StatsBox vector_stats(const vector<double> &a, const vector<double> &b){
    assert(a.size() == b.size());
    const int N = a.size();
    assert(a.size() != 0);
    assert(b.size() != 0);
    StatsBox result;
    result.min_a = a[0]; result.max_a = a[0];

    double sumsqr = 0.;
    for(int i=0; i<N; i++){
        sumsqr += (a[i] - b[i]) * (a[i] - b[i]);
        result.mean_a += a[i];
        result.mean_b += b[i];
        result.min_a = fmin(result.min_a, a[i]);
        result.max_a = fmax(result.max_a, a[i]);
    }

    result.mean_a /= N;
    result.mean_b /= N;
    result.diff_means = fabs(result.mean_a - result.mean_b);
    result.rmsd = sqrt(sumsqr / N);
    result.nrmsd = result.rmsd / (result.max_a - result.min_a);

    return result;
}
