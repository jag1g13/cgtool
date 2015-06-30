//
// Created by james on 20/03/15.
//
#include "small_functions.h"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <cassert>

#include <sys/stat.h>
#include <fstream>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

using std::string;
using std::cout;
using std::endl;
using std::vector;

bool file_exists(const string name){
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

long file_size(const string filename){
    struct stat buffer;
    int rc = stat(filename.c_str(), &buffer);
    return rc == 0 ? buffer.st_size : -1;
}

double start_timer(){
    struct timespec now;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  now.tv_sec = mts.tv_sec;
  now.tv_nsec = mts.tv_nsec;
#else
    if(clock_gettime(CLOCK_REALTIME, &now) == -1) throw std::runtime_error("No clock");
#endif

    return now.tv_sec + 1e-9 * now.tv_nsec;
}

double end_timer(const double start){
    return start_timer() - start;
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

void split_text_output(const string &name, const double start){
    // If time has passed, how much?  Ignore small times
    double time = end_timer(start);
    cout << endl;
    if(time > 0.1){
        cout << "--------------------" << endl;
        cout << time << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

void polar(const double cart[3], double polar[3]){
    throw std::logic_error("Not implemented yet");
}

uint32_t u4_from_buffer(char b[]){
    return uint8_t(b[3]) <<  0 | uint8_t(b[2]) <<  8
           | uint8_t(b[1]) << 16 | uint8_t(b[0]) << 24;
}

int get_xtc_num_frames(const string &xtcname){
    // Based on http://comments.gmane.org/gmane.science.biology.gromacs.devel/3901
    std::ifstream xtc;
    xtc.open(xtcname.c_str(), std::ios::binary);
    if(!xtc.is_open()) return -1;

    // Get frame size
    char header[92];
    xtc.read(header, 92);
    xtc.close();
    uint32_t frame_size = u4_from_buffer(header+88);

    const long total_size = file_size(xtcname);
    return total_size / frame_size;
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
