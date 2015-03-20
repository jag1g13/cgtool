//
// Created by james on 20/03/15.
//
#include "small_functions.h"

#include <stdexcept>
#include <iostream>
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
    if((float) (now - start) / CLOCKS_PER_SEC > 0.1){
        cout << "--------------------" << endl;
        cout << float(now - start) / (CLOCKS_PER_SEC * num_threads) << " seconds" << endl;
    }
    cout << "====================" << endl;
    cout << name << endl;
    cout << "--------------------" << endl;
}

