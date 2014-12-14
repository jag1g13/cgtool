#!/bin/bash

cd build
cmake ..
make clean
make all
ctest -V
#make test
