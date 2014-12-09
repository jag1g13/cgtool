#!/bin/bash

cd build
cmake ..
make clean
make
ctest -V
#make test
