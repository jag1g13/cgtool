#!/bin/bash

rm -r build/*
cd build
cmake ..
make
#ctest -V
make test
