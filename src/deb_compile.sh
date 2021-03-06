#!/bin/bash

gcc="g++"

flags="-g -O3 -std=c++14 -Wall -isystem -rdynamic"

ch_include=../chemfiles/install/include
ch_lib=../chemfiles/install/lib

compi="-I/usr/local/include -I/usr/include -I${ch_include}"
linker="-L${ch_lib} -L/usr/lib -L/usr/local/lib -lchemfiles -lnetcdf -lCGAL -lboost_program_options -lmpfr -lgmp -lm"


echo -e "\n $gcc main.cpp $flags $compi $linker -o ../ANA \n"

$gcc main.cpp $flags $compi $linker -o ../ANA

exit 0
