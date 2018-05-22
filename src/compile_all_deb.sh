#!/bin/bash

time g++ -c ANAwrite.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAread.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAutils.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAstatic.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAmd.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAndd.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c ANAPO.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include
time g++ -c main.cpp -g -O3 -std=c++14 -Wall -isystem -rdynamic -I/usr/local/include -I/usr/include -I../chemfiles/install/include

g++ *o  -L/usr/lib -L/usr/local/lib -L../chemfiles/install/lib -lchemfiles -lnetcdf -lCGAL -lboost_program_options -lmpfr -lgmp -lm -o ../ANA
