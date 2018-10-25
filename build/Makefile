CXX = /usr/bin/g++
CXXFLAGS = -g -O3 -std=c++14 -Wall
#LinOpt = -L/usr/lib -L/usr/local/lib -L../chemfiles/build/install/lib -lchemfiles -lnetcdf -lCGAL -lboost_program_options -lmpfr -lgmp -lm
LinOpt = ../chemfiles/build/install/lib/libchemfiles.a -lnetcdf -lCGAL -lboost_program_options -lmpfr -lgmp -lm
IncOpt = -I../chemfiles/build/install/include -I../include

SrcDir = ../src
IncDir = ../include/ANA

DepDir = dep
ObjDir = obj
BinDir = bin

Sources := $(wildcard ../src/*cpp)
Objects := $(Sources:$(SrcDir)/%.cpp=%.o)
Dependencies := $(Sources:$(SrcDir)/%.cpp=%.d)

.PHONY: build_msg all clean clean-fail
.DEFAULT_GOAL := all

all: build_msg ANA

ANA: $(Objects)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LinOpt)
	mkdir bin dep obj
	mv *.o obj
	mv *.d dep
	mv ANA bin

%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) $(IncOpt) -MMD -MP -c $< -o $@

clean-fail:
	rm *.o *.d ANA
	rm -r bin dep obj

clean:
	rm -r bin dep obj

build_msg:
	@echo "\n Building ANA \n"
