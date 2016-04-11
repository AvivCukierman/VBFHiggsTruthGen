# --------------------------------------------- #
# Makefile for Recluster code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS = -O2 -Wall -Wextra -std=c++0x -g

.PHONY: clean debug all

all: setup Pileup

setup:
	mkdir -p lib

Pileup:  lib/Pileup.so lib/PileupAnalysis.so lib/Configuration.so
	$(CXX) lib/Pileup.so lib/PileupAnalysis.so lib/Configuration.so -o $@ \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-I./include -L./lib \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs` \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

lib/Pileup.so: src/Pileup.C lib/PileupAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags` \
	-I./include -L./lib \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

lib/PileupAnalysis.so : src/PileupAnalysis.cc include/PileupAnalysis.h 
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags` \
	-I./include \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs` 

lib/Configuration.so : src/Configuration.cc include/Configuration.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags` \
	-I./include \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags --libs`

clean:
	rm -rf Pileup
	rm -rf lib
	rm -f *~

install:
	install Pileup -t ${HOME}/local/bin
	install setupPileup.sh -t ${HOME}/local/bin
	install scripts/Pileup.sh -t ${HOME}/local/bin
