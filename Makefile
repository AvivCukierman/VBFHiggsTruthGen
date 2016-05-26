# --------------------------------------------- #
# Makefile for Recluster code                        #
# Pascal Nef, March 6th 2014                    #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS = -O2 -Wall -Wextra -std=c++0x -g

.PHONY: clean debug all

all: setup Run 

setup:
	mkdir -p lib

Run:  lib/Run.so lib/Analysis.so lib/Configuration.so
	$(CXX) lib/Run.so lib/Analysis.so lib/Configuration.so -o $@ \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs` -lEG -lEGPythia8 \
	-I./include -L./lib \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs` \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -llhapdfdummy \
	-L$(BOOSTLIBLOCATION) -lboost_program_options 

lib/Run.so: src/Run.C lib/Analysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags` \
	-I./include -L./lib \
	-I$(PYTHIA8LOCATION)/include \
	-I $(BOOSTINCDIR) \
	`root-config --cflags` 

lib/Analysis.so : src/Analysis.cc include/Analysis.h 
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
	rm -rf Run 
	rm -rf lib
	rm -f *~

install:
	install Run -t ${HOME}/local/bin
	install setupPileup.sh -t ${HOME}/local/bin
	install scripts/Pileup.sh -t ${HOME}/local/bin
