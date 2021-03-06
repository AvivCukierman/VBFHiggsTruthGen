//-*-c++-*-

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <iostream>
#include <sstream>
#include "boost/program_options.hpp"
#include "Definitions.h"

using namespace std;
namespace po = boost::program_options;

class Configuration{
private:
  po::variables_map vm;
  int getSeed(int seed);
public:
  Configuration(int argc, char* argv[]);
  void ConfigurePythiaSignal(Pythia8::Pythia* hs);
  void ConfigurePythiaPileup(Pythia8::Pythia* pu);
  void print();
  
  string outName;
  string outDir;
  bool autoName;
  string procName;

  int    pileup;
  int    nEvents;
  int    fDebug;
  int    proc;
  int    seed;

  float  bunchsize;
  float  minEta;
  float  maxEta;  
  float  pThatmin;
  float  pThatmax;
  float  scalar_mass;

};

#endif
