//-*-c++-*-

#ifndef  Analysis_H
#define  Analysis_H

#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "TLorentzVector.h"

#include "Configuration.h"
#include "Definitions.h"
#include "ParticleInfo.h"

using namespace std;
using namespace fastjet;

double distance(double eta, double phi, double truthEta, double truthPhi);

class Analysis{
 private:
  int  ftest;
  bool  fDebug;
  string fOutName;
  string fOutDir;
  
  TFile *tF;
  TTree *tT;
  //unique_ptr<TimingDistribution> rnd;

  //these need to be ptrs because constructor must be called
  unique_ptr<JetDefinition> jetDef;
  unique_ptr<AreaDefinition> active_area;
  unique_ptr<GridMedianBackgroundEstimator> bge;
  unique_ptr<FunctionOfPseudoJet<double> > rescaling;
  unique_ptr<Selector> select_fwd;
  
  float bunchsize;
  float _minEta;
  float _maxEta;
  double _R;
  double _minP;
  
  // Tree Vars ---------------------------------------
  int fTEventNumber;
  int fTNPV;
  int Ntruthphotons;
  
  branch *gammapt;
  branch *gammaphi;
  branch *gammaeta;
  branch *gammam;

  branch *j0pt;
  branch *j0phi;
  branch *j0eta;
  branch *j0m;
  branch *j0id;

  Pythia8::Pythia *_pythiaHS;
  Pythia8::Pythia *_pythiaPU;
  
  void DeclareBranches();
  void ResetBranches();
  //void FillTruthTree(JetVector jets);
  bool Ignore(Pythia8::Particle &p);

  //Jet selection functions
  void selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets, vector<int> &jetIDs);
  
 public:
  Analysis (Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q);
  ~Analysis ();
  
  void AnalyzeEvent(int iEvt, int NPV);
  void Initialize(float minEta, float maxEta);
};

#endif

