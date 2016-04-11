//-*-c++-*-

#ifndef  PileupAnalysis_H
#define  PileupAnalysis_H

#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

#include "Configuration.h"
#include "Definitions.h"
#include "TimingInfo.h"

using namespace std;
using namespace fastjet;

typedef vector<float> timingBranch;

double distance(double eta, double phi, double truthEta, double truthPhi);

class TimingDistribution{
 private:
  float _bunchsize;
  int _seed;
  mt19937 rng;  // mt19937 is a standard mersenne_twister_engine

  double _phi;
  double _psi;
  double _phi_nums[2];
  double _psi_nums[2];

  double _gauss_norm;
  double _square_norm;

  double probability(double zpos, double time, distribution dtype);
  int randomSeed();

 public:
  TimingDistribution(float bunchsize=0.075, int seed=-1, double phi=0.0, double psi=0.0);
  void phi(double phi);
  void psi(double psi);
  double psi(){return _psi;};
  double phi(){return _phi;};
  pair<double,double> get(distribution dtype=gaussian);
  double uniform(double min=0, double max=1);
};

class PileupAnalysis{
 private:
  int  ftest;
  bool  fDebug;
  string fOutName;
  
  TFile *tF;
  TTree *tT;
  unique_ptr<TimingDistribution> rnd;

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

  distribution _dtype;
  double psi;
  double phi;
  
  // Tree Vars ---------------------------------------
  int fTEventNumber;
  int fTNPV;
  float fzvtxspread;
  float ftvtxspread;
  
  timingBranch *zpu;
  timingBranch *j0clpt;
  timingBranch *j0clphi;
  timingBranch *j0cleta;
  timingBranch *j0clm;
  timingBranch *j0cltruth;
  timingBranch *j0clpu;
  timingBranch *j0clcharge;
  timingBranch *j0clpdgid;

  timingBranch *truejpt;
  timingBranch *truejphi;
  timingBranch *truejeta;

  bool randomZ;
  bool randomT;
  bool smear;
  bool displace;

  Pythia8::Pythia *_pythiaHS;
  Pythia8::Pythia *_pythiaPU;
  
  void DeclareBranches();
  void ResetBranches();
  void FillTruthTree(JetVector jets);
  bool Ignore(Pythia8::Particle &p);

  //Jet selection functions
  void selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets);
  
 public:
  PileupAnalysis (Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q);
  ~PileupAnalysis ();
  
  void AnalyzeEvent(int iEvt, int NPV);
  void Initialize(float minEta, float maxEta, distribution dtype=gaussian,int seed=123);
  
  //settings (call before initialization)
  void Debug(int debug){fDebug = debug;}
  void Bunchsize(float bunchsize_){bunchsize=(bunchsize_>0)?bunchsize_:0;}
  void PileupMode(smearMode PU);
  void SignalMode(smearMode HS);
  void Phi(double phi);
  void Psi(double psi);
  void SetOutName(string outname){fOutName = outname;}
  
};

#endif

