#include "Analysis.h"

using namespace std;

double sgn(double val){
  if(val < 0)
    return -1;
  else if (val > 0)
    return 1;
  else
    return 0;
}

double distance(double eta, double phi,double truthEta, double truthPhi){
  double de=eta-truthEta;
  double dp=phi-truthPhi;
  return sqrt(pow(de,2)+pow(dp,2));
}

// Constructor 
PileupAnalysis::PileupAnalysis(Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q){

  if((pythiaHS != NULL) and (pythiaPU != NULL)){
    _pythiaHS=pythiaHS;
    _pythiaPU=pythiaPU;
  }
  else{
    cerr << "Invalid Pythia pointer passed to PileupAnalysis" << endl;
    exit(1);
  }

  fDebug=q.fDebug;
  if(fDebug) 
    cout << "PileupAnalysis::PileupAnalysis Start " << endl;

  ftest = 0;
  fOutName = q.outName;
  
  Bunchsize(q.bunchsize);
  PileupMode(q.PUmode);
  SignalMode(q.HSmode);
  Psi(q.psi);
  Phi(q.phi);

  if(fDebug) 
    cout << "PileupAnalysis::PileupAnalysis End " << endl;
}

void PileupAnalysis::PileupMode(smearMode PU){
  switch(PU){
  case Off:
    randomT=false;
    randomZ=false;
    break;
  case Z:
    randomT=false;
    randomZ=true;
    break;
  case T:
    randomT=true;
    randomZ=false;
    break;
  case ZT:
    randomT=true;
    randomZ=true;
  }
}

void PileupAnalysis::SignalMode(smearMode HS){
  switch(HS){
  case Off:
    smear=false;
    displace=false;
    break;
  case Z:
    smear=false;
    displace=true;
    break;
  case T:
    smear=true;
    displace=false;
    break;
  case ZT:
    smear=true;
    displace=true;
  }
}

void PileupAnalysis::Phi(double phi){
  this->phi= phi;
}

void PileupAnalysis::Psi(double psi){
  this->psi= psi;
}

// Destructor 
PileupAnalysis::~PileupAnalysis(){
  if(tT != NULL){
    tT->Write();
    tF->Close();
    delete tF;
  }

  if(zpu != NULL){
    delete zpu;
    delete j0clpt;
    delete j0clphi;
    delete j0cleta;
    delete j0clm;
    delete j0cltruth;
    delete j0clpu;
    delete j0clcharge;
    delete j0clpdgid;
    delete truejpt;
    delete truejphi;
    delete truejeta;
  }
}

// Begin method
void PileupAnalysis::Initialize(float minEta, float maxEta, distribution dtype, int seed){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("tree", "Event Tree for Timing");
   rnd.reset(new TimingDistribution(bunchsize,seed,phi,psi));
   _dtype=dtype;

   _minEta= (minEta > 0) ? minEta : 0;
   if(maxEta <= minEta){
     cerr << "Invalid Eta Limits " << minEta << " -> " << maxEta << "Passed to PileupAnalysis::Initialize" << endl;
     exit(20);
   }
   _maxEta= maxEta;

   const double R=0.4;
   const double grid_spacing(0.6);
   
   //suppress fastjet banner
   fastjet::ClusterSequence::set_fastjet_banner_stream(NULL);
   
   jetDef.reset(new JetDefinition(fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best));
   active_area.reset(new AreaDefinition(fastjet::active_area));
   bge.reset(new GridMedianBackgroundEstimator(_maxEta, grid_spacing));
   rescaling.reset(new BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05)); //function for rapidity rescaling of rho
   bge->set_rescaling_class(rescaling.get());
   select_fwd.reset(new Selector(SelectorAbsRapRange(_minEta,_maxEta)));
 
   DeclareBranches();
   
   zpu = new timingBranch();  

   j0clpt = new timingBranch();  
   j0clphi = new timingBranch();  
   j0cleta = new timingBranch();  
   j0clm = new timingBranch();  
   j0cltruth = new timingBranch();
   j0clpu = new timingBranch();
   j0clcharge = new timingBranch();
   j0clpdgid = new timingBranch();

   truejpt = new timingBranch();  
   truejphi = new timingBranch();  
   truejeta = new timingBranch();  

   ResetBranches();
   
   return;
}

// Analyze
void PileupAnalysis::AnalyzeEvent(int ievt, int NPV){

  if(fDebug) 
    cout << "PileupAnalysis::AnalyzeEvent Begin " << endl;
  
  // -------------------------
  if (!_pythiaHS->next()) return;
  if(fDebug) 
    cout << "PileupAnalysis::AnalyzeEvent Event Number " << ievt << endl;
  
  // reset branches 
  ResetBranches();
  
  // new event-----------------------
  fTEventNumber = ievt;
  JetVector particlesForJets;
  JetVector particlesForJets_np;
  
  //Pileup Loop
  
  fTNPV = NPV;
  std::pair<double,double> randomVariates;  
  randomVariates=rnd->get(_dtype);
  fzvtxspread = randomVariates.first;
  ftvtxspread = randomVariates.second;
  
  double zhs=0.0;
  double ths=0.0;
  if(displace) //randomly distribute "ideally measured" hard-scatter vertex
    zhs=fzvtxspread;
  if(smear)
    ths=ftvtxspread;

  static const double z0 = 3.5; // FIX THIS! DEFINE Z0 ONLY IN ONE PLACE! 

  //Loop over Pileup Events
  for (int iPU = 0; iPU < NPV; ++iPU) {
    
    //determine random vertex position in z-t space
    randomVariates=rnd->get(_dtype);
    double zvtx = 0;
    double tvtx = 0;
    if(randomZ)
      zvtx = randomVariates.first;
    if(randomT)
      tvtx = randomVariates.second;

    zpu->push_back(zvtx);
    
    //Loop over pileup particles
    for (int i = 0; i < _pythiaPU->event.size(); ++i) {

      if(Ignore(_pythiaPU->event[i]))
	continue;
      
      //Instantiate new pseudojet
      PseudoJet p(_pythiaPU->event[i].px(), 
		  _pythiaPU->event[i].py(), 
		  _pythiaPU->event[i].pz(),
		  _pythiaPU->event[i].e() ); 


      //extract event information
      double eta = p.rapidity();
      double sinheta = sinh(eta);
      
      //calculate eta from displacement (minEta pos)
      double zbase = z0*sgn(eta); //displace due to new location of Hard-Scatter Vertex
      double dz = zbase-zvtx;
      double hsEta = asinh((zbase-zhs)*sinheta/dz);
      
      //calculate time measured relative to if event was at 0
      double betaz=_pythiaPU->event[i].pz()/_pythiaPU->event[i].e();
      if(betaz > 1.001)
	cout << "Error: Invalid Beta value!!!" << endl;
      double time = fabs(dz/(betaz*LIGHTSPEED)); //plus random time
      time+=tvtx;
      
      double reftime = fabs((zbase-zhs)/(LIGHTSPEED*sinh(hsEta)/cosh(hsEta)));
      double corrtime = (time-reftime)*1e9;
      
      p.set_user_info(new TimingInfo(_pythiaPU->event[i].id(),_pythiaPU->event[i].charge(),
				     i,iPU,true,_pythiaPU->event[i].pT(),corrtime,time*1e9)); 
      particlesForJets.push_back(p); 
    }
    if (!_pythiaPU->next()) continue;
  }
  
  // Particle loop -----------------------------------------------------------
  for (int ip=0; ip<_pythiaHS->event.size(); ++ip){
    
    if(Ignore(_pythiaHS->event[ip]))
      continue;
    
    fastjet::PseudoJet p(_pythiaHS->event[ip].px(), 
			 _pythiaHS->event[ip].py(), 
			 _pythiaHS->event[ip].pz(),
			 _pythiaHS->event[ip].e() ); 
    
    double eta = p.rapidity();
    double sinheta = sinh(eta);

    //calculate eta from displacement (minEta pos)
    double zbase = z0*sgn(eta); //displace due to new location of Hard-Scatter Vertex    
    double dz = zbase-zhs;
    double betaz=_pythiaHS->event[ip].pz()/_pythiaHS->event[ip].e();
    if(betaz > 1.001)
      cout << "Error: Invalid Beta value!!!" << endl;
    double time = fabs(dz/(betaz*LIGHTSPEED)) + ths; //plus random time  

    double reftime = fabs((dz)/(LIGHTSPEED*sinheta/cosh(eta)));
    double corrtime = (time-reftime)*1e9;
    //0 for the primary vertex.
    p.set_user_info(new TimingInfo(_pythiaHS->event[ip].id(),_pythiaHS->event[ip].charge(),
				   ip,0, false,_pythiaHS->event[ip].pT(),corrtime,time*1e9));  
    
    particlesForJets.push_back(p);
    particlesForJets_np.push_back(p);
    
  } // end particle loop -----------------------------------------------

  JetVector selectedJets,selectedTruthJets;
  fastjet::ClusterSequenceArea clustSeq(particlesForJets, *jetDef, *active_area);
  selectJets(particlesForJets,clustSeq,selectedJets);

  fastjet::ClusterSequenceArea clustSeqTruth(particlesForJets_np, *jetDef, *active_area);
  selectedTruthJets = sorted_by_pt(clustSeqTruth.inclusive_jets(10.));

  for(unsigned int icl=0;icl<particlesForJets.size();++icl){
    j0clpt    ->push_back(particlesForJets[icl].pt());
    j0clphi   ->push_back(particlesForJets[icl].phi());
    j0cleta   ->push_back(particlesForJets[icl].eta());
    j0clm     ->push_back(particlesForJets[icl].m());
    j0clcharge->push_back(particlesForJets[icl].user_info<TimingInfo>().charge());
    j0clpdgid->push_back(particlesForJets[icl].user_info<TimingInfo>().pdg_id());
    j0cltruth->push_back(particlesForJets[icl].user_info<TimingInfo>().pileup() ? 0.0 : 1.0);
    j0clpu->push_back(particlesForJets[icl].user_info<TimingInfo>().pv());
  }//loop particles

  FillTruthTree(selectedTruthJets);
  
  tT->Fill();
  
  if(fDebug) 
    cout << "PileupAnalysis::AnalyzeEvent End " << endl;
  
  return;
}

void PileupAnalysis::selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets){
  try{
    bge->set_particles(particlesForJets);

    fastjet::Subtractor subtractor(bge.get());    
    JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
    JetVector subtractedJets = subtractor(inclusiveJets);

    JetVector allSelectedJets;
    allSelectedJets.clear();
    allSelectedJets = (*select_fwd)(subtractedJets);
    
    //select jets with pt > 20
    selectedJets.clear();
    for( auto ijet = allSelectedJets.begin(); ijet != allSelectedJets.end(); ++ijet){
      if(ijet->pt() >= 20)
	selectedJets.push_back(*ijet);
    }
  }
  catch(...){
    cerr << "Fastjet error caught in selectJets" << endl;
    exit(20);
  }
}

void PileupAnalysis::FillTruthTree(JetVector jets){  
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    truejpt->push_back(jets[ijet].pt());
    truejeta->push_back(jets[ijet].eta());
    truejphi->push_back(jets[ijet].phi());
  }
}

bool PileupAnalysis::Ignore(Pythia8::Particle &p){
  if (!p.isFinal() )      
    return true;
  switch(abs(p.id())){
  case 12:
  case 13:
  case 14:
  case 16:
    return true;
  default:
    return false;
  }
}

// declare branches
void PileupAnalysis::DeclareBranches(){
   // Event Properties 
  gROOT->ProcessLine("#include <vector>");
  
  tT->Branch("EventNumber",&fTEventNumber,"EventNumber/I");
  tT->Branch("NPV",&fTNPV,"NPV/I");
  tT->Branch("zvtxspread",&fzvtxspread,"zvtxspread/F");
  tT->Branch("tvtxspread",&ftvtxspread,"tvtxspread/F"); 
  tT->Branch("zpu","std::vector<float>",&zpu);

  tT->Branch("pt","std::vector<float>",&j0clpt);
  tT->Branch("phi","std::vector<float>",&j0clphi);
  tT->Branch("eta","std::vector<float>",&j0cleta);
  tT->Branch("m","std::vector<float>",&j0clm);
  tT->Branch("truth","std::vector<float>",&j0cltruth);
  tT->Branch("pu","std::vector<float>",&j0clpu);
  tT->Branch("charge","std::vector<float>",&j0clcharge);
  tT->Branch("pdgid","std::vector<float>",&j0clpdgid);

  tT->Branch("truejpt", "std::vector<float>",&truejpt);
  tT->Branch("truejphi","std::vector<float>",&truejphi);
  tT->Branch("truejeta","std::vector<float>",&truejeta);
  
  return;
}

// resets vars
void PileupAnalysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV = -1;
      fzvtxspread = -1;
      ftvtxspread = -1;

      zpu->clear();

      j0clpt->clear();
      j0clphi->clear();
      j0cleta->clear();
      j0clm->clear();
      j0cltruth->clear();
      j0clpu->clear();
      j0clcharge->clear();
      j0clpdgid->clear();
      truejpt->clear();
      truejphi->clear();
      truejeta->clear();
}

double TimingDistribution::probability(double zpos, double time, distribution dtype){
  double ampl=0;
  
  switch(dtype){
  case gaussian:
    return _gauss_norm*exp(-( pow(zpos,2) + pow(LIGHTSPEED*time,2) ) / (pow(_bunchsize,2)));
  case crabKissingGaussian:
    ampl=_gauss_norm*_phi_nums[1]*_psi_nums[1];
    return ampl*exp(-((pow(zpos,2)*_phi_nums[0]) + (pow(LIGHTSPEED*time,2)*_psi_nums[0])) / (pow(_bunchsize,2)));
  case pseudoRectangular:
    return _square_norm*exp(-(4*PI*PI/pow((tgamma(0.25)*_bunchsize),4))*(pow(zpos,4)+pow(LIGHTSPEED*time,4)+6*pow(LIGHTSPEED*time*zpos,2)));
  case crabKissingSquare:
    ampl=_square_norm*_phi_nums[1]*_psi_nums[1];
    return _square_norm*exp(-(4*PI*PI/pow((tgamma(0.25)*_bunchsize),4))*(pow(zpos,4)+pow(LIGHTSPEED*time,4)+6*pow(LIGHTSPEED*time*zpos,2)))*exp(-(pow(zpos*_phi,2) + pow(LIGHTSPEED*time*_psi,2)) / (pow(_bunchsize,2)));
  default:
    cerr << "Invalid RNG Distribution" << endl;
    exit(10);
  }
}

int TimingDistribution::randomSeed(){
  int timeSeed = time(NULL);
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729); 
}

TimingDistribution::TimingDistribution(float bunchsize, int seed, double phi, double psi) : _bunchsize(bunchsize){
  if(seed == -1){
    cout << "Timing Distribution Generating Random Seed" << endl;
    _seed=randomSeed();
  }
  else
    _seed=seed;

  _gauss_norm=LIGHTSPEED/(PI*_bunchsize*_bunchsize);
  _square_norm=pow(2,3.5)*LIGHTSPEED*PI/(pow(tgamma(0.25),4)*pow(_bunchsize,2));
  
  rng.seed(_seed);  
  this->phi(phi);
  this->psi(psi);
}

void TimingDistribution::phi(double phi){
  _phi=phi;
  _phi_nums[0]=1+pow(phi,2);
  _phi_nums[1]=sqrt(_phi_nums[0]);
}

void TimingDistribution::psi(double psi){
  _psi=psi;
  _psi_nums[0]=1+pow(psi,2);
  _psi_nums[1]=sqrt(_psi_nums[0]);
}

pair<double,double> TimingDistribution::get(distribution dtype){
  
  double maxprob = probability(0,0,dtype);
  double zpos;
  double time;
  
  while(true){
    zpos = uniform(-3.0*_bunchsize,3.0*_bunchsize);
    time = uniform(-3*_bunchsize/LIGHTSPEED,3*_bunchsize/LIGHTSPEED);
    if(probability(zpos,time,dtype) > maxprob*uniform())
      break;
  }
  return std::make_pair(zpos,time);
}

double TimingDistribution::uniform(double min, double max){
  static unsigned int range_min=rng.min();
  static unsigned int range_max=rng.max();

  unsigned int rn=rng();
  double drn=static_cast<double>(rn-range_min)/static_cast<double>(range_max-range_min);
  if((max != 1) and (min != 0))
    return (max-min)*drn+min;
  else
    return drn;
}
