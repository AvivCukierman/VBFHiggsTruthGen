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
Analysis::Analysis(Pythia8::Pythia *pythiaHS, Pythia8::Pythia *pythiaPU, Configuration q){

  if((pythiaHS != NULL) and (pythiaPU != NULL)){
    _pythiaHS=pythiaHS;
    _pythiaPU=pythiaPU;
  }
  else{
    cerr << "Invalid Pythia pointer passed to Analysis" << endl;
    exit(1);
  }

  fDebug=q.fDebug;
  if(fDebug) 
    cout << "Analysis::Analysis Start " << endl;

  ftest = 0;
  fOutName = q.outName;

  if(fDebug) 
    cout << "Analysis::Analysis End " << endl;
}

// Destructor 
Analysis::~Analysis(){
  if(tT != NULL){
    tT->Write();
    tF->Close();
    delete tF;
  }

  if(clpt != NULL){
    delete clpt;
    delete clphi;
    delete cleta;
    delete clm;
    delete clpdgid;

    delete j0pt;
    delete j0phi;
    delete j0eta;
    delete j0m;
  }
}

// Begin method
void Analysis::Initialize(float minEta, float maxEta, distribution dtype, int seed){
   // Declare TTree
   tF = new TFile(fOutName.c_str(), "RECREATE");
   tT = new TTree("tree", "Event Tree for VBFHiggs");
   //rnd.reset(new TimingDistribution(bunchsize,seed,phi,psi));
   _dtype=dtype;

   _minEta= (minEta > 0) ? minEta : 0;
   if(maxEta <= minEta){
     cerr << "Invalid Eta Limits " << minEta << " -> " << maxEta << "Passed to Analysis::Initialize" << endl;
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
   
   clpt = new branch();  
   clphi = new branch();  
   cleta = new branch();  
   clm = new branch();  
   clpdgid = new branch();

   j0pt = new branch();  
   j0phi = new branch();  
   j0eta = new branch();  
   j0m = new branch();  

   ResetBranches();
   
   return;
}

// Analyze
void Analysis::AnalyzeEvent(int ievt, int NPV){

  if(fDebug) 
    cout << "Analysis::AnalyzeEvent Begin " << endl;
  
  // -------------------------
  if (!_pythiaHS->next()) return;
  if(fDebug) 
    cout << "Analysis::AnalyzeEvent Event Number " << ievt << endl;
  
  // reset branches 
  ResetBranches();
  
  // new event-----------------------
  fTEventNumber = ievt;
  JetVector particlesForJets;
  //JetVector particlesForJets_np;
  
  fTNPV = NPV;

  //Pileup Loop
  //Loop over Pileup Events
  /*for (int iPU = 0; iPU < NPV; ++iPU) {
    
    //Loop over pileup particles
    for (int i = 0; i < _pythiaPU->event.size(); ++i) {

      if(Ignore(_pythiaPU->event[i]))
	continue;
      
      //Instantiate new pseudojet
      PseudoJet p(_pythiaPU->event[i].px(), 
		  _pythiaPU->event[i].py(), 
		  _pythiaPU->event[i].pz(),
		  _pythiaPU->event[i].e() ); 


      p.set_user_info(new ParticleInfo(_pythiaPU->event[i].id(),
				     i,iPU,true,_pythiaPU->event[i].pT())); 
      particlesForJets.push_back(p); 
    }
    if (!_pythiaPU->next()) continue;
  }*/
  
  // Particle loop -----------------------------------------------------------
  for (int ip=0; ip<_pythiaHS->event.size(); ++ip){
    
    if(Ignore(_pythiaHS->event[ip]))
      continue;
    
    fastjet::PseudoJet p(_pythiaHS->event[ip].px(), 
			 _pythiaHS->event[ip].py(), 
			 _pythiaHS->event[ip].pz(),
			 _pythiaHS->event[ip].e() ); 
    
    p.set_user_info(new ParticleInfo(_pythiaHS->event[ip].id(),
				   ip,0,false,_pythiaHS->event[ip].pT()));  
    
    particlesForJets.push_back(p);
    //particlesForJets_np.push_back(p);
    
  } // end particle loop -----------------------------------------------

  //JetVector selectedJets,selectedTruthJets;
  JetVector selectedJets;
  fastjet::ClusterSequenceArea clustSeq(particlesForJets, *jetDef, *active_area);
  selectJets(particlesForJets,clustSeq,selectedJets);

  //fastjet::ClusterSequenceArea clustSeqTruth(particlesForJets_np, *jetDef, *active_area);
  //selectedTruthJets = sorted_by_pt(clustSeqTruth.inclusive_jets(10.));

  for(unsigned int icl=0;icl<particlesForJets.size();++icl){
    clpt    ->push_back(particlesForJets[icl].pt());
    clphi   ->push_back(particlesForJets[icl].phi());
    cleta   ->push_back(particlesForJets[icl].eta());
    clm     ->push_back(particlesForJets[icl].m());
    clpdgid->push_back(particlesForJets[icl].user_info<ParticleInfo>().pdg_id());
  }//loop particles

  for(unsigned int ijet=0; ijet<selectedJets.size(); ++ijet){
    fastjet::PseudoJet jet = selectedJets[ijet];
    j0pt  ->push_back(jet.perp());
    j0eta ->push_back(jet.eta());
    j0phi ->push_back(jet.phi());
    j0m   ->push_back(jet.m());
  }

  //FillTruthTree(selectedTruthJets);
  
  tT->Fill();
  
  if(fDebug) 
    cout << "Analysis::AnalyzeEvent End " << endl;
  
  return;
}

void Analysis::selectJets(JetVector &particlesForJets, fastjet::ClusterSequenceArea &clustSeq, JetVector &selectedJets){
  try{
    bge->set_particles(particlesForJets);

    fastjet::Subtractor subtractor(bge.get());    
    JetVector inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
    JetVector subtractedJets = subtractor(inclusiveJets);

    JetVector allSelectedJets;
    allSelectedJets.clear();
    allSelectedJets = (*select_fwd)(subtractedJets);
    
    //select jets with pt > 10
    selectedJets.clear();
    for( auto ijet = allSelectedJets.begin(); ijet != allSelectedJets.end(); ++ijet){
      if(ijet->pt() >= 10)
	selectedJets.push_back(*ijet);
    }
  }
  catch(...){
    cerr << "Fastjet error caught in selectJets" << endl;
    exit(20);
  }
}

/*void Analysis::FillTruthTree(JetVector jets){  
  for (unsigned int ijet=0; ijet<jets.size(); ijet++){
    truejpt->push_back(jets[ijet].pt());
    truejeta->push_back(jets[ijet].eta());
    truejphi->push_back(jets[ijet].phi());
  }
}*/

bool Analysis::Ignore(Pythia8::Particle &p){
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
void Analysis::DeclareBranches(){
   // Event Properties 
  gROOT->ProcessLine("#include <vector>");
  
  tT->Branch("EventNumber",&fTEventNumber,"EventNumber/I");
  tT->Branch("NPV",&fTNPV,"NPV/I");
  //tT->Branch("zvtxspread",&fzvtxspread,"zvtxspread/F");
  //tT->Branch("tvtxspread",&ftvtxspread,"tvtxspread/F"); 
  //tT->Branch("zpu","std::vector<float>",&zpu);

  tT->Branch("clpt","std::vector<float>",&clpt);
  tT->Branch("clphi","std::vector<float>",&clphi);
  tT->Branch("cleta","std::vector<float>",&cleta);
  tT->Branch("clm","std::vector<float>",&clm);
  tT->Branch("clpdgid","std::vector<float>",&clpdgid);

  tT->Branch("j0pt","std::vector<float>",&j0pt);
  tT->Branch("j0phi","std::vector<float>",&j0phi);
  tT->Branch("j0eta","std::vector<float>",&j0eta);
  tT->Branch("j0m","std::vector<float>",&j0m);

  //tT->Branch("truth","std::vector<float>",&j0cltruth);
  //tT->Branch("pu","std::vector<float>",&j0clpu);
  //tT->Branch("charge","std::vector<float>",&j0clcharge);
  
  return;
}

// resets vars
void Analysis::ResetBranches(){
      // reset branches 
      fTEventNumber                 = -999;
      fTNPV = -1;

      clpt->clear();
      clphi->clear();
      cleta->clear();
      clm->clear();
      clpdgid->clear();

      j0pt->clear();
      j0phi->clear();
      j0eta->clear();
      j0m->clear();
}
