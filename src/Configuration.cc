#include "Configuration.h"

Configuration::Configuration(int argc, char* argv[]){

    HSmode =smearMode::ZT;
    PUmode =smearMode::ZT;
    useCK     =false;
    filterCharge=true;
    magfield = false;
    storeallparticles = false;
    int profile;

    po::options_description gen_desc("Allowed options");
    gen_desc.add_options()
      ("help", "produce help message")
      ("Debug",     po::value<int>(&fDebug) ->default_value(0) ,     "Debug flag")
      ("OutFile",   po::value<string>(&outName)->default_value("VBFHiggs.root"), "output file name")
      ("Seed",      po::value<int>(&seed)->default_value(-1), "Seed. -1 means random seed");

    po::options_description sim_flag("Simulation Flags");
    sim_flag.add_options()
      ("GenericFlag",   "Generic flag. Doesn't do anything.");

    po::options_description sim_desc("Simulation Settings");
    sim_desc.add_options()
      ("nevents",   po::value<int>(&nEvents)->default_value(1) ,    "Number of Events ")
      ("Pileup",    po::value<int>(&pileup)->default_value(0), "Number of Additional Interactions.")
      ("BunchSize", po::value<float>(&bunchsize)->default_value(0.075), "Size of Proton Bunches")
      ("MinEta",    po::value<float>(&minEta)->default_value(0.0), "Minimum (abs) Pseudorapidity for Particles")
      ("MaxEta",    po::value<float>(&maxEta)->default_value(5.0), "Minimum (abs) Pseudorapidity for Particles")
      ("Proc",      po::value<int>(&proc)->default_value(2), "Process:\n - 1: VBFHiggs->xxyy (x,y=gluon,photon)\n  - 2: VBF H->inv")
      //("pThatMin",  po::value<float>(&pThatmin)->default_value(100), "pThatMin for QCD")
      //("pThatMax",  po::value<float>(&pThatmax)->default_value(500), "pThatMax for QCD")
      ("ScalarMass", po::value<float>(&scalar_mass)->default_value(30), "a mass (GeV) [H->aa]");
    
    po::options_description desc;
    desc.add(gen_desc).add(sim_flag).add(sim_desc);
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")>0){
        cout << desc << "\n";
	exit(0);
    }

    print();

    cout << "\t";

    seed=getSeed(seed);
        
    cout << endl;
}

void Configuration::print(){

  cout << endl;
  cout << "=================================================================" << endl;
  cout << "=                        VBFHiggs Analysis                      =" << endl;
  cout << "=================================================================" << endl << endl;

  cout << "Settings:" << endl;
  for (po::variables_map::const_iterator itr=vm.begin();itr != vm.end();++itr){
    printf("%15s\t",itr->first.c_str());
    
    try { 
      cout << "= " << itr->second.as<double>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<float>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<int>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
    try { 
      cout << "= " << itr->second.as<std::string>() << std::endl;
      continue;
    } catch(...) {/* do nothing */ }
  }
}

void Configuration::ConfigurePythiaSignal(Pythia8::Pythia* hs){

  hs->readString("Print:quiet=on");
  hs->readString("Random:setSeed = on"); 
  std::stringstream ss; 
  ss << "Random:seed = " << seed;
  hs->readString(ss.str());
 
  //ggH2aa2qqgammagamma
  if(proc == 1){
    hs->readString("Higgs:useBSM = on");
    hs->readString("HiggsBSM:gg2H1 = on");
    ss.str("");
    ss << "35:m0 = " << H_MASS;
    hs->readString(ss.str());
    ss.str("");
    ss << "35:mWidth = " << H_WIDTH;
    hs->readString(ss.str());
    hs->readString("35:doForceWidth = on");
    hs->readString("35:onMode = off");
    hs->readString("35:oneChannel = 1 1 100 36 36");
    hs->readString("35:onIfAny = 36"); //h->aa
    hs->readString("36:onMode = off");
    hs->readString("36:oneChannel = 1 0.5 100 22 22"); //a->gammagamma
    hs->readString("36:addChannel = 1 0.5 100 21 21"); //a->gg
    hs->readString("36:onIfAny = 21 22");
    ss.str("");
    ss << "36:m0 = " << scalar_mass;
    hs->readString(ss.str());
    hs->readString("36:mWidth = 0.01"); //narrow width
    hs->readString("36:mMin = 29.5");
    hs->readString("36:mMax = 30.5");
    hs->readString("36:tau0 = 0"); //scalar lifetime
    hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
  }
  //ffbar2gammagamma
  else if(proc == 2){
    hs->readString("PromptPhoton:ffbar2gammagamma = on");
    //hs->readString("TimeShower:weakShower = on");
    hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
  }
  else if(proc == 3){
    hs->readString("Higgs:useBSM = on");
    hs->readString("HiggsBSM:ff2H1ff(t:ZZ) = on");
    hs->readString("HiggsBSM:ff2H1ff(t:WW) = on");
    ss.str("");
    ss << "35:m0 = " << H_MASS;
    hs->readString(ss.str());
    ss.str("");
    ss << "35:mWidth = " << H_WIDTH;
    hs->readString(ss.str());
    hs->readString("35:doForceWidth = on");
    hs->readString("35:onMode = off");
    hs->readString("35:oneChannel = 1 1 100 36 36");
    hs->readString("35:onIfAny = 36"); //h->aa
    hs->readString("36:onMode = off");
    hs->readString("36:oneChannel = 1 0.5 100 22 22"); //a->gammagamma
    hs->readString("36:addChannel = 1 0.5 100 21 21"); //a->gg
    hs->readString("36:onIfAny = 21 22");
    ss.str("");
    ss << "36:m0 = " << scalar_mass;
    hs->readString(ss.str());
    hs->readString("36:mWidth = 0.01"); //narrow width
    hs->readString("36:mMin = 29.5");
    hs->readString("36:mMax = 30.5");
    hs->readString("36:tau0 = 0"); //scalar lifetime
    hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
  }
  else if(proc == 4){
    hs->readString("HiggsSM:ff2Hff(t:ZZ) = on");
    hs->readString("HiggsSM:ff2Hff(t:WW) = on");
    hs->readString("25:m0 = 125");
    hs->readString("25:onMode = off");
    hs->readString("25:onIfAny = 23");
    hs->readString("23:onMode = off");
    hs->readString("23:onIfAny = 12 14 16");
    hs->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */); //this has to be the last line!
  }
  else{ 
    throw std::invalid_argument("received invalid 'process'");
  }
  return;
}

void Configuration::ConfigurePythiaPileup(Pythia8::Pythia* pu){
  //Setup the pileup
  pu->readString("Random:setSeed = on");   
  std::stringstream ss;
  ss.clear(); 
  ss.str(""); 
  ss << "Random:seed = " << seed+1; 
  pu->readString(ss.str());
  pu->readString("Print:quiet=on");
  pu->readString("SoftQCD:nonDiffractive = on");
  pu->readString("HardQCD:all = off");
  pu->readString("PhaseSpace:pTHatMin  = .1");
  pu->readString("PhaseSpace:pTHatMax  = 20000");
  pu->init(2212 /* p */, 2212 /* p */, 14000. /* TeV */);
  return;
}

int Configuration::getSeed(int seed){                                                      
  if (seed > -1) return seed;
  int timeSeed = time(NULL);                                                                 
  return abs(((timeSeed*181)*((getpid()-83)*359))%104729);
}
