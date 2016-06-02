from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram
import os
import numpy
from optparse import OptionParser
from sys import stdout,argv
import ROOT as r
import pdb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
#plt.style.use('atlas')
import matplotlib.mlab as mlab

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="Directory containing input files",type=str, default="data")
parser.add_option("--inputFile", help="Input file name",type=str, default="VBFHiggs.root")
parser.add_option("--submitDir", help="Directory containing output files",type=str, default="output")
parser.add_option("--plotDir", help="Directory containing plots",type=str, default="plots")
parser.add_option("--numEvents", help="How many events to include (set to -1 for all events)",type=int, default=100000)
parser.add_option("--identifier", help="Identify dataset",type=str, default="VBFHiggs")

# Root configuration
parser.add_option("--jetpt", help="jet pT branch name",type=str, default="j0pt")
parser.add_option("--jeteta", help="jet eta branch name",type=str, default="j0eta")
parser.add_option("--jetphi", help="jet phi branch name",type=str, default="j0phi")
parser.add_option("--jetm", help="jet m branch name",type=str, default="j0m")

parser.add_option("--gammapt", help="gamma pT branch name",type=str, default="gammapt")
parser.add_option("--gammaeta", help="gamma eta branch name",type=str, default="gammaeta")
parser.add_option("--gammaphi", help="gamma phi branch name",type=str, default="gammaphi")
parser.add_option("--gammam", help="gamma m branch name",type=str, default="gammam")

# object configuration
parser.add_option("--minjetpt", help="min pt cut on jets", type=float, default=0)
parser.add_option("--mingammapt", help="min pt cut on photons", type=float, default=0)

(options, args) = parser.parse_args()

if not os.path.exists(options.inputDir): raise OSError(options.inputDir +' does not exist. This is where the input Root files go.')
if not os.path.exists(options.submitDir):
  print '== Making folder '+options.submitDir+' =='
  os.makedirs(options.submitDir)
if not os.path.exists(options.plotDir):
  print '== Making folder '+options.plotDir+' =='
  os.makedirs(options.plotDir)

def readRoot():
  import glob

  filenames = glob.glob(options.inputDir+'/'+options.inputFile)
  if len(filenames) == 0: raise OSError('Can\'t find file '+options.inputDir+'/'+options.inputFile) 
  filename = filenames[0]
  print '== Reading in '+filename+' =='

  tree = r.TChain('tree')
  tree.Add(filename) 

  # make sure the branches are compatible between the two
  branches = set(i.GetName() for i in tree.GetListOfBranches())

  # required:
  if options.jetpt not in branches: raise RuntimeError(options.jetpt+' branch does not exist. This is the branch containing jet pTs.')
  else: print '== \''+options.jetpt+'\' branch is being read as jet pTs =='
  if options.jeteta not in branches: raise RuntimeError(options.jeteta+' branch does not exist. This is the branch containing jet etas.')
  else: print '== \''+options.jeteta+'\' branch is being read as jet etas =='
  if options.jetphi not in branches: raise RuntimeError(options.jetphi+' branch does not exist. This is the branch containing jet phis.')
  else: print '== \''+options.jetphi+'\' branch is being read as jet phis =='
  if options.jetm not in branches: raise RuntimeError(options.jetm+' branch does not exist. This is the branch containing jet ms.')
  else: print '== \''+options.jetm+'\' branch is being read as jet ms =='

  if options.gammapt not in branches: raise RuntimeError(options.gammapt+' branch does not exist. This is the branch containing gamma pTs.')
  else: print '== \''+options.gammapt+'\' branch is being read as gamma pTs =='
  if options.gammaeta not in branches: raise RuntimeError(options.gammaeta+' branch does not exist. This is the branch containing gamma etas.')
  else: print '== \''+options.gammaeta+'\' branch is being read as gamma etas =='
  if options.gammaphi not in branches: raise RuntimeError(options.gammaphi+' branch does not exist. This is the branch containing gamma phis.')
  else: print '== \''+options.gammaphi+'\' branch is being read as gamma phis =='
  if options.gammam not in branches: raise RuntimeError(options.gammam+' branch does not exist. This is the branch containing gamma ms.')
  else: print '== \''+options.gammam+'\' branch is being read as gamma ms =='

  nentries = tree.GetEntries()

  jetpts = []
  jetetas = []
  jetphis = []
  jetms = []

  gammapts = []
  gammaetas = []
  gammaphis = []
  gammams = []

  for jentry in xrange(nentries):
    if jentry>options.numEvents and options.numEvents>0: continue
    tree.GetEntry(jentry)

    if not jentry%1000:
      stdout.write('== \r%d events read ==\n'%jentry)
      stdout.flush() 
    
    jpts = array([j for j in getattr(tree,options.jetpt)])
    jetas = array([j for j in getattr(tree,options.jeteta)])
    jphis = array([j for j in getattr(tree,options.jetphi)])
    jms = array([j for j in getattr(tree,options.jetm)])

    jetas = jetas[jpts>options.minjetpt]
    jphis = jphis[jpts>options.minjetpt]
    jms = jms[jpts>options.minjetpt]
    jpts = jpts[jpts>options.minjetpt]

    gpts = array([g for g in getattr(tree,options.gammapt)])
    getas = array([g for g in getattr(tree,options.gammaeta)])
    gphis = array([g for g in getattr(tree,options.gammaphi)])
    gms = array([g for g in getattr(tree,options.gammam)])

    getas = getas[gpts>options.mingammapt]
    gphis = gphis[gpts>options.mingammapt]
    gms = gms[gpts>options.mingammapt]
    gpts = gpts[gpts>options.mingammapt]

    jetpts.append(jpts)
    jetetas.append(jetas)
    jetphis.append(jphis)
    jetms.append(jms)

    gammapts.append(gpts)
    gammaetas.append(getas)
    gammaphis.append(gphis)
    gammams.append(gms)

  return jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams

def plot(jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams):
  jetmults = [len(j) for j in jetpts]
  gammamults = [len(g) for g in gammapts]

  n,bins,patches = plt.hist(jetmults,normed=True,bins=max(jetmults),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Jet Multiplicity')
  plt.ylabel('a.u.')
  plt.savefig(options.plotDir+'/jetmultiplicity_'+options.identifier+'.png')
  plt.close()

  n,bins,patches = plt.hist(gammamults,normed=True,bins=max(gammamults),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Photon Multiplicity')
  plt.ylabel('a.u.')
  plt.savefig(options.plotDir+'/gammamultiplicity_'+options.identifier+'.png')
  plt.close()

  jetamasses = []
  jetamasses_incl = []
  jetavecs = []
  for ept,eeta,ephi,em,emult in zip(jetpts,jetetas,jetphis,jetms,jetmults):
    if(emult<2):
      jetamasses_incl.append(-1)
    else:
      v1 = r.TLorentzVector()
      v2 = r.TLorentzVector()

      v1.SetPtEtaPhiM(ept[0],eeta[0],ephi[0],em[0])
      v2.SetPtEtaPhiM(ept[1],eeta[1],ephi[1],em[1])

      va = v1 + v2
      jetavecs.append(va)
      jetamasses.append(va.M())
      jetamasses_incl.append(va.M())

  n,bins,patches = plt.hist(jetamasses,normed=True,bins=20,facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel(r'$a\rightarrow jj$ mass')
  plt.ylabel('a.u.')
  plt.savefig(options.plotDir+'/ajjmass_'+options.identifier+'.png')
  plt.close()

  gammaamasses = []
  gammaamasses_incl = []
  gammaavecs = []
  for ept,eeta,ephi,em,emult in zip(gammapts,gammaetas,gammaphis,gammams,gammamults):
    if(emult<2):
      gammaamasses_incl.append(-1)
    else:
      v1 = r.TLorentzVector()
      v2 = r.TLorentzVector()

      v1.SetPtEtaPhiM(ept[0],eeta[0],ephi[0],em[0])
      v2.SetPtEtaPhiM(ept[1],eeta[1],ephi[1],em[1])

      va = v1 + v2
      gammaavecs.append(va)
      gammaamasses.append(va.M())
      gammaamasses_incl.append(va.M())

  n,bins,patches = plt.hist(gammaamasses,normed=True,bins=20,facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel(r'$a\rightarrow \gamma\gamma$ mass')
  plt.ylabel('a.u.')
  plt.savefig(options.plotDir+'/agamgammass_'+options.identifier+'.png')
  plt.close()




jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams = readRoot()
plot(jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams)

