from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram
from pickle import dump
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
parser.add_option("--inputFile", help="Input file name",type=str, default="")
parser.add_option("--submitDir", help="Directory containing output files",type=str, default="output")
parser.add_option("--plotDir", help="Directory containing plots",type=str, default="plots")
parser.add_option("--numEvents", help="How many events to include (set to -1 for all events)",type=int, default=-1)
parser.add_option("--identifier", help="Identify dataset",type=str, default="VBFHiggs")

# Root configuration
parser.add_option("--jetpt", help="jet pT branch name",type=str, default="j0pt")
parser.add_option("--jeteta", help="jet eta branch name",type=str, default="j0eta")
parser.add_option("--jetphi", help="jet phi branch name",type=str, default="j0phi")
parser.add_option("--jetm", help="jet m branch name",type=str, default="j0m")
parser.add_option("--jetid", help="jet id branch name",type=str, default="j0id")

parser.add_option("--gammapt", help="gamma pT branch name",type=str, default="gammapt")
parser.add_option("--gammaeta", help="gamma eta branch name",type=str, default="gammaeta")
parser.add_option("--gammaphi", help="gamma phi branch name",type=str, default="gammaphi")
parser.add_option("--gammam", help="gamma m branch name",type=str, default="gammam")

parser.add_option("--Ntruthphotons", help="Ntruthphotons branch name",type=str, default="Ntruthphotons")

# object configuration
parser.add_option("--minjetpt", help="min pt cut on jets", type=float, default=0)
parser.add_option("--mingammapt", help="min pt cut on photons", type=float, default=0)
parser.add_option("--maxjeteta", help="max eta cut on jets", type=float, default=2.5)
parser.add_option("--maxgammaeta", help="max eta cut on photons", type=float, default=2.5)

(options, args) = parser.parse_args()

HIGGS_MASS = 125

if not os.path.exists(options.inputDir): raise OSError(options.inputDir +' does not exist. This is where the input Root files go.')
if not os.path.exists(options.submitDir):
  print '== Making folder '+options.submitDir+' =='
  os.makedirs(options.submitDir)
if not os.path.exists(options.plotDir):
  print '== Making folder '+options.plotDir+' =='
  os.makedirs(options.plotDir)

def dr(eta1,phi1,eta2,phi2):
  v1 = r.TLorentzVector()
  v2 = r.TLorentzVector()

  #pt,m don't matter for deltaR
  v1.SetPtEtaPhiM(1,eta1,phi1,1)
  v2.SetPtEtaPhiM(1,eta2,phi2,1)

  return v1.DeltaR(v2)

def dphi(phi1,phi2):
  v1 = r.TLorentzVector()
  v2 = r.TLorentzVector()

  #pt,m,eta don't matter for deltaR
  v1.SetPtEtaPhiM(1,0,phi1,1)
  v2.SetPtEtaPhiM(1,0,phi2,1)

  return v1.DeltaR(v2)

def mindr(etas1,phis1,etas2,phis2,same):
  if not len(etas1)==len(phis1): raise RuntimeError('etas and phis different lengths')
  if not len(etas2)==len(phis2): raise RuntimeError('etas and phis different lengths')
  if same and len(etas1)<2: raise RuntimeError('not enough etas and phis to calculate mindr')
  cmindr = float('inf')
  for i,(eta1,phi1) in enumerate(zip(etas1,phis1)):
    if same:
      for eta2,phi2 in zip(etas1[i+1:],phis1[i+1:]):
        cdr = dr(eta1,phi1,eta2,phi2)
        cmindr = min(cmindr,cdr)
    else:
      for eta2,phi2 in zip(etas2,phis2):
        cdr = dr(eta1,phi1,eta2,phi2)
        cmindr = min(cmindr,cdr)
  return cmindr

def combine_vecs(pts,etas,phis,ms):
    v1 = r.TLorentzVector()
    v2 = r.TLorentzVector()

    v1.SetPtEtaPhiM(pts[0],etas[0],phis[0],ms[0])
    v2.SetPtEtaPhiM(pts[1],etas[1],phis[1],ms[1])

    va = v1 + v2
    return va

def readRoot():
  import glob
  global cutflow

  if len(options.inputFile)>0:
    filenames = glob.glob(options.inputDir+'/'+options.inputFile)
    if len(filenames) == 0: raise OSError('Can\'t find file '+options.inputDir+'/'+options.inputFile) 
  else:
    filenames = glob.glob(options.inputDir+'/*.root')
    if len(filenames) == 0: raise OSError('Can\'t find files in '+options.inputDir) 
  tree = r.TChain('tree')
  for filename in filenames:
    statinfo = os.stat(filename)
    if statinfo.st_size < 10000: continue #sometimes batch jobs fail
    print '== Reading in '+filename+' =='

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
  if options.jetid not in branches: raise RuntimeError(options.jetid+' branch does not exist. This is the branch containing jet ids.')
  else: print '== \''+options.jetid+'\' branch is being read as jet ids =='

  if options.gammapt not in branches: raise RuntimeError(options.gammapt+' branch does not exist. This is the branch containing gamma pTs.')
  else: print '== \''+options.gammapt+'\' branch is being read as gamma pTs =='
  if options.gammaeta not in branches: raise RuntimeError(options.gammaeta+' branch does not exist. This is the branch containing gamma etas.')
  else: print '== \''+options.gammaeta+'\' branch is being read as gamma etas =='
  if options.gammaphi not in branches: raise RuntimeError(options.gammaphi+' branch does not exist. This is the branch containing gamma phis.')
  else: print '== \''+options.gammaphi+'\' branch is being read as gamma phis =='
  if options.gammam not in branches: raise RuntimeError(options.gammam+' branch does not exist. This is the branch containing gamma ms.')
  else: print '== \''+options.gammam+'\' branch is being read as gamma ms =='

  if options.Ntruthphotons not in branches: raise RuntimeError(options.Ntruthphotons+' branch does not exist. This is the branch containing the number of truth photons in the hardest process.')
  else: print '== \''+options.Ntruthphotons+'\' branch is being read as number of truth photons in the hardest process =='

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
    cutflow[0]+=1 #entries read in -> cross section
    tree.GetEntry(jentry)

    if not jentry%1000:
      stdout.write('== \r%d events read ==\n'%jentry)
      stdout.flush() 

    Ntruthphotons = getattr(tree,options.Ntruthphotons)
    if not Ntruthphotons == 2: continue #only events with 2 truth photons
    cutflow[1]+=1

    treejpts = getattr(tree,options.jetpt)
    treejetas = getattr(tree,options.jeteta)
    treejphis = getattr(tree,options.jetphi)
    treejms = getattr(tree,options.jetm)
    treejids = getattr(tree,options.jetid)
    
    jpts = []
    jetas = []
    jphis = []
    jms = []

    for jpt,jeta,jphi,jm,jid in zip(treejpts,treejetas,treejphis,treejms,treejids):
      if jid==22: continue #don't count jets formed from photons
      jpts.append(jpt)
      jetas.append(jeta)
      jphis.append(jphi)
      jms.append(jm)

    jpts = array(jpts)
    jetas = array(jetas)
    jphis = array(jphis)
    jms = array(jms)

    good_indices = all([jpts>options.minjetpt,abs(jetas)<options.maxjeteta],axis=0)
    jetas = jetas[good_indices]
    jphis = jphis[good_indices]
    jms = jms[good_indices]
    jpts = jpts[good_indices]

    gpts = array([g for g in getattr(tree,options.gammapt)])
    getas = array([g for g in getattr(tree,options.gammaeta)])
    gphis = array([g for g in getattr(tree,options.gammaphi)])
    gms = array([g for g in getattr(tree,options.gammam)])

    good_indices = all([gpts>options.mingammapt,abs(getas)<options.maxgammaeta],axis=0)
    getas = getas[good_indices]
    gphis = gphis[good_indices]
    gms = gms[good_indices]
    gpts = gpts[good_indices]

    jetpts.append(jpts)
    jetetas.append(jetas)
    jetphis.append(jphis)
    jetms.append(jms)

    gammapts.append(gpts)
    gammaetas.append(getas)
    gammaphis.append(gphis)
    gammams.append(gms)

  return jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams

def initial_cuts(jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams):
  global cutflow,hists
  jetmults = [len(j) for j in jetpts]
  gammamults = [len(g) for g in gammapts]

  n,bins,patches = plt.hist(jetmults,normed=True,bins=range(0,max(jetmults)+1,1),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Jet Multiplicity')
  plt.ylabel('a.u.')
  plotname = 'jetmultiplicity'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(gammamults,normed=True,bins=range(0,max(gammamults)+1,1),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Photon Multiplicity')
  plt.ylabel('a.u.')
  plotname = 'gammamultiplicity'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  jjmindrs = []
  ggmindrs = []
  gjmindrs = []
  good_indices = []

  newjetpts = []
  newjetetas = []
  newjetphis = []
  newjetms = []

  newgammapts = []
  newgammaetas = []
  newgammaphis = []
  newgammams = []

  for i,(jpt,jeta,jphi,jm,jmult,gpt,geta,gphi,gm,gmult) in enumerate(zip(jetpts,jetetas,jetphis,jetms,jetmults,gammapts,gammaetas,gammaphis,gammams,gammamults)):
    if not gmult==2:
      continue
    cutflow[2]+=1
    if jmult<2:
      continue
    cutflow[3]+=1
    jjmindr = mindr(jeta,jphi,jeta,jphi,True)
    jjmindrs.append(jjmindr)
    ggmindr = mindr(geta,gphi,geta,gphi,True)
    ggmindrs.append(ggmindr)
    gjmindr = mindr(geta,gphi,jeta,jphi,False)
    gjmindrs.append(gjmindr)
    if jjmindr>0.4 and ggmindr>0.4 and gjmindr>0.4:
      newjetpts.append(array([jpt[0],jpt[1]]))
      newjetetas.append(array([jeta[0],jeta[1]]))
      newjetphis.append(array([jphi[0],jphi[1]]))
      newjetms.append(array([jm[0],jm[1]]))

      newgammapts.append(array([gpt[0],gpt[1]]))
      newgammaetas.append(array([geta[0],geta[1]]))
      newgammaphis.append(array([gphi[0],gphi[1]]))
      newgammams.append(array([gm[0],gm[1]]))
  cutflow[4] = len(newjetpts)
  
  gam1pts = [newgammapt[0] for newgammapt in newgammapts]
  n,bins,patches = plt.hist(gam1pts,normed=True,bins=range(0,105,5),facecolor='b',histtype='stepfilled')
  plt.xlim(0,100)
  plt.ylim(0,max(n)*1.1)
  plt.xlabel('Leading $\gamma$ $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'gam1pt'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  gam2pts = [newgammapt[1] for newgammapt in newgammapts]
  n,bins,patches = plt.hist(gam2pts,normed=True,bins=range(0,105,5),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,100)
  plt.xlabel('Subleading $\gamma$ $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'gam2pt'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(jjmindrs,normed=True,bins=numpy.arange(0,3.15,.15),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum jet-jet $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'jjmindr'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(ggmindrs,normed=True,bins=numpy.arange(0,3.15,.15),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum $\gamma-\gamma$ $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'ggmindr'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(gjmindrs,normed=True,bins=numpy.arange(0,3.15,.15),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum jet-$\gamma$ $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'gjmindr'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  return array(newjetpts),array(newjetetas),array(newjetphis),array(newjetms),array(newgammapts),array(newgammaetas),array(newgammaphis),array(newgammams)

def final_cuts(jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams):
  global cutflow,hists
  jetdphis = []
  jetamasses = []
  jetavecs = []
  for ept,eeta,ephi,em in zip(jetpts,jetetas,jetphis,jetms):
    jetdphis.append(dphi(ephi[0],ephi[1]))
    va = combine_vecs(ept,eeta,ephi,em)
    jetavecs.append(va)
    jetamasses.append(va.M())

  binwidth = 5
  n,bins,patches = plt.hist(jetamasses,normed=True,bins=numpy.arange(0, 60 + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,60)
  plt.xlabel(r'$M_{jj}$')
  plt.ylabel('a.u.')
  plotname = 'mjj'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(jetdphis,normed=True,bins=numpy.arange(0,4.2,0.2),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.5)
  plt.xlabel(r'$\Delta\phi_{jj}$')
  plt.ylabel('a.u.')
  plotname = 'dphijj'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  gammadphis = []
  gammaamasses = []
  gammaavecs = []
  for ept,eeta,ephi,em in zip(gammapts,gammaetas,gammaphis,gammams):
    gammadphis.append(dphi(ephi[0],ephi[1]))
    va = combine_vecs(ept,eeta,ephi,em)
    gammaavecs.append(va)
    gammaamasses.append(va.M())

  binwidth = 5
  n,bins,patches = plt.hist(gammaamasses,normed=True,bins=numpy.arange(0, 60 + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,60)
  plt.xlabel(r'$M_{\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'mgamgam'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  n,bins,patches = plt.hist(gammadphis,normed=True,bins=numpy.arange(0,4.2,0.2),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.5)
  plt.xlabel(r'$\Delta\phi_{\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'dphigamgam'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  jetdphis = array(jetdphis)
  gammadphis = array(gammadphis)
  jetamasses = array(jetamasses)
  gammaamasses = array(gammaamasses)

  data = abs(jetamasses-gammaamasses)
  binwidth = 5
  n,bins,patches = plt.hist(data,normed=True,bins=numpy.arange(0, 50 + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,50)
  plt.xlabel(r'$|M_{jj}-M_{\gamma\gamma}|$')
  plt.ylabel('a.u.')
  plotname = 'amassdiff'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  good_indices = all([jetdphis<1.0,gammadphis<1.3,abs(jetamasses-gammaamasses)<15],axis=0)
  cutflow[5] = sum(good_indices)

  hmasses = array([(va1+va2).M() for va1,va2 in zip(gammaavecs,jetavecs)])
  binwidth = 5
  n,bins,patches = plt.hist(hmasses,normed=True,bins=numpy.arange(0, 200 + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,200)
  plt.xlabel(r'$M_{jj\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'mh'
  plt.savefig(options.plotDir+'/'+plotname+'_'+options.identifier+'.png')
  plt.close()
  hists[plotname] = (n,bins)

  good_indices = all([jetdphis<1.0,gammadphis<1.3,abs(jetamasses-gammaamasses)<15,abs(hmasses-HIGGS_MASS)<25],axis=0)
  cutflow[6] = sum(good_indices)

  return

hists = {}
cutflow = [0,0,0,0,0,0,0]
jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams = readRoot()
ijetpts,ijetetas,ijetphis,ijetms,igammapts,igammaetas,igammaphis,igammams = initial_cuts(jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams)
final_cuts(ijetpts,ijetetas,ijetphis,ijetms,igammapts,igammaetas,igammaphis,igammams)
with open(options.submitDir+'/cutflow_'+options.identifier+'.p','wb') as outfile:
  dump(cutflow,outfile)
with open(options.submitDir+'/hists_'+options.identifier+'.p','wb') as outfile:
  dump(hists,outfile)
