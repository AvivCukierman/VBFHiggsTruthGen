from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram,cosh,append
import json
from collections import namedtuple
from pickle import dump,load
import os
import numpy
from optparse import OptionParser
from sys import stdout,argv
import ROOT as r
import pdb
from helper_functions import *

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
parser.add_option("-r","--root", help="force reading in root files",action="store_true", default=False)
parser.add_option("--analysis_label", help="identifier of analysis",type=str, default="")

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
#parser.add_option("--minjetpt", help="min pt cut on any jet in event", type=float, default=0)
#parser.add_option("--minjet1pt", help="min pt cut on leading jet", type=float, default=0)
#parser.add_option("--minjet2pt", help="min pt cut on subleading jet", type=float, default=0)
#parser.add_option("--mingammapt", help="min pt cut on any photon in event", type=float, default=0)
#parser.add_option("--mingamma1pt", help="min pt cut on leading photon", type=float, default=0)
#parser.add_option("--mingamma2pt", help="min pt cut on subleading photon", type=float, default=0)
#parser.add_option("--maxjeteta", help="max eta cut on jets", type=float, default=2.5)
#parser.add_option("--maxVBFjeteta", help="max eta cut on jets for VBF analysis", type=float, default=4.5)
#parser.add_option("--maxgammaeta", help="max eta cut on photons", type=float, default=2.5)
#parser.add_option("-s","--smear", help="apply smearing to object pTs",action="store_true", default=False)
#parser.add_option("-m","--minmjj", help="min mjj (for VBF jets)",type=float, default=-1)

(options, args) = parser.parse_args()


HIGGS_MASS = 125

if not os.path.exists(options.inputDir): raise OSError(options.inputDir +' does not exist. This is where the input Root files go.')
if not os.path.exists(options.submitDir):
  print '== Making folder '+options.submitDir+' =='
  os.makedirs(options.submitDir)
if not os.path.exists(options.plotDir):
  print '== Making folder '+options.plotDir+' =='
  os.makedirs(options.plotDir)

analysis_label=options.analysis_label
analysis_file_path = options.submitDir+'/'+analysis_label+'.json'
if not os.path.exists(analysis_file_path):  raise OSError(analysis_file+' does not exist. Contains analysis details.')
with open(analysis_file_path) as analysis_file:
  analysis = json.load(analysis_file)

#read in options
options.minjetpt = analysis['minjetpt']
options.mingammapt = analysis['mingammapt']
options.minjet1pt = analysis['minjet1pt']
options.minjet2pt = analysis['minjet2pt']
options.minjet3pt = analysis['minjet3pt']
options.minjet4pt = analysis['minjet4pt']
options.mingamma1pt = analysis['mingamma1pt']
options.mingamma2pt = analysis['mingamma2pt']
options.minmjj = analysis['mjj']
options.minHT = analysis['HT']
options.minggmindr = analysis['ggmindr']
options.maxjeteta = analysis['maxjeteta']
options.maxVBFjeteta = analysis['maxVBFjeteta']
options.maxgammaeta = analysis['maxgammaeta']
smear = analysis['s']

identifier = options.identifier
if not smear: identifier+='_nosmear'

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
  jetids = []

  VBF_jetpts = []
  VBF_jetetas = []
  VBF_jetphis = []
  VBF_jetms = []

  gammapts = []
  gammaetas = []
  gammaphis = []
  gammams = []

  j1pts = []
  j2pts = []

  jj1ms = []

  g1pts = []
  g2pts = []

  hpts = []
  for jentry in xrange(nentries):
    if jentry>options.numEvents and options.numEvents>0: continue
    cutflow[0]+=1 #entries read in -> cross section
    tree.GetEntry(jentry)

    if not jentry%1000:
      stdout.write('== \r%d events read ==\n'%jentry)
      stdout.flush() 

    Ntruthphotons = getattr(tree,options.Ntruthphotons)
    if not Ntruthphotons == 2: continue #only events with 2 truth photons

    hpt = 0

    treejpts = getattr(tree,options.jetpt)
    treejetas = getattr(tree,options.jeteta)
    treejphis = getattr(tree,options.jetphi)
    treejms = getattr(tree,options.jetm)
    treejids = getattr(tree,options.jetid)
    
    jpts = array([])
    jetas = array([])
    jphis = array([])
    jms = array([])
    jids = array([])

    for jpt,jeta,jphi,jm,jid in zip(treejpts,treejetas,treejphis,treejms,treejids):
      #if jid==22: continue #don't count jets formed from photons
      jpts = append(jpts,jpt)
      jetas = append(jetas,jeta)
      jphis = append(jphis,jphi)
      jms = append(jms,jm)
      jids = append(jids,jid)

    if smear:
      smeared_jpts = jet_smear(jpts,jetas,jms)
      #jpts = smeared_jpts
      if len(jpts)>0: jpts,jetas,jphis,jms,jids = array(zip(*sorted(zip(smeared_jpts,jetas,jphis,jms,jids),reverse=True)))

    good_indices = jpts>options.minjetpt
    jpts = jpts[good_indices]
    jetas = jetas[good_indices]
    jphis = jphis[good_indices]
    jms = jms[good_indices]
    jids = jids[good_indices]

    '''jjms = []
    if len(jpts)>=2:
      for ij1 in range(len(jpts)):
        if jpts[ij1]<options.minjetpt: continue
        for ij2 in range(ij1+1,len(jpts)):
          if jpts[ij2]<options.minjetpt: continue
          vjj = combine_vecs([jpts[ij1],jpts[ij2]],[jetas[ij1],jetas[ij2]],[jphis[ij1],jphis[ij2]],[jms[ij1],jms[ij2]])
          jjms.append(vjj.M())
    if len(jjms)>0: jj1ms.append(max(jjms))
    else: jj1ms.append(-1)'''

    gpts = array([g for g in getattr(tree,options.gammapt)])
    getas = array([g for g in getattr(tree,options.gammaeta)])
    gphis = array([g for g in getattr(tree,options.gammaphi)])
    gms = array([g for g in getattr(tree,options.gammam)])

    if smear:
      smeared_gpts = gamma_smear(gpts,getas)
      #gpts = smeared_gpts
      if len(gpts)>0: gpts,getas,gphis,gms = array(zip(*sorted(zip(smeared_gpts,getas,gphis,gms),reverse=True)))

    good_indices = gpts>options.mingammapt
    gpts = gpts[good_indices]
    getas = getas[good_indices]
    gphis = gphis[good_indices]
    gms = gms[good_indices]

    jetpts.append(jpts)
    jetetas.append(jetas)
    jetphis.append(jphis)
    jetms.append(jms)
    jetids.append(jids)

    gammapts.append(gpts)
    gammaetas.append(getas)
    gammaphis.append(gphis)
    gammams.append(gms)

  '''binwidth = 1
  n,bins,patches = plt.hist(all_jetmults,normed=True,bins=numpy.arange(0, max(all_jetmults) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,8)
  plt.ylim(0,max(n)*1.1)
  plt.xlabel('Jet multiplicity (>20 GeV, no $\eta$ cut)')
  plt.ylabel('a.u.')
  plotname = 'all_jetmult'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(all_jetmults)
  hists[plotname]['Nfinal'] = sum(array(all_jetmults))'''

  dump(jetpts,open(options.submitDir+'/'+'jetpts_'+identifier+'.p',"wb"))
  dump(jetetas,open(options.submitDir+'/'+'jetetas_'+identifier+'.p',"wb"))
  dump(jetphis,open(options.submitDir+'/'+'jetphis_'+identifier+'.p',"wb"))
  dump(jetms,open(options.submitDir+'/'+'jetms_'+identifier+'.p',"wb"))
  dump(jetids,open(options.submitDir+'/'+'jetids_'+identifier+'.p',"wb"))
  dump(gammapts,open(options.submitDir+'/'+'gammapts_'+identifier+'.p',"wb"))
  dump(gammaetas,open(options.submitDir+'/'+'gammaetas_'+identifier+'.p',"wb"))
  dump(gammaphis,open(options.submitDir+'/'+'gammaphis_'+identifier+'.p',"wb"))
  dump(gammams,open(options.submitDir+'/'+'gammams_'+identifier+'.p',"wb"))

  return [jetpts,jetetas,jetphis,jetms,jetids,gammapts,gammaetas,gammaphis,gammams]

def read():
  try:
    if options.root: raise IOError
    print '<< Trying to read in pickle data >>'
    jetpts = load(open(options.submitDir+'/'+'jetpts_'+identifier+'.p',"rb"))
    jetetas = load(open(options.submitDir+'/'+'jetetas_'+identifier+'.p',"rb"))
    jetphis = load(open(options.submitDir+'/'+'jetphis_'+identifier+'.p',"rb"))
    jetms = load(open(options.submitDir+'/'+'jetms_'+identifier+'.p',"rb"))
    jetids = load(open(options.submitDir+'/'+'jetids_'+identifier+'.p',"rb"))
    gammapts = load(open(options.submitDir+'/'+'gammapts_'+identifier+'.p',"rb"))
    gammaetas = load(open(options.submitDir+'/'+'gammaetas_'+identifier+'.p',"rb"))
    gammaphis = load(open(options.submitDir+'/'+'gammaphis_'+identifier+'.p',"rb"))
    gammams = load(open(options.submitDir+'/'+'gammams_'+identifier+'.p',"rb"))
    print '<< Pickle data read in>>'
    data = [jetpts,jetetas,jetphis,jetms,jetids,gammapts,gammaetas,gammaphis,gammams]
  except IOError:
    print '<< Reading in Root data >>'
    data = readRoot()

  cutflow[1]=len(data[0])
  return data

def trigger(data):
  [jetpts,jetetas,jetphis,jetms,jetids,gammapts,gammaetas,gammaphis,gammams] = data

  newjetpts = []
  newjetetas = []
  newjetphis = []
  newjetms = []

  VBF_jetpts = []
  VBF_jetetas = []
  VBF_jetphis = []
  VBF_jetms = []

  newgammapts = []
  newgammaetas = []
  newgammaphis = []
  newgammams = []

  j1pts = []
  j2pts = []
  g1pts = []
  g2pts = []
  HTs = []

  print '<< Applying triggers >>'
  for jpts,jetas,jphis,jms,jids,gpts,getas,gphis,gms in zip(jetpts,jetetas,jetphis,jetms,jetids,gammapts,gammaetas,gammaphis,gammams):
    passtrigger = True

    triggerjpts = jpts[abs(jetas)<options.maxVBFjeteta]
    jettriggers = [options.minjet1pt,options.minjet2pt,options.minjet3pt,options.minjet4pt]
    for jeti in range(4):
      if jettriggers[jeti]>0:
        if len(triggerjpts)<=jeti or triggerjpts[jeti]<jettriggers[jeti]: passtrigger=False
    if len(triggerjpts)==0: j1pts.append(0)
    else: j1pts.append(jpts[0])
    if len(triggerjpts)<2: j2pts.append(0)
    else: j2pts.append(jpts[1])

    HT = sum(triggerjpts)
    HTs.append(HT)
    if HT<options.minHT: passtrigger = False

    triggergpts = gpts[abs(getas)<options.maxgammaeta]
    gammatriggers = [options.mingamma1pt,options.mingamma2pt]
    for gammai in range(2):
      if gammatriggers[gammai]>0:
        if len(triggergpts)<=gammai or triggergpts[gammai]<gammatriggers[gammai]: passtrigger=False
    if len(triggergpts)==0: g1pts.append(0)
    else: g1pts.append(gpts[0])
    if len(triggergpts)<2: g2pts.append(0)
    else: g2pts.append(gpts[1])


    if not passtrigger: continue

    newjpts = jpts[all([abs(jetas)<options.maxjeteta,jids!=22],axis=0)]
    newjetas = jetas[all([abs(jetas)<options.maxjeteta,jids!=22],axis=0)]
    newjphis = jphis[all([abs(jetas)<options.maxjeteta,jids!=22],axis=0)]
    newjms = jms[all([abs(jetas)<options.maxjeteta,jids!=22],axis=0)]
    newjetpts.append(newjpts)
    newjetetas.append(newjetas)
    newjetphis.append(newjphis)
    newjetms.append(newjms)

    VBF_jpts = jpts[all([abs(jetas)<options.maxVBFjeteta,jids!=22],axis=0)]
    VBF_jetas = jetas[all([abs(jetas)<options.maxVBFjeteta,jids!=22],axis=0)]
    VBF_jphis = jphis[all([abs(jetas)<options.maxVBFjeteta,jids!=22],axis=0)]
    VBF_jms = jms[all([abs(jetas)<options.maxVBFjeteta,jids!=22],axis=0)]
    VBF_jetpts.append(VBF_jpts)
    VBF_jetetas.append(VBF_jetas)
    VBF_jetphis.append(VBF_jphis)
    VBF_jetms.append(VBF_jms)

    newgpts = gpts[abs(getas)<options.maxgammaeta]
    newgetas = getas[abs(getas)<options.maxgammaeta]
    newgphis = gphis[abs(getas)<options.maxgammaeta]
    newgms = gms[abs(getas)<options.maxgammaeta]
    newgammapts.append(newgpts)
    newgammaetas.append(newgetas)
    newgammaphis.append(newgphis)
    newgammams.append(newgms)

  print '<< Trigger efficiency: '+str(float(len(newjetpts))/len(jetpts))+'  >>'

  binwidth = 5
  n,bins,patches = plt.hist(g1pts,normed=True,bins=numpy.arange(0, max(g1pts) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,200)
  plt.ylim(0,max(n)*1.1)
  plt.xlabel('Leading $\gamma$ $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'gam1pt'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(g1pts)
  hists[plotname]['Nfinal'] = sum(array(g1pts)>options.mingamma1pt)

  n,bins,patches = plt.hist(g2pts,normed=True,bins=numpy.arange(0, max(g2pts) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,100)
  plt.xlabel('Subleading $\gamma$ $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'gam2pt'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(g2pts)
  hists[plotname]['Nfinal'] = sum(array(g2pts)>options.mingamma2pt)

  n,bins,patches = plt.hist(j1pts,normed=True,bins=numpy.arange(0, max(j1pts) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,200)
  plt.ylim(0,max(n)*1.1)
  plt.xlabel('Leading jet $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'jet1pt'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(j1pts)
  hists[plotname]['Nfinal'] = sum(array(j1pts)>options.minjet1pt)

  n,bins,patches = plt.hist(j2pts,normed=True,bins=numpy.arange(0, max(j2pts) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,100)
  plt.xlabel('Subleading jet $p_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'jet2pt'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(j2pts)
  hists[plotname]['Nfinal'] = sum(array(j2pts)>options.minjet2pt)

  binwidth = 10
  n,bins,patches = plt.hist(HTs,normed=True,bins=numpy.arange(0, max(HTs) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,500)
  plt.xlabel('$H_T$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'HT'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(HTs)
  hists[plotname]['Nfinal'] = sum(array(HTs)>options.minHT)


  return newjetpts,newjetetas,newjetphis,newjetms,VBF_jetpts,VBF_jetetas,VBF_jetphis,VBF_jetms,newgammapts,newgammaetas,newgammaphis,newgammams

def initial_cuts(data):
  print '<< Applying initial cuts >>'
  [jetpts,jetetas,jetphis,jetms,VBF_jetpts,VBF_jetetas,VBF_jetphis,VBF_jetms,gammapts,gammaetas,gammaphis,gammams] = data
  global cutflow,hists
  VBF_jetmults = [len(j) for j in VBF_jetpts]
  jetmults = [len(j) for j in jetpts]
  gammamults = [len(g) for g in gammapts]

  n,bins,patches = plt.hist(jetmults,normed=True,bins=range(0,max(jetmults)+1,1),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Jet Multiplicity ($|\eta|<2.5$)')
  plt.ylabel('a.u.')
  plotname = 'jetmultiplicity'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(jetmults)
  hists[plotname]['Nfinal'] = sum(array(jetmults)>=2)

  n,bins,patches = plt.hist(VBF_jetmults,normed=True,bins=range(0,max(VBF_jetmults)+1,1),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Jet Multiplicity ($|\eta|<4.5$)')
  plt.ylabel('a.u.')
  plotname = 'fjetmultiplicity'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(VBF_jetmults)
  hists[plotname]['Nfinal'] = sum(array(VBF_jetmults)>=0)

  n,bins,patches = plt.hist(gammamults,normed=True,bins=range(0,max(gammamults)+1,1),facecolor='b',histtype='stepfilled')
  #note: last bin is max+max-1
  plt.xlabel('Photon Multiplicity')
  plt.ylabel('a.u.')
  plotname = 'gammamultiplicity'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(gammamults)
  hists[plotname]['Nfinal'] = sum(array(gammamults)>=2)

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

  gammagammams = []
  VBF_jj1ms = []

  for i,(jpt,jeta,jphi,jm,jmult,VBF_jpt,VBF_jeta,VBF_jphi,VBF_jm,gpt,geta,gphi,gm,gmult) in enumerate(zip(jetpts,jetetas,jetphis,jetms,jetmults,VBF_jetpts,VBF_jetetas,VBF_jetphis,VBF_jetms,gammapts,gammaetas,gammaphis,gammams,gammamults)):
    if gmult<2:
      continue
    if gpt[0]<options.mingamma1pt or gpt[1]<options.mingamma2pt: continue
    cutflow[2]+=1
    gammagammap4 = combine_vecs([gpt[0],gpt[1]],[geta[0],geta[1]],[gphi[0],gphi[1]],[gm[0],gm[1]])
    gammagammams.append(gammagammap4.M())

    if jmult<2:
      continue
    if jpt[0]<options.minjet1pt or jpt[1]<options.minjet2pt: continue
    cutflow[3]+=1
    jjmindr = mindr(jeta,jphi,jeta,jphi,True)
    jjmindrs.append(jjmindr)
    ggmindr = mindr(geta,gphi,geta,gphi,True)
    ggmindrs.append(ggmindr)
    gjmindr = mindr(geta,gphi,jeta,jphi,False)
    gjmindrs.append(gjmindr)
    if jjmindr>0.4 and ggmindr>options.minggmindr and gjmindr>0.4:
      newjetpts.append(array([jpt[0],jpt[1]]))
      newjetetas.append(array([jeta[0],jeta[1]]))
      newjetphis.append(array([jphi[0],jphi[1]]))
      newjetms.append(array([jm[0],jm[1]]))

      newgammapts.append(array([gpt[0],gpt[1]]))
      newgammaetas.append(array([geta[0],geta[1]]))
      newgammaphis.append(array([gphi[0],gphi[1]]))
      newgammams.append(array([gm[0],gm[1]]))

      '''jjms = []'''
      VBF_jjms = inv_masses(VBF_jpt,VBF_jeta,VBF_jphi,VBF_jm)
      VBF_jj1ms.append(max(VBF_jjms))

  print '<< Initial cuts efficiency: '+str(float(len(newjetpts))/len(jetpts))+' >>'
  cutflow[4] = len(newjetpts)
  
  binwidth = 5
  n,bins,patches = plt.hist(gammagammams,normed=True,bins=numpy.arange(0, max(gammagammams) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,110)
  plt.ylim(0,max(n)*1.1)
  plt.xlabel('$m_{\gamma\gamma}$ (GeV) [No cuts]')
  plt.ylabel('a.u.')
  plotname = 'mgamgam_nocuts'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(gammagammams)
  hists[plotname]['Nfinal'] = len(gammagammams)

  binwidth = 0.15
  n,bins,patches = plt.hist(jjmindrs,normed=True,bins=numpy.arange(0, max(jjmindrs) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum jet-jet $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'jjmindr'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(jjmindrs)
  hists[plotname]['Nfinal'] = sum(array(jjmindrs)>0.4)

  n,bins,patches = plt.hist(ggmindrs,normed=True,bins=numpy.arange(0, max(ggmindrs) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum $\gamma-\gamma$ $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'ggmindr'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(ggmindrs)
  hists[plotname]['Nfinal'] = sum(array(ggmindrs)>options.minggmindr)

  n,bins,patches = plt.hist(gjmindrs,normed=True,bins=numpy.arange(0, max(gjmindrs) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.0)
  plt.xlabel('Minimum jet-$\gamma$ $\Delta R$')
  plt.ylabel('a.u.')
  plotname = 'gjmindr'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(gjmindrs)
  hists[plotname]['Nfinal'] = sum(array(gjmindrs)>0.4)

  binwidth = 20
  n,bins,patches = plt.hist(VBF_jj1ms,normed=True,bins=numpy.arange(0, max(VBF_jj1ms) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,1000)
  plt.xlabel('Max $m_{jj}$ (GeV)')
  plt.ylabel('a.u.')
  plotname = 'mjj_max'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = sum(array(VBF_jj1ms)>0)
  hists[plotname]['Nfinal'] = sum(array(VBF_jj1ms)>options.minmjj)

  #data = namedtuple('data',['jpts','jetas','jphis','jms','gpts','getas','gphis','gms',
  return [array(newjetpts),array(newjetetas),array(newjetphis),array(newjetms),array(newgammapts),array(newgammaetas),array(newgammaphis),array(newgammams),array(VBF_jj1ms)]

def final_cuts(data):
  print '<< Applying final cuts >>'
  [jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams,VBF_jj1ms] = data
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
  n,bins,patches = plt.hist(jetamasses,normed=True,bins=numpy.arange(0, max(jetamasses) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,60)
  plt.xlabel(r'$M_{jj}$')
  plt.ylabel('a.u.')
  plotname = 'mjj'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(jetamasses)
  hists[plotname]['Nfinal'] = len(jetamasses)

  binwidth = 0.2
  n,bins,patches = plt.hist(jetdphis,normed=True,bins=numpy.arange(0, max(jetdphis) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.5)
  plt.xlabel(r'$\Delta\phi_{jj}$')
  plt.ylabel('a.u.')
  plotname = 'dphijj'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(jetdphis)
  hists[plotname]['Nfinal'] = sum(array(jetdphis)<1.0)

  gammadphis = []
  gammaamasses = []
  gammaavecs = []
  for ept,eeta,ephi,em in zip(gammapts,gammaetas,gammaphis,gammams):
    gammadphis.append(dphi(ephi[0],ephi[1]))
    va = combine_vecs(ept,eeta,ephi,em)
    gammaavecs.append(va)
    gammaamasses.append(va.M())

  binwidth = 5
  n,bins,patches = plt.hist(gammaamasses,normed=True,bins=numpy.arange(0, max(gammaamasses) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,110)
  plt.xlabel(r'$M_{\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'mgamgam'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(gammaamasses)
  hists[plotname]['Nfinal'] = sum(array(gammaamasses)>-1)

  binwidth = 0.2
  n,bins,patches = plt.hist(gammadphis,normed=True,bins=numpy.arange(0, max(gammadphis) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.5)
  plt.xlabel(r'$\Delta\phi_{\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'dphigamgam'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(gammadphis)
  hists[plotname]['Nfinal'] = sum(array(gammadphis)<1.3)

  jetdphis = array(jetdphis)
  gammadphis = array(gammadphis)
  jetamasses = array(jetamasses)
  gammaamasses = array(gammaamasses)

  data = abs(jetamasses-gammaamasses)
  binwidth = 5
  n,bins,patches = plt.hist(data,normed=True,bins=numpy.arange(0, max(data) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,50)
  plt.xlabel(r'$|M_{jj}-M_{\gamma\gamma}|$')
  plt.ylabel('a.u.')
  plotname = 'amassdiff'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(data)
  hists[plotname]['Nfinal'] = sum(array(data)<15)

  good_indices = abs(jetamasses-gammaamasses)<15
  binwidth = 0.2
  n,bins,patches = plt.hist(jetdphis[good_indices],normed=True,bins=numpy.arange(0, max(jetdphis[good_indices]) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.xlim(0,3.5)
  plt.xlabel(r'$\Delta\phi_{jj}\cdot|M_{jj}-M_{\gamma\gamma}|<15.0$ GeV')
  plt.ylabel('a.u.')
  plotname = 'dphijj_amassdiff'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(jetdphis[good_indices])
  hists[plotname]['Nfinal'] = sum(array(jetdphis[good_indices])<1.0)

  binwidth = 20
  n,bins,patches = plt.hist(VBF_jj1ms[good_indices],normed=True,bins=numpy.arange(0, max(VBF_jj1ms[good_indices]) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
  plt.ylim(0,max(n)*1.1)
  plt.xlim(0,1000)
  plt.xlabel('Max $M_{jj}$ (GeV) $\cdot|M_{jj}-M_{\gamma\gamma}|<15.0$ GeV')
  plt.ylabel('a.u.')
  plotname = 'mjj_max_amassdiff'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = sum(array(VBF_jj1ms[good_indices])>0)
  hists[plotname]['Nfinal'] = sum(array(VBF_jj1ms[good_indices])>options.minmjj)

  good_indices = all([jetdphis<1.0,gammadphis<1.3,abs(jetamasses-gammaamasses)<15,VBF_jj1ms>options.minmjj],axis=0)
  cutflow[5] = sum(good_indices)

  hmasses = array([(va1+va2).M() for va1,va2 in zip(gammaavecs,jetavecs)])
  binwidth = 20
  if sum(good_indices)>0:
    n,bins,patches = plt.hist(hmasses[good_indices],normed=True,bins=numpy.arange(0, max(hmasses[good_indices]) + binwidth, binwidth),facecolor='b',histtype='stepfilled')
    plt.ylim(0,max(n)*1.1)
  else:
    n,bins,patches = plt.hist(hmasses[good_indices],normed=True,bins=numpy.arange(0, 100 + binwidth, binwidth),facecolor='b',histtype='stepfilled')
    plt.ylim(0,1)
  plt.xlim(0,200)
  plt.xlabel(r'$M_{jj\gamma\gamma}$')
  plt.ylabel('a.u.')
  plotname = 'mh'
  plt.savefig(options.plotDir+'/'+analysis_label+'_'+plotname+'_'+identifier+'.png')
  plt.close()
  hists[plotname] = {}
  hists[plotname]['hists'] = (n/sum(n),bins)
  hists[plotname]['Ninitial'] = len(hmasses[good_indices])
  hists[plotname]['Nfinal'] = sum(abs(hmasses[good_indices]-HIGGS_MASS)<25)

  good_indices = all([jetdphis<1.0,gammadphis<1.3,abs(jetamasses-gammaamasses)<15,VBF_jj1ms>options.minmjj,abs(hmasses-HIGGS_MASS)<25],axis=0)
  cutflow[6] = sum(good_indices)

  return

hists = {}
cutflow = [0,0,0,0,0,0,0]
#jetpts,jetetas,jetphis,jetms,gammapts,gammaetas,gammaphis,gammams = readRoot()
#ijetpts,ijetetas,ijetphis,ijetms,igammapts,igammaetas,igammaphis,igammams,iVBF_jj1ms = initial_cuts()
#final_cuts(ijetpts,ijetetas,ijetphis,ijetms,igammapts,igammaetas,igammaphis,igammams,iVBF_jj1ms)
final_cuts(initial_cuts(trigger(read())))
hists['Ninitial'] = cutflow[1]
hists['Nfinal'] = cutflow[6]
with open(options.submitDir+'/cutflow_'+analysis_label+'_'+identifier+'.p','wb') as outfile:
  dump(cutflow,outfile)
with open(options.submitDir+'/hists_'+analysis_label+'_'+identifier+'.p','wb') as outfile:
  dump(hists,outfile)
