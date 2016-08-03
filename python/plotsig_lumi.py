from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram,maximum
from math import floor,log10
import pickle
import json
import os
import numpy
from optparse import OptionParser
from sys import stdout,argv
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
parser.add_option("--submitDir", help="Directory containing output files",type=str, default="output")
parser.add_option("--identifiers", help="File containing identifiers of data to compare",type=str, default="bkg_sig_test.json")
parser.add_option("--plotDir", help="Directory containing plots",type=str, default="plots")
#parser.add_option("--lumi", help="Luminosity (fb^-1)",type=float, default=300)
#parser.add_option("-s","--smear", help="apply smearing to object pTs",action="store_true", default=False)
#parser.add_option("-v","--VBF", help="VBF analysis",action="store_true", default=False)
parser.add_option("--VH", help="VH analysis",action="store_true", default=False)
#parser.add_option("--bkgcolor", help="color of background lines", default='r')
#parser.add_option("--sigcolor", help="color of signal lines", default='b')

(options, args) = parser.parse_args()

hists = {}
keys = []
finalname = ''
try: identifiers = json.load(open(options.submitDir+'/'+options.identifiers))
except IOError: 
  try: identifiers = json.load(open(options.identifiers)) 
  except IOError:
    raise IOError('Can\'t find '+options.identifiers)

for analysis,color in zip(['VBF0','GGH0'],['b','g']):
  finalname = '_'+analysis+finalname
  #labels = json.load(open(options.submitDir+'/labels.json'))

  Nsig = 0
  Nbkg = 0
  m=30
  for i,identifier in enumerate(identifiers):
    sig = identifier['sig']
    if sig:
      names = identifier['identifier']
      name = analysis+'_'+names[0]+str(m)+names[1]
      if analysis=='GGH0' or analysis=='GGH1':
        if 'ggH' in names[0]:
          namearr = str.split(str(name),'_')
          name = analysis
          for npart in namearr[1:]:
            name=name+'_'+npart
            if 'ggH' in npart: name=name+'_all'
          print name
    else:
      name = identifier['identifier']
      name = analysis+'_'+name

    print name
    hists[name] = pickle.load(open(options.submitDir+'/hists_'+name+'.p','rb'))
    if len(keys)==0: keys = hists[name].keys()
    else:
      if not hists[name].keys() == keys: raise RuntimeError('Keys in hists don\'t match.')

    Ntotal = float(hists[name]['Ninitial'])
    Nfinal = float(hists[name]['Nfinal'])
    xsec = identifier['xsec']
    if 'GGH_ggH_a' in name: xsec/=0.3295 # effective cross section of new sample
    if 'GGH0_ggH_a' in name: xsec/=0.3295 # effective cross section of new sample
    if 'GGH1_ggH_a' in name: xsec/=0.3295 # effective cross section of new sample
    norm = xsec*identifier['filter'] # lumi = 1 ifb, BR = 1
    N = Nfinal*norm/Ntotal
    if sig: Nsig+=N
    else: Nbkg+=N

  missing_backgrounds_range = [1.2,2.0] # missing backgrounds make it 20% more
  BRs = [[],[]]
  for mbkgi,missing_backgrounds in enumerate(missing_backgrounds_range):
    lumis = array(range(1,310))
    Nsig_arr = Nsig*lumis
    Nbkg_arr = Nbkg*lumis*missing_backgrounds #account for missing backgrounds
    bkgerrs = Nbkg_arr*sqrt(0.1**2+1./Nbkg_arr)
    BRs[mbkgi] = maximum(5.0*bkgerrs/Nsig_arr,10.0/Nsig_arr) #require at least 10 signal events

    if analysis == 'VBF0': analysis_label='VBF'
    if analysis == 'GGH0': analysis_label='GGH'
  plt.semilogy(lumis,BRs[0],label = analysis_label+' ($m_a = $'+str(m)+' GeV)')
  plt.fill_between(lumis,BRs[0],BRs[1],facecolor=color,alpha=0.15)

if options.VH:
  Nsig_arr = 50.*lumis/300/0.04
  Nbkg_arr = 100.*lumis/300
  bkgerrs = maximum(Nbkg_arr*0.1,sqrt(Nbkg_arr))
  VHBRs = maximum(5.0*bkgerrs/Nsig_arr,10.0/Nsig_arr) #require at least 10 signal events
  plt.semilogy(lumis,VHBRs,label = 'VH ($m_a = 30$ GeV)',ls='--')
  marker_style = dict(color='r', markersize=20, fillstyle='full')
  plt.semilogy(300,0.04,marker="*",**marker_style)
plt.ylim(0.01,1)
plt.xlim(0,max(lumis))
plt.xlabel('Luminosity (ifb)')
plt.ylabel('BR($H\\rightarrow \gamma\gamma jj$) Required for Significance$\ge$5.0')
plt.legend(loc='upper right')
plt.savefig(options.plotDir+'/BR_lumi'+finalname)
plt.close()
