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
parser.add_option("--lumi", help="Luminosity (fb^-1)",type=float, default=300)
#parser.add_option("-s","--smear", help="apply smearing to object pTs",action="store_true", default=False)
#parser.add_option("-v","--VBF", help="VBF analysis",action="store_true", default=False)
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

#for analysis in ['VBF0','VBF1','VBF2','VBF3']:
#for analysis in ['VBF0','VBF1','VBF4','VBF5']:
for analysis in ['VBF0','VBF1','GGH0','GGH1']:
  finalname = '_'+analysis+finalname
  #labels = json.load(open(options.submitDir+'/labels.json'))

  BRs = []
  ms = [10,20,30,40,50,60]
  for m in ms:
    Nsig = 0
    Nbkg = 0
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

    #missing_backgrounds = [1.2,2] # missing backgrounds make it 20-100% more
    missing_backgrounds = 1.2
    Nsig_norm = Nsig*options.lumi
    Nbkg_norm = Nbkg*options.lumi*missing_backgrounds #account for missing backgrounds
    bkgerr = Nbkg_norm*sqrt(0.1**2+1/float(Nbkg_norm))
    if Nsig_norm == 0: BR = 10 
    else: BR = maximum(5.0*bkgerr/Nsig_norm,10.0/Nsig_norm) #require at least 10 signal events
    BRs.append(BR)

  if analysis == 'VBF0': analysis_label='100\% trigger'
  if analysis == 'VBF1': analysis_label='g35g25'
  if analysis == 'VBF2': analysis_label='g25+4j35'
  if analysis == 'VBF3': analysis_label='g25+2j35+HT300'
  if analysis == 'VBF4': analysis_label='100\% trigger, no photon isolation'
  if analysis == 'VBF5': analysis_label='g35g25, no photon isolation'
  if analysis == 'GGH0': analysis_label='100\% trigger, no $m_{jj}$ cut'
  if analysis == 'GGH1': analysis_label='g35g25, no $m_{jj}$ cut'
  plt.semilogy(ms,BRs,label = analysis_label)

plt.semilogy([0],[1],ls = ' ',label='Luminosity: '+str(int(options.lumi))+'ifb')
plt.ylim(0.01,1)
plt.xlim(0,max(ms))
plt.xlabel('$m_a$ (GeV)')
plt.ylabel('BR($H\\rightarrow \gamma\gamma jj$) Required for Significance$\ge$5.0')
plt.legend(loc='upper right')
plt.savefig(options.plotDir+'/BR_ma'+'_lumi'+str(int(options.lumi))+finalname+'.png')
plt.close()

