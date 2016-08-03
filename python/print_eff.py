from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram,maximum
import csv
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

round_to_n = lambda x, n: round(x, -int(floor(log10(abs(x)))) + (n - 1))

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
with open('eggs.csv', 'wb') as csvfile:
  writer = csv.writer(csvfile, delimiter=' ')
  for analysis in ['VBF0','VBF1','VBF2','VBF3','GGH0','GGH1']:
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
              #print name
        else:
          name = identifier['identifier']
          name = analysis+'_'+name

        #print name
        hists[name] = pickle.load(open(options.submitDir+'/hists_'+name+'.p','rb'))
        if len(keys)==0: keys = hists[name].keys()
        else:
          if not hists[name].keys() == keys: raise RuntimeError('Keys in hists don\'t match.')

        Ntotal = float(hists[name]['Ninitial'])
        Nfinal = float(hists[name]['Nfinal'])
        if Nfinal>0:
          print name+'\t'+str(round_to_n(float(Nfinal)/Ntotal,3))
          writer.writerow([name,str(round_to_n(float(Nfinal)/Ntotal,3))])
        else:
          print name+'\t'+'0'
          writer.writerow([name,'0'])
