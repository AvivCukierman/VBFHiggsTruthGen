from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram
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

(options, args) = parser.parse_args()

hists = {}
keys = []
finalname = ''
identifiers = json.load(open(options.submitDir+'/'+options.identifiers))
for identifier in identifiers:
  name = identifier['identifier']
  finalname+='_'+name
  hists[name] = pickle.load(open(options.submitDir+'/hists_'+name+'.p','rb'))
  if len(keys)==0: keys = hists[name].keys()
  else:
    if not hists[name].keys() == keys: raise RuntimeError('Keys in hists don\'t match.')

labels = json.load(open(options.submitDir+'/labels.json'))

for key in keys:
  for identifier in identifiers:
    name = identifier['identifier']
    color = identifier['color']
    label = identifier['name']
    n = hists[name][key][0]
    bins = hists[name][key][1]
    n=numpy.insert(n,0,0)
    #plt.bar(bins[0:len(bins)-1],n,width=bins[1]-bins[0],fill=False,color=color,label=label)
    plt.plot(bins,n,ls='steps',color=color,label=label)
    plt.xlim(min(bins),max(bins))

  binsize = bins[1]-bins[0]
  plt.xlabel(labels[key]['xlabel'])
  cut = labels[key]['cut']
  greater = labels[key]['greater']
  if not greater<0:
    axes = plt.gca()
    ymin,ymax = axes.get_ylim()
    plt.plot([cut,cut],[ymin,ymax],color='g',linestyle='--')
    for x in numpy.arange(0,0.51,0.05):
      if greater: x = -x
      plt.fill_between([cut+x*binsize,cut+(x+0.05)*binsize],[ymin,ymin],[ymax,ymax],alpha=0.3*(1-abs(x)/0.5),color='g',linewidth=0.0)
    plt.ylim(ymin,ymax)
  plt.legend(loc='upper right')
  plt.savefig(options.plotDir+'/'+key+'_compare'+finalname)
  plt.close()


  


