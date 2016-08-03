from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram
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
parser.add_option("-s","--nosmear", help="don't apply smearing to object pTs",action="store_true", default=False)
parser.add_option("--analysis_label", help="Name of analysis", type=str, default="VBF0")
parser.add_option("--bkgcolor", help="color of background lines", default='r')
parser.add_option("--sigcolor", help="color of signal lines", default='b')

(options, args) = parser.parse_args()

analysis_label = options.analysis_label
with open(options.submitDir+'/'+analysis_label+'.json') as analysis_file:
  analysis = json.load(analysis_file)

hists = {}
keys = []
finalname = ''
try: identifiers = json.load(open(options.submitDir+'/'+options.identifiers))
except IOError: 
  try: identifiers = json.load(open(options.identifiers)) 
  except IOError:
    raise IOError('Can\'t find '+options.identifiers)

for identifier in identifiers:
  name = identifier['identifier']
  finalname+='_'+name
  if options.nosmear: name+='_nosmear'
  name = analysis_label+'_'+name
  #print name
  hists[name] = pickle.load(open(options.submitDir+'/hists_'+name+'.p','rb'))
  if len(keys)==0: keys = hists[name].keys()
  else:
    if not hists[name].keys() == keys: raise RuntimeError('Keys in hists don\'t match.')
if options.nosmear: finalname+='_nosmear'
finalname = '_'+options.analysis_label+finalname
labels = json.load(open(options.submitDir+'/labels.json'))

round_to_n = lambda x, n: round(x, -int(floor(log10(abs(x)))) + (n - 1))

for key in keys:
  if key=='Ninitial' or key=='Nfinal': continue
  loc = [0.87 for p in range(2)]
  ymaxs = [float('-inf'),float('-inf')]
  ymins = [float('inf'),float('inf')]
  setlogs = [False,False]
  Nsofar = [{},{}] 
  for i,identifier in enumerate(identifiers):
    name = identifier['identifier']
    if options.nosmear: name+='_nosmear'
    name = analysis_label+'_'+name
    if 'ls' in identifier: ls = identifier['ls']
    else: ls = '-'
    label = identifier['name']
    sig = identifier['sig']
    #color = identifier['color']
    color = options.sigcolor if sig else options.bkgcolor
    n = hists[name][key]['hists'][0]
    bins = hists[name][key]['hists'][1]
    n=numpy.insert(n,0,0)
    Ntotal = float(hists[name]['Ninitial'])
    Nin = float(hists[name][key]['Ninitial'])
    Nin_err = sqrt(Nin)
    Nfin = float(hists[name][key]['Nfinal'])
    Nfin_err = sqrt(Nfin)
    xsec = identifier['xsec']
    if 'GGH0_ggH_a' in name: xsec/=0.3295 # effective cross section of old sample
    norm = options.lumi*xsec*identifier['filter']*identifier['BR']
    for ib,b in enumerate(bins):
      if b not in Nsofar[sig]: Nsofar[sig][b] = n[ib]*Nin*norm/Ntotal
      else: Nsofar[sig][b] += n[ib]*Nin*norm/Ntotal
    #plt.bar(bins[0:len(bins)-1],n,width=bins[1]-bins[0],fill=False,color=color,label=label)
    if 'xlim' in labels[key]:
      minxlim = labels[key]['xlim'][0]
      maxxlim = labels[key]['xlim'][1]
    else:
      minxlim = min(bins)
      maxxlim = max(bins)
    for p in [0,1]:
      plt.figure(p)
      if p==0:
        if 'mjj_max' or 'HT' in key:
          plt.semilogy(bins,n,drawstyle='steps',ls=ls,color=color,label=label)
          setlogs[0]=True
        else:
          plt.plot(bins,n,drawstyle='steps',ls=ls,color=color,label=label)
        n = array(n)
        nplotted = n[all([bins<=maxxlim,bins>=minxlim,n>0],axis=0)]
        ymaxs[0] = max(max(nplotted),ymaxs[0])
        ymins[0] = min(min(nplotted),ymins[0])
      if p==1:
        plt.plot(bins,n*Nin*norm/Ntotal,drawstyle='steps',ls=ls,color=color,label=label)
        nplotted = n[all([bins<=maxxlim,bins>=minxlim,n>0],axis=0)]*Nin*norm/Ntotal
        ymaxs[1] = max(max(nplotted),ymaxs[1])
        ymins[1] = min(min(nplotted),ymins[1])
      plt.figtext(0.71,loc[p],label+':')
      duplicated_events = 1.
      if Nin>0:
        plt.figtext(0.75,loc[p]-0.05,'N before cuts: ' + str(round_to_n(Nin*norm/Ntotal,3))+'$\pm$'+str(round_to_n(sqrt(Nin)*sqrt(duplicated_events)*norm/Ntotal,3)))
      else:
        plt.figtext(0.75,loc[p]-0.05,'N before cuts: 0')
      if Nfin>0:
        plt.figtext(0.75,loc[p]-0.1,'N after cuts: ' + str(round_to_n(Nfin*norm/Ntotal,3))+'$\pm$'+str(round_to_n(sqrt(Nfin)*sqrt(duplicated_events)*norm/Ntotal,3)))
      else:
        plt.figtext(0.75,loc[p]-0.1,'N after cuts: 0')
      loc[p] = loc[p]-0.15
      if identifier['BR']<1:
        plt.figtext(0.75,loc[p],'Assumed BR: ' + str(round_to_n(identifier['BR'],3)))
        loc[p] = loc[p]-0.05

  for sig in [0,1]:
    plt.figure(1)
    xs = Nsofar[sig].keys()
    ys = [Nsofar[sig][xx] for xx in xs]
    xs, ys = zip(*sorted(zip(xs, ys)))
    label = 'Total Signal' if sig else 'Total Background'
    color = options.sigcolor if sig else options.bkgcolor
    plt.plot(xs,ys,drawstyle='steps',ls='-',color=color,label=label)
    ys = array(ys)
    xs = array(xs)
    ysplotted = ys[all([xs<=maxxlim,xs>=minxlim],axis=0)]
    ymaxs[1] = max(max(ysplotted),ymaxs[1])
    ymins[1] = min(min(ysplotted),ymins[1])

  for p in range(2):
    plt.figure(p)
    plt.xlim(minxlim,maxxlim)
    ax = plt.gca()
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.74, box.height])

  for p in range(2):
    plt.figure(p)
    axes = plt.gca()
    plt.figtext(0.71,loc[p],'Luminosity: ' +str(options.lumi)+' ifb')
    plt.figtext(0.71,loc[p]-0.05,'Mass units: GeV')
    #if p==1: 
      #plt.figtext(0.71,loc[p]-0.1,'Backgrounds, signals added')

  for i in range(2):
    plt.figure(i)
    if setlogs[i]:
      plt.ylim(ymins[i],ymaxs[i]*10)
    else:
      plt.ylim(0,ymaxs[i]*1.4)
    binsize = bins[1]-bins[0]
    plt.xlabel(labels[key]['xlabel'])
    cuts = []
    greaters = []
    if 'cuts' in labels[key]:
      cuts = [c['cut'] for c in labels[key]['cuts']]
      greaters = [c['greater'] for c in labels[key]['cuts']]
    if 'cut' in labels[key]:
      cut = labels[key]['cut']
      if key=='jet1pt': cut = max(analysis['minjetpt'],analysis['minjet1pt'])
      if key=='jet2pt': cut = max(analysis['minjetpt'],analysis['minjet2pt'])
      if key=='gam1pt': cut = max(analysis['mingammapt'],analysis['mingamma1pt'])
      if key=='gam2pt': cut = max(analysis['mingammapt'],analysis['mingamma2pt'])
      if key=='mjj_max': cut = analysis['mjj']
      if key=='HT': cut = analysis['HT']
      if key=='jetmultiplicity': cut = min(sum([analysis['minjet1pt']>0,analysis['minjet2pt']>0,analysis['minjet3pt']>0,analysis['minjet3pt']>0]),2)
      if key=='ggmindr': cut = analysis['ggmindr']
      #if key=='ggmindr'
      cuts.append(cut)
      greaters.append(labels[key]['greater'])
    for cut,greater in zip(cuts,greaters):
      if not greater<0:
        axes = plt.gca()
        ymin,ymax = axes.get_ylim()
        xmin,xmax = axes.get_xlim()
        plt.plot([cut,cut],[ymin,ymax],color='g',linestyle='--')
        for x in numpy.arange(0,0.51,0.05):
          if greater: x = -x
          plt.fill_between([cut+x*binsize,cut+(x+0.05)*binsize],[ymin,ymin],[ymax,ymax],alpha=0.3*(1-abs(x)/0.5),color='g',linewidth=0.0)
        plt.ylim(ymin,ymax)
        plt.xlim(xmin,xmax)
    plt.legend(loc='upper right')
    if i==0:
      plt.ylabel('Fraction of Events')
      plt.savefig(options.plotDir+'/'+key+'_compare'+finalname)
    if i==1:
      plt.ylabel('Number of Events')
      plt.savefig(options.plotDir+'/'+key+'_normalized_compare'+finalname)
    plt.close()


  


