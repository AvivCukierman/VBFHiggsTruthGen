from numpy import array,mean,std,all,sqrt,asarray,average,save,histogram,cosh
import numpy
import ROOT as r

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

def jet_smear(pts,etas,ms):
  Ps = pts*cosh(etas)
  Es = sqrt(Ps**2+ms**2)
  deltaEs = 0.8*sqrt(Es)
  smearedEs = numpy.random.normal(loc=Es,scale=deltaEs)
  smearedPs = numpy.where(smearedEs>ms,sqrt(smearedEs**2-ms**2),0)
  #if sum(smearedPs==0)/len(smearedPs) > 0.1: print '==A lot of jet pTs are being set to 0 in the smearing=='
  smearedpts = smearedPs/cosh(etas)
  return smearedpts

def gamma_smear(pts,etas):
  Es = pts*cosh(etas) #E=P
  deltaEs = 0.1*sqrt(Es)+0.007*Es
  smearedEs = numpy.random.normal(loc=Es,scale=deltaEs)
  smearedpts = smearedEs/cosh(etas)
  return smearedpts

def inv_masses(pts,etas,phis,ms):
  len_list = [len(pts),len(etas),len(phis),len(ms)]
  if len_list[1:] != len_list[:-1]: raise RuntimeError('pTs, etas, phis, and ms should all be same length')
  if len_list[0]<2: raise RuntimeError('Have to have at least two objects')

  inv_ms = []
  for i1 in range(len(pts)):
    for i2 in range(i1+1,len(pts)):
      v = combine_vecs([pts[i1],pts[i2]],[etas[i1],etas[i2]],[phis[i1],phis[i2]],[ms[i1],ms[i2]])
      inv_ms.append(v.M())
  return inv_ms
