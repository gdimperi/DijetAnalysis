#!usr/bin/python
import ROOT
from ROOT import TFile, TTree, TH1, TH1F, TCanvas, TMath, gROOT
from math import *

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--sample",action="store",type="string",dest="sample",default='signal_M1000')
parser.add_option("--trigger",action="store",type="string",dest="trigger",default='signal')

(options, args) = parser.parse_args()
sample = options.sample
trigger = options.trigger

trigger_options = ['signal','ref','refSig']
if trigger not in trigger_options:
  print 'WARNING: the requested trigger option ('+trigger+') does not exist !!!'
  print 'Available trigger options are: '+str(trigger_options)
  print 'Forcing the default option: \"signal\"'
  trigger = 'signal'

## ---- CERN -------
PATH = '/cmshome/gdimperi//Dijet/CMSSW_7_1_0_pre9/src/CMSROMA/DijetAnalysis/prod/'


inputf  = PATH+'dijetTree_'+sample+'.root'
outputf = 'dijetHisto_'+sample+'_'+trigger+'.root'

gROOT.Reset()
#---- get the tree ---------------------
inf = TFile.Open(inputf)
events = inf.Get('dijets/events')

#---- create the output file -----------
outf = TFile(outputf,'RECREATE')

#---- get the trigger histogram --------
hTrig = inf.Get('dijets/TriggerPass')

#---- list of variables ----------------
eventVarName    = ['ht','mjj','dEtajj','dPhijj','metOvSumEt']
eventVarMin     = [0,0,0,0,0]
eventVarMax     = [3000,5000,5,pi,1]
eventVarBins    = [300,5000,100,100,100]
jetVarName      = ['jetTau21','jetPt','jetEta','jetPhi','jetMass','jetMassPruned']  
jetVarMin       = [0,0,-5,-pi,0,0]
jetVarMax       = [1,2000,5,pi,1000,1000]
jetVarBins      = [100,200,100,100,200,200]
histEventVar    = []
histJetVar      = []
histEventVarSub = []
histJetVarSub   = []
#---- create the histograms ------------
k = 0
for i in eventVarName:
  h = TH1F('h_'+i,'h_'+i,eventVarBins[k],eventVarMin[k],eventVarMax[k]) 
  h.Sumw2()
  histEventVar.append(h)
  hSub = TH1F('h_sub_'+i,'h_sub_'+i,eventVarBins[k],eventVarMin[k],eventVarMax[k]) 
  hSub.Sumw2()
  histEventVarSub.append(hSub)
  k += 1

k = 0
for i in jetVarName:
  h = TH1F('h_'+i,'h_'+i,jetVarBins[k],jetVarMin[k],jetVarMax[k]) 
  h.Sumw2()
  histJetVar.append(h)
  hSub = TH1F('h_sub_'+i,'h_sub_'+i,jetVarBins[k],jetVarMin[k],jetVarMax[k]) 
  hSub.Sumw2()
  histJetVarSub.append(hSub)
  k += 1  

#---- read the tree & fill histosgrams -
N = events.GetEntriesFast()
print 'Reading the input file: '+inf.GetName()
print 'Number of events: '+str(N)
d = 0

for i in xrange(N):
  events.GetEntry(i)
  #---- progress of the reading --------
  fraction = 10.*i/(1.*N)
  if TMath.FloorNint(fraction) > d:
    print str(10*TMath.FloorNint(fraction))+'%' 
  d = TMath.FloorNint(fraction)

  cut_trigger      = True 
  cut_mass         = True
  cut_dEtajj       = events.dEtajjCA8 < 1.3
  cut_leptonVeto   = events.jetElfCA8[0]<0.7 and events.jetElfCA8[1]<0.7 and events.jetMufCA8[0]<0.7 and events.jetMufCA8[1]<0.7
  cut_eta          = fabs(events.jetEtaCA8[0])<2.5 and fabs(events.jetEtaCA8[1])<2.5
  cut_pt           = events.jetPtCA8[1]>40
  cut_substructure = (
    events.jetMassPrunedCA8[0] > 60 and 
    events.jetMassPrunedCA8[0] < 100 and 
    events.jetMassPrunedCA8[1] > 60 and 
    events.jetMassPrunedCA8[1] < 100 and 
    events.jetTau2CA8[0]/events.jetTau1CA8[0] < 0.5 and 
    events.jetTau2CA8[1]/events.jetTau1CA8[1] < 0.5)

  if trigger == 'signal':
    cut_trigger = events.triggerResult[0] or events.triggerResult[1] or events.triggerResult[2] or events.triggerResult[3]
    cut_mass = events.mjjCA8 > 890
  elif trigger == 'ref':
    cut_trigger = events.triggerResult[4]
    cut_mass = events.mjjCA8 > 0
  elif trigger == 'refSig':
    cut_trigger = events.triggerResult[4] and (events.triggerResult[0] or events.triggerResult[1] or events.triggerResult[2] or events.triggerResult[3])
    cut_mass = events.mjjCA8 > 0
    
#################################################

# ---------------- test plots  (to be changed) ---------------------

###################################################
# just turn off all the cuts

  cut_trigger      = True 
  cut_mass         = True
  cut_dEtajj       = True
  cut_leptonVeto   = True
  cut_eta          = True
  cut_pt           = True
  cut_substructure = True
 

#--------------------------------------------------------------

  if (cut_trigger and cut_dEtajj and cut_mass and cut_leptonVeto and cut_eta and cut_pt):
    #---- set the event variables ----
    eventVar = [events.htCA8,events.mjjCA8,events.dEtajjCA8,events.dPhijjCA8,events.metSig]
    for k in xrange(len(histEventVar)):
      histEventVar[k].Fill(eventVar[k])
      if cut_substructure:
        histEventVarSub[k].Fill(eventVar[k])
    for j in xrange(2):
      #---- set the jet variables ----
      jetVar = [events.jetTau2CA8[j]/events.jetTau1CA8[j],events.jetPtCA8[j],events.jetEtaCA8[j],events.jetPhiCA8[j],events.jetMassCA8[j],events.jetMassPrunedCA8[j]]
      for k in xrange(len(histJetVar)):
        histJetVar[k].Fill(jetVar[k])
        if cut_substructure:
          histJetVarSub[k].Fill(jetVar[k])
      
#---- write output file ----------------
print 'Writing output file: '+outf.GetName()
outf.cd()
hTrig.Write()
outf.Write()
outf.Close()
inf.Close() 
  
