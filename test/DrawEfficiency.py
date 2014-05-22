#!usr/bin/python
import sys, math
import ROOT
from ROOT import TFile, TH1F, TCanvas, TLegend, TGraphErrors, gROOT, gPad
from array import array
from setTDRStyle import setTDRStyle

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

mass = [1000,1500,2000,2500,3000]

## ---- CERN -------
PATH = 'root://eoscms//eos/cms/store/cmst3/group/das2014/EXODijetsLE/'
## ---- FNAL -------
# PATH = '/eos/uscms/store/user/cmsdas/2014/EXODijetsLE/'

y1  = []
y2  = []
x   = []
ey1 = []
ey2 = []
ex  = []

for im in mass:
  inf = TFile.Open(PATH+'dijetHisto_RS'+str(im)+'_signal.root')
  h   = inf.Get('TriggerPass')
  h1  = inf.Get('h_mjj')
  h2  = inf.Get('h_sub_mjj')
  N   = h.GetBinContent(1)
  N1  = h1.GetEntries()
  N2  = h2.GetEntries()
  p1  = N1/N
  p2  = N2/N
  y1.append(p1)
  y2.append(p2)
  x.append(im*0.001)
  ey1.append(math.sqrt(p1*(1-p1)/N))
  ey2.append(math.sqrt(p2*(1-p2)/N))
  ex.append(0.0)
  
k = 0
print 'Mass (TeV)    simple     substructure'
while k < len(mass):
  print '  '+str(mass[k]*0.001)+'          '+str(round(y1[k],3))+'       '+str(round(y2[k],3))
  k+=1
vy1  = array('d',y1)
vy2  = array('d',y2) 
vx   = array('d',x)
vey1 = array('d',ey1)
vey2 = array('d',ey2) 
vex  = array('d',ex)

g1 = TGraphErrors(len(vx),vx,vy1,vex,vey1)
g2 = TGraphErrors(len(vx),vx,vy2,vex,vey2)

g1.SetLineColor(ROOT.kBlue)
g2.SetLineColor(ROOT.kRed)
g1.SetMarkerColor(ROOT.kBlue)
g2.SetMarkerColor(ROOT.kRed)
g1.SetMarkerStyle(20)
g2.SetMarkerStyle(21)
g1.SetLineWidth(2)
g2.SetLineWidth(2)

can = TCanvas('SignalEfficiency','SignalEfficiency',900,600)
g1.GetXaxis().SetTitle('Resonance Mass (TeV)')
g1.GetYaxis().SetTitle('Efficiency')
g1.GetYaxis().SetNdivisions(505)
g1.GetYaxis().SetRangeUser(0,0.5)
g1.Draw('ALPE')
g2.Draw('sameLPE')

leg = TLegend(0.2,0.75,0.5,0.9)
leg.AddEntry(g1,'Simple selection','LP')
leg.AddEntry(g2,'Substructure selection','LP')
leg.SetFillColor(0)
leg.SetTextFont(42)
leg.SetTextSize(0.04)
leg.Draw()

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]











