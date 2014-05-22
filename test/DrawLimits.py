#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TTree, TH1F, TCanvas, TLegend, TGraph, TGraphAsymmErrors, gROOT, gPad
from array import array
from setTDRStyle import setTDRStyle
import optparse

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--useSub",action="store_true",default=False,dest="useSub")
(options, args) = parser.parse_args()useSub = options.useSub

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

mass = [1000,1500,2000,2500,3000]
xsection = [4.254e-1,3.298e-2,4.083e-3,6.191e-4,1.010e-4]
x = []
exl = []
exh = []
y1 = []
y2 = []
yExp = []
yObs = []
ey1l = []
ey1h = []
ey2l = []
ey2h = []  

k = 0

for im in mass:
  filename = 'higgsCombineNormal.Asymptotic.mH'+str(im)+'.root' 
  if useSub:
    filename = 'higgsCombineSub.Asymptotic.mH'+str(im)+'.root' 
  print filename
  inf = TFile.Open(filename)
  tr = inf.Get('limit')
  y = []
  N = tr.GetEntriesFast()
  for i in xrange(N): 
    tr.GetEntry(i)
    y.append(tr.limit)
  
  x.append(im*0.001)
  exl.append(0.0)
  exh.append(0.0)
  y1.append(y[2]*xsection[k])
  y2.append(y[2]*xsection[k])
  yExp.append(y[2]*xsection[k])
  yObs.append(y[5]*xsection[k])

  ey1l.append((y[2]-y[1])*xsection[k])
  ey2l.append((y[2]-y[0])*xsection[k])
  ey1h.append((y[3]-y[2])*xsection[k])
  ey2h.append((y[4]-y[2])*xsection[k])
  
  k+=1
  
vxs = array('d',xsection)
vx = array('d',x)
vyObs = array('d',yObs)
vyExp = array('d',yExp)
vy1 = array('d',y1)
vy2 = array('d',y2)
vexl = array('d',exl)
vexh = array('d',exh)
vey1l = array('d',ey1l)
vey1h = array('d',ey1h)
vey2l = array('d',ey2l)
vey2h = array('d',ey2h)

g1   = TGraphAsymmErrors(len(vx),vx,vy1,vexl,vexh,vey1l,vey1h)
g2   = TGraphAsymmErrors(len(vx),vx,vy2,vexl,vexh,vey2l,vey2h)
gObs = TGraph(len(vx),vx,vyObs)
gExp = TGraph(len(vx),vx,vyExp)
gxs  = TGraph(len(vx),vx,vxs)

g1.SetFillColor(ROOT.kYellow)
g1.SetLineColor(ROOT.kYellow)
g2.SetFillColor(ROOT.kGreen)
g2.SetLineColor(ROOT.kGreen)
gExp.SetLineWidth(2)
gExp.SetLineStyle(9)
gObs.SetLineWidth(2)
gObs.SetLineStyle(1)
gObs.SetLineColor(ROOT.kBlue+1)
gObs.SetMarkerColor(ROOT.kBlue+1)
gObs.SetMarkerStyle(21)
gObs.SetMarkerSize(1.5)
gxs.SetLineStyle(5)
gxs.SetLineWidth(3)
gxs.SetLineColor(ROOT.kRed)
    
canvasName = 'Limits'
if useSub:
  canvasName = 'Limits_Sub'
can = TCanvas(canvasName,canvasName,900,600)
gPad.SetLogy()
g2.GetXaxis().SetTitle('Resonance Mass (TeV)')
g2.GetYaxis().SetTitle('#sigma #times BR (X#rightarrow WW) (pb)')
g2.GetYaxis().CenterTitle(ROOT.kTRUE)
g2.GetYaxis().SetNdivisions(510)
g2.GetYaxis().SetRangeUser(1e-3,2)
g2.Draw('AE3')
g1.Draw('sameE3')
gExp.Draw('sameL')
gObs.Draw('sameLP')
gxs.Draw('sameL')

leg = TLegend(0.65,0.65,0.9,0.9)
if useSub:
  leg.SetHeader('Substructure Selection')
else:
  leg.SetHeader('Simle Dijet Selection')
leg.AddEntry(gObs,'Observed','LP')
leg.AddEntry(gExp,'Expected','L')
leg.AddEntry(g1,'Expected #pm 1 #sigma','F')
leg.AddEntry(g2,'Expected #pm 2 #sigma','F')
leg.AddEntry(gxs,'G_{RS}#rightarrow WW','L')
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
