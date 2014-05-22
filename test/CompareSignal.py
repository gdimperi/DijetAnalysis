#!usr/bin/python
import ROOT
from ROOT import TFile, TH1F, TMath, TCanvas, TLegend, gROOT, gPad
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='mjj')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
xtitle = options.xtitle
rebin = options.rebin

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

mass = [1000,1500,2000,2500,3000]
color = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kGreen+1,ROOT.kMagenta]
hist      = []
## ---- CERN -------
PATH = 'root://eoscms//eos/cms/store/cmst3/group/das2014/EXODijetsLE/'
## ---- FNAL -------
# PATH = '/eos/uscms/store/user/cmsdas/2014/EXODijetsLE/'
#---- open the files --------------------
i = 0
for im in mass:
  inf = TFile.Open(PATH+'dijetHisto_RS'+str(im)+'_signal.root')
  print inf.GetName()
  
  h = inf.Get('h_'+var)
  Nev = h.Integral()
  wt = 1.0/Nev
  h.Scale(wt)
  h.Rebin(rebin)
  h.SetDirectory(0)
  h.SetLineColor(color[i])
  h.SetLineWidth(2)
  hist.append(h)
   
  i += 1

#----- Drawing -----------------------
can = TCanvas('can_Signal_'+var,'can_Signal_'+var,900,600)
can.SetRightMargin(0.2)
hist[0].GetYaxis().SetNdivisions(505)
hist[0].GetXaxis().SetTitle(xtitle)
hist[0].GetXaxis().SetRangeUser(xmin,xmax)
hist[0].SetMaximum(1.2*TMath.Max(hist[0].GetBinContent(hist[0].GetMaximumBin()),hist[4].GetBinContent(hist[4].GetMaximumBin())))
hist[0].Draw('hist')
hist[1].Draw('same hist')
hist[2].Draw('same hist')
hist[3].Draw('same hist')
hist[4].Draw('same hist')
#gPad.RedrawAxis()
    
leg = TLegend(0.81,0.65,0.96,0.9)
leg.AddEntry(hist[0],'1.0 TeV','L')
leg.AddEntry(hist[1],'1.5 TeV','L')
leg.AddEntry(hist[2],'2.0 TeV','L')
leg.AddEntry(hist[3],'2.5 TeV','L')
leg.AddEntry(hist[4],'3.0 TeV','L')
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
