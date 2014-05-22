#!usr/bin/python
import ROOT
from ROOT import TFile, TH1F, TCanvas, TMath, TEfficiency, TF1, TLegend, TLine, gROOT, gPad
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--set",action="store",type="string",dest="set",default='data')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=500)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1500)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=10)

(options, args) = parser.parse_args()

set = options.set
xmin = options.xmin
xmax = options.xmax
rebin = options.rebin

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

## ---- CERN -------
PATH = 'root://eoscms//eos/cms/store/cmst3/group/das2014/EXODijetsLE/'
## ---- FNAL -------
# PATH = '/eos/uscms/store/user/cmsdas/2014/EXODijetsLE/'

infRef = TFile.Open(PATH+'dijetHisto_'+set+'_ref.root')
infSig = TFile.Open(PATH+'dijetHisto_'+set+'_refSig.root')

hRef = TH1F(infRef.Get('h_mjj'))
hRef.Rebin(rebin)
hRef.SetDirectory(0)
hSig = TH1F(infSig.Get('h_mjj'))
hSig.Rebin(rebin)
hSig.SetDirectory(0)

can1 = TCanvas('can_mjj_trig','can_mjj_trig',900,600)
hRef.GetXaxis().SetTitle('m_{jj} (GeV)')
hRef.GetYaxis().SetTitle('Events')
hRef.GetXaxis().SetRangeUser(xmin,xmax)
hRef.SetFillColor(ROOT.kGray)
hRef.Draw('hist')
hSig.Draw('sameE')
leg = TLegend(0.6,0.7,0.9,0.9)
leg.AddEntry(hSig,'Ref. + Sig. trigger','LP')
leg.AddEntry(hRef,'Ref. trigger','F')
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.05)
leg.Draw()

can2 = TCanvas('can_TrigEff','can_TrigEff',900,600)
eff = TEfficiency('efficiency','efficiency',hRef.GetNbinsX(),0,2000)
eff.SetPassedHistogram(hSig,'f')
eff.SetTotalHistogram(hRef,'f')
eff.SetLineColor(ROOT.kBlack)
eff.SetMarkerColor(ROOT.kBlack)
eff.SetMarkerStyle(20)

fit = TF1('fit','(1-[2])+[2]*TMath::Erf((x-[0])/[1])',xmin,xmax)
fit.SetParameters(800,100,0.3)
fit.SetLineColor(ROOT.kBlue)
eff.Fit(fit,'RQ')

line = TF1('line','1',xmin,xmax)
line.GetXaxis().SetTitle('m_{jj} (GeV)')
line.GetYaxis().SetTitle('Trigger Efficiency')
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(2)
line.SetMinimum(0.3)
line.SetMaximum(1.1)
line.Draw()
eff.Draw('samePE1')

gPad.Update()
x0 = fit.GetX(0.99)
cut = TLine(x0,gPad.GetFrame().GetY1(),x0,gPad.GetFrame().GetY2())
cut.SetLineColor(ROOT.kRed)
cut.SetLineStyle(2)
cut.Draw()

print '99% efficiency point: '+str(round(x0,1))+' GeV'

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
