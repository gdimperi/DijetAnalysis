#!usr/bin/python
import ROOT
from ROOT import TFile, TH1F, THStack, TCanvas, TMath, gROOT, gPad
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='mjj')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
xtitle = options.xtitle
rebin = options.rebin
logy = options.logy

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

fileNames = ['QCD500','QCD1000','RS2000','data']
xsections = [8426.,204.,4.083e-3,1.]
colorF    = [ROOT.kBlue-9,ROOT.kBlue-8,ROOT.kWhite,ROOT.kBlack]
colorL    = [ROOT.kBlack,ROOT.kBlack,ROOT.kRed,ROOT.kBlack]
hist      = []
LUMI      = 19800.

## ---- CERN -------
PATH = 'root://eoscms//eos/cms/store/cmst3/group/das2014/EXODijetsLE/'
## ---- FNAL -------
# PATH = '/eos/uscms/store/user/cmsdas/2014/EXODijetsLE/'

#---- open the files --------------------
i_f = 0
for f in fileNames:
  inf = TFile.Open(PATH+'dijetHisto_'+f+'_signal.root')
  print inf.GetName()
  
  Nev = inf.Get('TriggerPass').GetBinContent(1)
  wt = 1.0
  if i_f < 3:
    wt = LUMI*xsections[i_f]/Nev
  
  h = inf.Get('h_'+var)
  h.Scale(wt)
  h.Rebin(rebin)
  h.SetDirectory(0)
  h.SetFillColor(colorF[i_f])
  h.SetLineColor(colorL[i_f])
  h.SetMarkerColor(colorL[i_f])
  hist.append(h)
   
  i_f += 1

NQCD = hist[0].Integral()+hist[1].Integral()
NDAT = hist[3].Integral()
kFactor = NDAT/NQCD
print kFactor

hist[0].Scale(kFactor)
hist[1].Scale(kFactor)

histQCD = hist[0].Clone('histQCD')
histQCD.Add(hist[1])

hsQCD = THStack('QCD','QCD')

hsQCD.Add(hist[0])
hsQCD.Add(hist[1])

#----- Drawing -----------------------
can = TCanvas('can_'+var,'can_'+var,900,600)
if logy:
  gPad.SetLogy()
hAux = hist[3].Clone('aux')
hAux.Reset()
hAux.GetXaxis().SetRangeUser(xmin,xmax)
hAux.GetXaxis().SetTitle(xtitle)
hAux.SetMaximum(1.2*TMath.Max(hist[3].GetBinContent(hist[3].GetMaximumBin()),histQCD.GetBinContent(histQCD.GetMaximumBin())))
hAux.SetMinimum(0.01)
hAux.Draw()
hsQCD.Draw('same hist')
hist[3].Draw('same E')
hist[2].Draw('same hist')
gPad.RedrawAxis()
    
#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
