#!usr/bin/python
import ROOT, math
from ROOT import TFile, TH2F, TCanvas, TPad
from ROOT import gROOT, gPad
from array import array 
from ROOT import RooArgSet,  RooBernstein, RooFit, RooMsgService
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

RooMsgService.instance().setSilentMode(ROOT.kTRUE)
RooMsgService.instance().setStreamStatus(0,ROOT.kFALSE)
RooMsgService.instance().setStreamStatus(1,ROOT.kFALSE)

filename = 'RS1500_workspace.root'
sigXS = 3.298e-2
sigEff = 0.307
lumi = 19800
if useSub:
  filename = 'RS1500_sub_workspace.root'
  sigEff = 0.066

#--- access the workspace and the objects it contains -----
inf = TFile.Open(filename)
w = inf.Get('w')
bkg_pdf = w.pdf('background')
sig_pdf = w.pdf('signal')
x = w.var('mjj')
bkg_hist = w.data('data_obs')
#--- compute the total signal yield ----
sigYield = sigXS*sigEff*lumi
#--- define the minimum and maximum mass cuts -----
xmin = [1300,1350,1400,1450,1500,1550]
xmax = [1550,1600,1650,1700,1750,1800,1850]
vxmin = array('d',xmin)
vxmax = array('d',xmax)
hSig = TH2F('significance','signicificance',len(vxmin)-1,vxmin,len(vxmax)-1,vxmax)
i=0
for m1 in xmin[:len(vxmin)-1]:
  j=0
  for m2 in xmax[:len(vxmax)-1]: 
    #--- define the mass window ------
    x.setRange('window',m1,m2)
    #--- compute the fraction of the background and signal PDF in the window -----
    argset = RooArgSet(x)
    fB = bkg_pdf.createIntegral(argset,RooFit.NormSet(argset),RooFit.Range('window')).getVal()
    fS = sig_pdf.createIntegral(argset,RooFit.NormSet(argset),RooFit.Range('window')).getVal()
    #--- compute the yields ----------
    B = fB*bkg_hist.sumEntries()
    S = fS*sigYield
    #--- compute the significance -------
    sign = round(S/math.sqrt(B+S),2)
    print str(m1)+'-'+str(m2)+' '+str(round(B,1))+' '+str(round(S,1))+' '+str(sign)
    hSig.SetBinContent(i+1,j+1,sign)
    j+=1
  i+=1 

canvasName = '2DSign'
if useSub:
  canvasName = '2DSign_sub'

can = TCanvas(canvasName,canvasName,600,600)
gPad.SetGridx()
gPad.SetGridy()
hSig.GetXaxis().SetNdivisions(505)
hSig.GetXaxis().SetTitle('Minimum m_{jj} (GeV)')
hSig.GetYaxis().SetTitle('Maximum m_{jj} (GeV)')
hSig.Draw('text col')

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]





