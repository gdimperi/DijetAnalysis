#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TH1F, TCanvas, TPad
from ROOT import gROOT, gPad 
from ROOT import RooRealVar, RooDataHist, RooPlot, RooArgList, RooArgSet,  RooBernstein, RooCBShape, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService
from setTDRStyle import setTDRStyle

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-s","--fitSig",action="store_true",default=False,dest="fitSig")
parser.add_option("-d","--fitDat",action="store_true",default=False,dest="fitDat")
parser.add_option("-m","--mass",action="store",type="int",dest="mass",default=2000)


parser.add_option("--lumi",action="store",type="float",dest="lumi",default=19.8)
parser.add_option("--sigEff",action="store",type="float",dest="sigEff",default=0.2941)
parser.add_option("--sigXS",action="store",type="float",dest="sigXS",default=0.004083e3)

parser.add_option("--useSub",action="store_true",default=False,dest="useSub")

(options, args) = parser.parse_args()

mass = options.mass
fitSig = options.fitSig
fitDat = options.fitDat
useSub = options.useSub

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')

RooMsgService.instance().setSilentMode(ROOT.kTRUE)
RooMsgService.instance().setStreamStatus(0,ROOT.kFALSE)
RooMsgService.instance().setStreamStatus(1,ROOT.kFALSE)

# -----------------------------------------
# get histograms
## ---- CERN -------
PATH = 'root://eoscms//eos/cms/store/cmst3/group/das2014/EXODijetsLE/'
## ---- FNAL -------
# PATH = '/eos/uscms/store/user/cmsdas/2014/EXODijetsLE/'

filenameSig = PATH+'dijetHisto_RS'+str(mass)+'_signal.root'
filenameDat = PATH+'dijetHisto_data_signal.root'

inputHistName = 'h_mjj'

if useSub:
 inputHistName = 'h_sub_mjj'

infSig = TFile.Open(filenameSig)
hSig   = infSig.Get(inputHistName)
hSig.Rebin(20)

infDat = TFile.Open(filenameDat)
hDat   = infDat.Get(inputHistName)
#hDat.Rebin(20)

# -----------------------------------------
# define observable
x = RooRealVar('mjj','mjj',900,4500)

if fitSig: 

    # define parameters for signal fit
    m = RooRealVar('mean','mean',float(mass),float(mass)-200,float(mass)+200)
    s = RooRealVar('sigma','sigma',0.1*float(mass),0,10000)
    a = RooRealVar('alpha','alpha',1,-10,10)
    n = RooRealVar('n','n',1,0,100)
    sig = RooCBShape('sig','sig',x,m,s,a,n)        

    p  = RooRealVar('p','p',1,0,5)
    x0 = RooRealVar('x0','x0',1000,100,5000)

    bkg = RooGenericPdf('bkg','1/(exp(pow(@0/@1,@2))+1)',RooArgList(x,x0,p))

    fsig= RooRealVar('fsig','fsig',0.5,0.,1.)
    signal = RooAddPdf('signal','signal',sig,bkg,fsig)

    # -----------------------------------------
    # fit signal
    canSname = 'can_Mjj'+str(mass)
    if useSub:
      canSname = 'can_Sub_Mjj'+str(mass)
    canS = TCanvas(canSname,canSname,900,600)
    #gPad.SetLogy() 

    roohistSig = RooDataHist('roohist','roohist',RooArgList(x),hSig)

    signal.fitTo(roohistSig)
    frame = x.frame()
    roohistSig.plotOn(frame)
    signal.plotOn(frame)
    signal.plotOn(frame,RooFit.Components('bkg'),RooFit.LineColor(ROOT.kRed),RooFit.LineWidth(2),RooFit.LineStyle(ROOT.kDashed))
    frame.GetXaxis().SetRangeUser(900,4500)
    frame.GetXaxis().SetTitle('m_{jj} (GeV)')
    frame.Draw()

    parsSig = signal.getParameters(roohistSig)
    parsSig.setAttribAll('Constant', True)

if fitDat: 

    # -----------------------------------------
    # define parameters for background
    NBINS = 180
    p1 = RooRealVar('p1','p1',7,1,10)
    p2 = RooRealVar('p2','p2',5,1,10)
    p3 = RooRealVar('p3','p3',0.03,0.01,0.07)

    background = RooGenericPdf('background','(pow(1-@0/8000,@1)/pow(@0/8000,@2+@3*log(@0/8000)))',RooArgList(x,p1,p2,p3))
    roohistBkg = RooDataHist('roohist','roohist',RooArgList(x),hDat)
    res = background.fitTo(roohistBkg)

    # -----------------------------------------
    # plot background
    canBname = 'can_Mjj_Data'
    if useSub:
      canBname = 'can_Sub_Mjj_Data'
    canB = TCanvas(canBname,canBname,900,600)
    gPad.SetLogy() 
    canB.cd(1).SetBottomMargin(0.4)

    frame1 = x.frame()
    frame2 = x.frame()
    roohistBkg.plotOn(frame1,RooFit.Binning(NBINS))
    background.plotOn(frame1)
    hpull = frame1.pullHist()
    frame2.addPlotable(hpull,'p')

    frame1.SetMinimum(0.5)
    frame1.GetXaxis().SetTitle('')
    if useSub:
      frame1.GetXaxis().SetRangeUser(900,2500)
    frame1.GetXaxis().SetLabelSize(0.0)
    frame1.GetYaxis().SetTickLength(0.06)
    frame1.Draw()

    pad = TPad('pad','pad',0.,0.,1.,1.)
    pad.SetTopMargin(0.6)
    pad.SetFillColor(0)
    pad.SetFillStyle(0)
    pad.Draw()
    pad.cd(0)
    frame2.SetMinimum(-5)
    frame2.SetMaximum(5)
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetXaxis().SetTitleOffset(0.9)
    frame2.GetYaxis().SetTitleOffset(0.8)
    frame2.GetYaxis().SetTickLength(0.06)
    frame2.GetYaxis().SetTitleSize(0.05)
    frame2.GetYaxis().SetLabelSize(0.03)
    frame2.GetYaxis().SetTitle('(Data-Fit)/Error')
    frame2.GetXaxis().SetTitle('m_{jj} (GeV)')
    if useSub:
      frame2.GetXaxis().SetRangeUser(900,2500)
    frame2.Draw();

    parsBkg = background.getParameters(roohistBkg)
    parsBkg.setAttribAll('Constant', True)

if fitSig and fitDat:
    
    # -----------------------------------------
    # write everything to a workspace to make a datacard
    dcFN = 'RS'+str(mass)+'_datacard.txt'
    wsFN = 'RS'+str(mass)+'_workspace.root'
    
    if useSub:
      dcFN = 'RS'+str(mass)+'_sub_datacard.txt'
      wsFN = 'RS'+str(mass)+'_sub_workspace.root'

    nObs = roohistBkg.sumEntries();
    
    w = RooWorkspace('w','workspace')
    getattr(w,'import')(signal)
    getattr(w,'import')(background)
    getattr(w,'import')(roohistBkg,RooFit.Rename("data_obs"))  
    w.Print()
    w.writeToFile(wsFN)
    
    # -----------------------------------------
    # write a datacard
    LUMI = options.lumi
    signalCrossSection = options.sigXS
    signalEfficiency = options.sigEff
    ExpectedSignalRate = signalCrossSection*LUMI*signalEfficiency

    datacard = open(dcFN,'w')
    datacard.write('imax 1\n')
    datacard.write('jmax 1\n')
    datacard.write('kmax *\n')
    datacard.write('---------------\n')
    datacard.write('shapes * * '+wsFN+' w:$PROCESS\n')
    datacard.write('---------------\n')
    datacard.write('bin 1\n')    
    datacard.write('observation '+str(nObs)+'\n')
    datacard.write('------------------------------\n')
    datacard.write('bin          1          1\n')          
    datacard.write('process      signal     background\n')
    datacard.write('process      0          1\n')          
    datacard.write('rate         '+str(ExpectedSignalRate)+'         '+str(nObs)+'\n')
    datacard.write('------------------------------\n')      
                        
#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
