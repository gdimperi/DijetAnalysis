#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "analysis/DijetAnalysis/plugins/DijetTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

DijetTreeProducer::DijetTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>             ("jets");
  srcJetsPruned_      = cfg.getParameter<edm::InputTag>             ("jetsPruned");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcVrtx_            = cfg.getParameter<edm::InputTag>             ("vtx");
  srcPU_              = cfg.getUntrackedParameter<edm::InputTag>    ("pu",edm::InputTag(""));
  mjjMin_             = cfg.getParameter<double>                    ("mjjMin");
  ptMin_              = cfg.getParameter<double>                    ("ptMin");
  dEtaMax_            = cfg.getParameter<double>                    ("dEtaMax");
  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"));
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }

  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);
  
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("mjj"                  ,&mjj_               ,"mjj_/F");
  outTree_->Branch("dEtajj"               ,&dEtajj_            ,"dEtajj_/F");
  outTree_->Branch("dPhijj"               ,&dPhijj_            ,"dPhijj_/F"); 
  //------------------------------------------------------------------
  pt_             = new std::vector<float>;
  jec_            = new std::vector<float>;
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  mass_           = new std::vector<float>;
  massPruned_     = new std::vector<float>;
  tau1_           = new std::vector<float>;
  tau2_           = new std::vector<float>;
  dR_             = new std::vector<float>;
  energy_         = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetJec"               ,"vector<float>"     ,&jec_);
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetMassPruned"        ,"vector<float>"     ,&massPruned_);
  outTree_->Branch("jetTau1"              ,"vector<float>"     ,&tau1_);
  outTree_->Branch("jetTau2"              ,"vector<float>"     ,&tau2_);
  outTree_->Branch("jetDR"                ,"vector<float>"     ,&dR_); 
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);   
  //------------------------------------------------------------------
  triggerResult_ = new std::vector<bool>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  //------------------- MC ---------------------------------
  outTree_->Branch("npu"                  ,&npu_               ,"npu_/I");
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::endJob() 
{  
  delete triggerResult_;
  delete pt_;
  delete jec_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete massPruned_;
  delete tau1_;
  delete tau2_;
  delete dR_;
  delete energy_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<edm::View<pat::Jet> > jetsPruned;
  iEvent.getByLabel(srcJetsPruned_,jetsPruned);
  edm::View<pat::Jet> pat_jets_pruned = *jetsPruned;

  edm::Handle<edm::View<MET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(srcVrtx_,recVtxs);
  
  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    iEvent.getByLabel(srcPU_,PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
    for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() == 0) {
        npu_ = PUI->getTrueNumInteractions();
      }
    }
  }// if MC
  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  if (triggerCache_.setEvent(iEvent,iSetup)) {
    for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++) {
      bool result(false);
      if (vtriggerSelector_[itrig]) {
        if (triggerCache_.configurationUpdated()) {
          vtriggerSelector_[itrig]->init(triggerCache_);
        }
        result = (*(vtriggerSelector_[itrig]))(triggerCache_);
      }
      if (result) {
        triggerPassHisto_->Fill(vtriggerAlias_[itrig].c_str(),1);
      }
      triggerResult_->push_back(result);
    }
  }
     
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  
  if (cut_vtx) {
    nJets_ = 0;
    float ht(0.0);
    vector<TLorentzVector> vP4;
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin();ijet != pat_jets.end(); ++ijet) { 
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm    = ijet->chargedHadronMultiplicity();
      int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
      float eta  = fabs(ijet->eta());
      float pt   = ijet->pt();
      bool idL   = (npr>1 && phf<0.99 && nhf<0.99);
      bool idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      if (idT && pt > ptMin_) {
        ht += pt;
        vP4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
        chf_           ->push_back(chf);
        nhf_           ->push_back(nhf);
        phf_           ->push_back(phf);
        elf_           ->push_back(elf);
        muf_           ->push_back(muf);
        jec_           ->push_back(1./ijet->jecFactor(0));
        pt_            ->push_back(pt);
        phi_           ->push_back(ijet->phi());
        eta_           ->push_back(ijet->eta());
        mass_          ->push_back(ijet->mass());
        energy_        ->push_back(ijet->energy());
        tau1_          ->push_back(ijet->userFloat("tau1"));
        tau2_          ->push_back(ijet->userFloat("tau2"));
          
        //---- match with the pruned jet collection -----
        double dRmin(1000);
        double auxm(0.0);
        for(edm::View<pat::Jet>::const_iterator ijetpr = pat_jets_pruned.begin();ijetpr != pat_jets_pruned.end(); ++ijetpr) { 
          float dR = deltaR(ijet->eta(),ijet->phi(),ijetpr->eta(),ijetpr->phi());
          if (dR < dRmin) {
            auxm = ijetpr->mass();
            dRmin = dR;
          } 
        } 
        massPruned_->push_back(auxm);
        dR_->push_back(dRmin);
        nJets_++;
      }// matching with pruned jets
    }// jet loop  
    ht_     = ht;
    met_    = (*met)[0].et();
    if ((*met)[0].sumEt() > 0) {
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
    }
    nVtx_   = recVtxs->size();
    run_    = iEvent.id().run();
    evt_    = iEvent.id().event();
    lumi_   = iEvent.id().luminosityBlock();
    if (nJets_ > 1) { 
      mjj_    = (vP4[0]+vP4[1]).M();
      dEtajj_ = fabs((*eta_)[0]-(*eta_)[1]); 
      dPhijj_ = fabs(deltaPhi((*phi_)[0],(*phi_)[1]));
      if (mjj_ > mjjMin_ && dEtajj_ < dEtaMax_) {
        outTree_->Fill();     
      }
    }// if nJets > 1
  }// if vtx
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nJets_          = -999;
  met_            = -999;
  metSig_         = -999;
  ht_             = -999;
  mjj_            = -999; 
  dEtajj_         = -999; 
  dPhijj_         = -999;
  pt_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  mass_           ->clear();
  massPruned_     ->clear();
  tau1_           ->clear();
  tau2_           ->clear();
  dR_             ->clear();
  energy_         ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  jec_            ->clear();
  triggerResult_  ->clear();
  //----- MC -------
  npu_ = -999;
}
//////////////////////////////////////////////////////////////////////////////////////////
DijetTreeProducer::~DijetTreeProducer() 
{
}

DEFINE_FWK_MODULE(DijetTreeProducer);
