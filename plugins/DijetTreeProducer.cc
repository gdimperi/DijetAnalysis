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

#include "CMSROMA/DijetAnalysis/plugins/DijetTreeProducer.h"
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
  srcJetsAK4_         = cfg.getParameter<edm::InputTag>             ("jetsAK4");
  srcJetsAK8_         = cfg.getParameter<edm::InputTag>             ("jetsAK8");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcVrtx_            = cfg.getParameter<edm::InputTag>             ("vtx");
  srcPU_              = cfg.getUntrackedParameter<edm::InputTag>    ("pu",edm::InputTag(""));
  ptMinAK4_           = cfg.getParameter<double>                    ("ptMinAK4");
  ptMinAK8_           = cfg.getParameter<double>                    ("ptMinAK8");
  //mjjMin_             = cfg.getParameter<double>                    ("mjjMin");
  //dEtaMax_            = cfg.getParameter<double>                    ("dEtaMax");
  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"),consumesCollector());
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
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("nJetsAK4"             ,&nJetsAK4_          ,"nJetsAK4_/I");
  outTree_->Branch("htAK4"                ,&htAK4_             ,"htAK4_/F");
  outTree_->Branch("mjjAK4"               ,&mjjAK4_            ,"mjjAK4_/F");
  outTree_->Branch("dEtajjAK4"            ,&dEtajjAK4_         ,"dEtajjAK4_/F");
  outTree_->Branch("dPhijjAK4"            ,&dPhijjAK4_         ,"dPhijjAK4_/F"); 
  outTree_->Branch("nJetsAK8"             ,&nJetsAK8_          ,"nJetsAK8_/I");
  outTree_->Branch("htAK8"                ,&htAK8_             ,"htAK8_/F");
  outTree_->Branch("mjjAK8"               ,&mjjAK8_            ,"mjjAK8_/F");
  outTree_->Branch("dEtajjAK8"            ,&dEtajjAK8_         ,"dEtajjAK8_/F");
  outTree_->Branch("dPhijjAK8"            ,&dPhijjAK8_         ,"dPhijjAK8_/F"); 
  //------------------------------------------------------------------
  ptAK4_             = new std::vector<float>;
  jecAK4_            = new std::vector<float>;
  etaAK4_            = new std::vector<float>;
  phiAK4_            = new std::vector<float>;
  massAK4_           = new std::vector<float>;
  energyAK4_         = new std::vector<float>;
  chfAK4_            = new std::vector<float>;
  nhfAK4_            = new std::vector<float>;
  phfAK4_            = new std::vector<float>;
  mufAK4_            = new std::vector<float>;
  elfAK4_            = new std::vector<float>;
  idLAK4_            = new std::vector<int>;
  idTAK4_            = new std::vector<int>;
  //massPrunedAK4_     = new std::vector<float>;
  //tau1AK4_           = new std::vector<float>;
  //tau2AK4_           = new std::vector<float>;
  //dRAK4_             = new std::vector<float>;
  outTree_->Branch("jetPtAK4"                ,"vector<float>"     ,&ptAK4_);
  outTree_->Branch("jetJecAK4"               ,"vector<float>"     ,&jecAK4_);
  outTree_->Branch("jetEtaAK4"               ,"vector<float>"     ,&etaAK4_);
  outTree_->Branch("jetPhiAK4"               ,"vector<float>"     ,&phiAK4_);
  outTree_->Branch("jetMassAK4"              ,"vector<float>"     ,&massAK4_);
  outTree_->Branch("jetEnergyAK4"            ,"vector<float>"     ,&energyAK4_);
  outTree_->Branch("jetChfAK4"               ,"vector<float>"     ,&chfAK4_);
  outTree_->Branch("jetNhfAK4"               ,"vector<float>"     ,&nhfAK4_);
  outTree_->Branch("jetPhfAK4"               ,"vector<float>"     ,&phfAK4_);
  outTree_->Branch("jetMufAK4"               ,"vector<float>"     ,&mufAK4_);
  outTree_->Branch("jetElfAK4"               ,"vector<float>"     ,&elfAK4_);   
  outTree_->Branch("idLAK4"                  ,"vector<int>"      ,&idLAK4_);   
  outTree_->Branch("idTAK4"                  ,"vector<int>"      ,&idTAK4_);   
  //outTree_->Branch("jetMassPrunedAK4"        ,"vector<float>"     ,&massPrunedAK4_);
  //outTree_->Branch("jetTau1AK4"              ,"vector<float>"     ,&tau1AK4_);
  //outTree_->Branch("jetTau2AK4"              ,"vector<float>"     ,&tau2AK4_);
  //outTree_->Branch("jetDRAK4"                ,"vector<float>"     ,&dRAK4_); 

  ptAK8_             = new std::vector<float>;
  jecAK8_            = new std::vector<float>;
  etaAK8_            = new std::vector<float>;
  phiAK8_            = new std::vector<float>;
  massAK8_           = new std::vector<float>;
  energyAK8_         = new std::vector<float>;
  chfAK8_            = new std::vector<float>;
  nhfAK8_            = new std::vector<float>;
  phfAK8_            = new std::vector<float>;
  mufAK8_            = new std::vector<float>;
  elfAK8_            = new std::vector<float>;
  idLAK8_            = new std::vector<int>;
  idTAK8_            = new std::vector<int>;
  //massPrunedAK8_     = new std::vector<float>;
  //tau1AK8_           = new std::vector<float>;
  //tau2AK8_           = new std::vector<float>;
  //dRAK8_             = new std::vector<float>;
  outTree_->Branch("jetPtAK8"                ,"vector<float>"     ,&ptAK8_);
  outTree_->Branch("jetJecAK8"               ,"vector<float>"     ,&jecAK8_);
  outTree_->Branch("jetEtaAK8"               ,"vector<float>"     ,&etaAK8_);
  outTree_->Branch("jetPhiAK8"               ,"vector<float>"     ,&phiAK8_);
  outTree_->Branch("jetMassAK8"              ,"vector<float>"     ,&massAK8_);
  outTree_->Branch("jetEnergyAK8"            ,"vector<float>"     ,&energyAK8_);
  outTree_->Branch("jetChfAK8"               ,"vector<float>"     ,&chfAK8_);
  outTree_->Branch("jetNhfAK8"               ,"vector<float>"     ,&nhfAK8_);
  outTree_->Branch("jetPhfAK8"               ,"vector<float>"     ,&phfAK8_);
  outTree_->Branch("jetMufAK8"               ,"vector<float>"     ,&mufAK8_);
  outTree_->Branch("jetElfAK8"               ,"vector<float>"     ,&elfAK8_);   
  outTree_->Branch("idLAK8"                  ,"vector<int>"      ,&idLAK8_);   
  outTree_->Branch("idTAK8"                  ,"vector<int>"      ,&idTAK8_);   
  //outTree_->Branch("jetMassPrunedAK8"        ,"vector<float>"     ,&massPrunedAK8_);
  //outTree_->Branch("jetTau1AK8"              ,"vector<float>"     ,&tau1AK8_);
  //outTree_->Branch("jetTau2AK8"              ,"vector<float>"     ,&tau2AK8_);
  //outTree_->Branch("jetDRAK8"                ,"vector<float>"     ,&dRAK8_); 
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

  delete ptAK4_;
  delete jecAK4_;
  delete etaAK4_;
  delete phiAK4_;
  delete massAK4_;
  delete energyAK4_;
  delete chfAK4_;
  delete nhfAK4_;
  delete phfAK4_;
  delete mufAK4_;
  delete elfAK4_;
  delete idLAK4_;
  delete idTAK4_;
  //delete massPrunedAK4_;
  //delete tau1AK4_;
  //delete tau2AK4_;
  //delete dRAK4_;

  delete ptAK8_;
  delete jecAK8_;
  delete etaAK8_;
  delete phiAK8_;
  delete massAK8_;
  delete energyAK8_;
  delete chfAK8_;
  delete nhfAK8_;
  delete phfAK8_;
  delete mufAK8_;
  delete elfAK8_;
  delete idLAK8_;
  delete idTAK8_;
  //delete massPrunedAK8_;
  //delete tau1AK8_;
  //delete tau2AK8_;
  //delete dRAK8_;
  
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<edm::View<pat::Jet> > jetsAK4;
  iEvent.getByLabel(srcJetsAK4_,jetsAK4);
  edm::View<pat::Jet> pat_jetsAK4 = *jetsAK4;

  edm::Handle<edm::View<pat::Jet> > jetsAK8;
  iEvent.getByLabel(srcJetsAK8_,jetsAK8);
  edm::View<pat::Jet> pat_jetsAK8 = *jetsAK8;

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

    // Event
    met_    = (*met)[0].et();
    if ((*met)[0].sumEt() > 0) {
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
    }
    nVtx_   = recVtxs->size();
    run_    = iEvent.id().run();
    evt_    = iEvent.id().event();
    lumi_   = iEvent.id().luminosityBlock();
    
    // AK4
    nJetsAK4_ = 0;
    float htAK4(0.0);
    vector<TLorentzVector> vP4AK4;
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jetsAK4.begin();ijet != pat_jetsAK4.end(); ++ijet) { 
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm    = ijet->chargedHadronMultiplicity();
      int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
      float eta  = fabs(ijet->eta());
      float pt   = ijet->pt();
      int idL   = (npr>1 && phf<0.99 && nhf<0.99);
      int idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      if (pt > ptMinAK4_) {
        htAK4 += pt;
        nJetsAK4_++;

        vP4AK4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
        chfAK4_           ->push_back(chf);
        nhfAK4_           ->push_back(nhf);
        phfAK4_           ->push_back(phf);
        elfAK4_           ->push_back(elf);
        mufAK4_           ->push_back(muf);
        jecAK4_           ->push_back(1./ijet->jecFactor(0));
        ptAK4_            ->push_back(pt);
        phiAK4_           ->push_back(ijet->phi());
        etaAK4_           ->push_back(ijet->eta());
        massAK4_          ->push_back(ijet->mass());
        energyAK4_        ->push_back(ijet->energy());
	idLAK4_           ->push_back(idL);
	idTAK4_           ->push_back(idT);
        //tau1AK4_          ->push_back(ijet->userFloat("tau1"));
        //tau2AK4_          ->push_back(ijet->userFloat("tau2"));

	//---- match with the pruned jet collection -----
        // double dRmin(1000);
        // double auxm(0.0);
        // for(edm::View<pat::Jet>::const_iterator ijetpr = pat_jetsAK8.begin();ijetpr != pat_jetsAK8.end(); ++ijetpr) { 
        //   float dR = deltaR(ijet->eta(),ijet->phi(),ijetpr->eta(),ijetpr->phi());
        //   if (dR < dRmin) {
        //     auxm = ijetpr->mass();
        //     dRmin = dR;
        //   } 
        // } 
        // massPruned_->push_back(auxm);
        // dR_->push_back(dRmin);

      }// matching with pruned jets
    }// jet loop  
    htAK4_     = htAK4;
    if (nJetsAK4_ > 1) { //assuming jets are ordered by pt in the pat collection
      mjjAK4_    = (vP4AK4[0]+vP4AK4[1]).M();
      dEtajjAK4_ = fabs((*etaAK4_)[0]-(*etaAK4_)[1]); 
      dPhijjAK4_ = fabs(deltaPhi((*phiAK4_)[0],(*phiAK4_)[1]));
    }


    // AK8
    nJetsAK8_ = 0;
    float htAK8(0.0);
    vector<TLorentzVector> vP4AK8;
    for(edm::View<pat::Jet>::const_iterator ijet = pat_jetsAK8.begin();ijet != pat_jetsAK8.end(); ++ijet) { 
      double chf = ijet->chargedHadronEnergyFraction();
      double nhf = ijet->neutralHadronEnergyFraction() + ijet->HFHadronEnergyFraction();
      double phf = ijet->photonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double elf = ijet->electronEnergy()/(ijet->jecFactor(0) * ijet->energy());
      double muf = ijet->muonEnergy()/(ijet->jecFactor(0) * ijet->energy());
      int chm    = ijet->chargedHadronMultiplicity();
      int npr    = ijet->chargedMultiplicity() + ijet->neutralMultiplicity(); 
      float eta  = fabs(ijet->eta());
      float pt   = ijet->pt();
      int idL   = (npr>1 && phf<0.99 && nhf<0.99);
      int idT   = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      if (pt > ptMinAK8_) {
        htAK8 += pt;
        nJetsAK8_++;

        vP4AK8.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
        chfAK8_           ->push_back(chf);
        nhfAK8_           ->push_back(nhf);
        phfAK8_           ->push_back(phf);
        elfAK8_           ->push_back(elf);
        mufAK8_           ->push_back(muf);
        jecAK8_           ->push_back(1./ijet->jecFactor(0));
        ptAK8_            ->push_back(pt);
        phiAK8_           ->push_back(ijet->phi());
        etaAK8_           ->push_back(ijet->eta());
        massAK8_          ->push_back(ijet->mass());
        energyAK8_        ->push_back(ijet->energy());
	idLAK8_           ->push_back(idL);
	idTAK8_           ->push_back(idT);
        //tau1AK8_          ->push_back(ijet->userFloat("tau1"));
        //tau2AK8_          ->push_back(ijet->userFloat("tau2"));

	//---- match with the pruned jet collection -----
        // double dRmin(1000);
        // double auxm(0.0);
        // for(edm::View<pat::Jet>::const_iterator ijetpr = pat_jetsAK8.begin();ijetpr != pat_jetsAK8.end(); ++ijetpr) { 
        //   float dR = deltaR(ijet->eta(),ijet->phi(),ijetpr->eta(),ijetpr->phi());
        //   if (dR < dRmin) {
        //     auxm = ijetpr->mass();
        //     dRmin = dR;
        //   } 
        // } 
        // massPruned_->push_back(auxm);
        // dR_->push_back(dRmin);

      }// matching with pruned jets
    }// jet loop  
    htAK8_     = htAK8;
    if (nJetsAK8_ > 1) { //assuming jets are ordered by pt in the pat collection
      mjjAK8_    = (vP4AK8[0]+vP4AK8[1]).M();
      dEtajjAK8_ = fabs((*etaAK8_)[0]-(*etaAK8_)[1]); 
      dPhijjAK8_ = fabs(deltaPhi((*phiAK8_)[0],(*phiAK8_)[1]));
    }


    //---- Fill Tree ---
    //if (mjjAK4_ > mjjMin_ && dEtajjAK4_ < dEtaMax_) {
    outTree_->Fill();     
    //}
    //------------------

  }// if vtx
}//end analyze for each event

//////////////////////////////////////////////////////////////////////////////////////////
void DijetTreeProducer::initialize()
{
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  met_            = -999;
  metSig_         = -999;
  nJetsAK4_          = -999;
  htAK4_             = -999;
  mjjAK4_            = -999; 
  dEtajjAK4_         = -999; 
  dPhijjAK4_         = -999;
  ptAK4_             ->clear();
  etaAK4_            ->clear();
  phiAK4_            ->clear();
  massAK4_           ->clear();
  energyAK4_         ->clear();
  chfAK4_            ->clear();
  nhfAK4_            ->clear();
  phfAK4_            ->clear();
  elfAK4_            ->clear();
  mufAK4_            ->clear();
  jecAK4_            ->clear();
  jecAK4_            ->clear();
  idLAK4_            ->clear();
  idTAK4_            ->clear();
  //massPrunedAK4_     ->clear();
  //tau1AK4_           ->clear();
  //tau2AK4_           ->clear();
  //dRAK4_             ->clear();

  nJetsAK8_          = -999;
  htAK8_             = -999;
  mjjAK8_            = -999; 
  dEtajjAK8_         = -999; 
  dPhijjAK8_         = -999;
  ptAK8_             ->clear();
  etaAK8_            ->clear();
  phiAK8_            ->clear();
  massAK8_           ->clear();
  energyAK8_         ->clear();
  chfAK8_            ->clear();
  nhfAK8_            ->clear();
  phfAK8_            ->clear();
  elfAK8_            ->clear();
  mufAK8_            ->clear();
  jecAK8_            ->clear();
  jecAK8_            ->clear();
  idLAK8_            ->clear();
  idTAK8_            ->clear();
  //massPrunedAK8_     ->clear();
  //tau1AK8_           ->clear();
  //tau2AK8_           ->clear();
  //dRAK8_             ->clear();

  triggerResult_     ->clear();
  //----- MC -------
  npu_ = -999;
}
//////////////////////////////////////////////////////////////////////////////////////////
DijetTreeProducer::~DijetTreeProducer() 
{
}

DEFINE_FWK_MODULE(DijetTreeProducer);
