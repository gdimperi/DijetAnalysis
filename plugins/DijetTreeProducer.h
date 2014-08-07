#ifndef DijetTreeProducer_h
#define DijetTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TH1F.h"

class DijetTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit DijetTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~DijetTreeProducer();

  private:  
    void initialize();
    //---- configurable parameters --------   
    double ptMinAK4_,ptMinAK8_;//mjjMin_,,dEtaMax_;
    edm::InputTag srcJetsAK4_,srcJetsAK8_,srcMET_,srcPU_,srcVrtx_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    //---- TRIGGER -------------------------
    triggerExpression::Data triggerCache_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_,*puHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_;
    int   nJetsAK4_, nJetsAK8_;
    float rho_,met_,metSig_;
    float htAK4_,mjjAK4_,dEtajjAK4_,dPhijjAK4_;
    float htAK8_,mjjAK8_,dEtajjAK8_,dPhijjAK8_;
    std::vector<bool> *triggerResult_;
    //---- jet variables --------------
    std::vector<float> *ptAK4_,*jecAK4_,*etaAK4_,*phiAK4_,*massAK4_,*energyAK4_,*chfAK4_,*nhfAK4_,*phfAK4_,*elfAK4_,*mufAK4_;// *massPruned_, *dR_,*tau1_,*tau2_ ;
    std::vector<int> *idLAK4_,*idTAK4_;
    std::vector<float> *ptAK8_,*jecAK8_,*etaAK8_,*phiAK8_,*massAK8_,*energyAK8_,*chfAK8_,*nhfAK8_,*phfAK8_,*elfAK8_,*mufAK8_;// *massPruned_, *dR_,*tau1_,*tau2_ ;
    std::vector<int> *idLAK8_,*idTAK8_;
    //---- MC variables ---------------
    int npu_;
};

#endif
