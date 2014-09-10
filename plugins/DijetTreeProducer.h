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
    double ptMinAK4_,ptMinAK8_,ptMinCA8_;//mjjMin_,,dEtaMax_;
    edm::InputTag srcJetsAK4_,srcJetsAK8_, srcJetsCA8_,srcMET_,srcPU_,srcVrtx_, srcGenInfo_, srcGenJetsAK4_, srcGenJetsAK8_, srcGenJetsCA8_;
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
    int   nJetsAK4_, nJetsAK8_, nJetsCA8_, nGenJetsAK4_, nGenJetsAK8_, nGenJetsCA8_;
    float rho_,met_,metSig_;
    float htAK4_,mjjAK4_,dEtajjAK4_,dPhijjAK4_;
    float htAK8_,mjjAK8_,dEtajjAK8_,dPhijjAK8_;
    float htCA8_,mjjCA8_,dEtajjCA8_,dPhijjCA8_;
    std::vector<bool> *triggerResult_;
    //---- jet variables --------------
    std::vector<float> *ptAK4_,*jecAK4_,*etaAK4_,*phiAK4_,*massAK4_,*energyAK4_,*chfAK4_,*nhfAK4_,*phfAK4_,*elfAK4_,*mufAK4_;
    std::vector<int> *idLAK4_,*idTAK4_;
    //std::vector<float> *cutbasedJetId_, *fullJetId_, *fullJetDiscriminant_;
    std::vector<float> *ptGenAK4_,*etaGenAK4_,*phiGenAK4_,*massGenAK4_,*energyGenAK4_;

    std::vector<float> *ptAK8_,*jecAK8_,*etaAK8_,*phiAK8_,*massAK8_,*energyAK8_,*chfAK8_,*nhfAK8_,*phfAK8_,*elfAK8_,*mufAK8_, *massPrunedAK8_, *dR_AK8_,*tau1AK8_,*tau2AK8_, *tau3AK8_ ;
    std::vector<int> *idLAK8_,*idTAK8_;
    std::vector<float> *ptGenAK8_,*etaGenAK8_,*phiGenAK8_,*massGenAK8_,*energyGenAK8_;

    std::vector<float> *ptCA8_,*jecCA8_,*etaCA8_,*phiCA8_,*massCA8_,*energyCA8_,*chfCA8_,*nhfCA8_,*phfCA8_,*elfCA8_,*mufCA8_, *massPrunedCA8_, *dR_CA8_,*tau1CA8_,*tau2CA8_, *tau3CA8_ ;
    std::vector<int> *idLCA8_,*idTCA8_;
    std::vector<float> *ptGenCA8_,*etaGenCA8_,*phiGenCA8_,*massGenCA8_,*energyGenCA8_;

    //---- MC variables ---------------
    double npu_; 
    float ptHat_; 
    int processID_;
    double weight_;
};

#endif
