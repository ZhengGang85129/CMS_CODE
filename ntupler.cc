// -*- C++ -*-
//
// Package:    ntuple_generator/ntupler
// Class:      ntupler
// 
/**\class ntupler ntupler.cc ntuple_generator/ntupler/plugins/ntupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  ZHENG-GANG CHEN
//         Created:  Fri, 07 Aug 2020 08:50:01 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
//////////////////
//-- Frame Works//
//////////////////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

////////////////
// -- Else -- //
////////////////
#include "DataFormats/Math/interface/LorentzVector.h"
//#include "/afs/cern.ch/user/z/zhenggan/work_space/test/CMSSW_11_1_0_pre4/src/project_001/ntupler/plugins/constituents.h"
////////////////
// -- ROOT --///
///////////////
#include <vector>
#include <cmath>
#include <TTree.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/foreach.hpp>
///////////////////////
// -- GenParticle -- //
///////////////////////
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
///////////////////////
// -- HGCAL LAYER -- //
///////////////////////

#include "DataFormats/CaloRecHit/interface/CaloCluster.h" 
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//////////////////
// --CaloJets --//
//////////////////
#include "DataFormats/JetReco/interface/CaloJet.h"
//////////////////
// -- PFJets --///
//////////////////
#include "DataFormats/JetReco/interface/PFJet.h"
//Particle Flow Candidate
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
using namespace std;
using namespace hgcal;
class ntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ntupler(const edm::ParameterSet&);
      ~ntupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<vector<reco::GenParticle>> GenParticlesToken_;
      edm::EDGetTokenT<vector<reco::CaloJet>> CaloJetsToken_;
      edm::EDGetTokenT<vector<reco::PFJet>> AK4PFJetsToken_;
      edm::EDGetTokenT<vector<reco::CaloCluster>> HGCalToken_;
      TTree *mytree = new TTree("Tree","Tree");
      //      reco::GenParticle *genparticle = new reco::GenParticle();
//      reco::CaloJet *Jet = new reco::CaloJet();
      static const int EVTSIZE    =  100000;
      static const int n_constituents   = 300;
      ////////////////////
      //-- Genparticle--//
      ////////////////////
      Int_t genparticle_track ; 
      Int_t genparticle_charge[EVTSIZE];
      Int_t genparticle_pdgId[EVTSIZE];
      Int_t genparticle_status[EVTSIZE];
      Double_t genparticle_p[EVTSIZE];
      Double_t genparticle_E[EVTSIZE];
      Double_t genparticle_Et[EVTSIZE];
      Double_t genparticle_mass[EVTSIZE];
      Double_t genparticle_mt[EVTSIZE];
      Double_t genparticle_px[EVTSIZE];
      Double_t genparticle_py[EVTSIZE];
      Double_t genparticle_pz[EVTSIZE];
      Double_t genparticle_pt[EVTSIZE];
      Double_t genparticle_phi[EVTSIZE];
      Double_t genparticle_theta[EVTSIZE];
      Double_t genparticle_eta[EVTSIZE];
      Double_t genparticle_rapidity[EVTSIZE];
      /////////////////////
      ///-- RecoJets --////
      /////////////////////
      Int_t recojets_track ;
      Int_t recojets_charge[EVTSIZE];
      Int_t recojets_pdgId[EVTSIZE];
      Double_t recojets_p[EVTSIZE];
      Double_t recojets_E[EVTSIZE];
      Double_t recojets_Et[EVTSIZE];
      Double_t recojets_mass[EVTSIZE];
      Double_t recojets_mt[EVTSIZE];
      Double_t recojets_px[EVTSIZE];
      Double_t recojets_py[EVTSIZE];
      Double_t recojets_pz[EVTSIZE];
      Double_t recojets_pt[EVTSIZE];
      Double_t recojets_phi[EVTSIZE];
      Double_t recojets_theta[EVTSIZE];
      Double_t recojets_eta[EVTSIZE];
      Double_t recojets_rapidity[EVTSIZE];
      /////////////////////
      ///-- akPFJets --////
      /////////////////////
      Int_t ak4pfjets_track ;
      Int_t ak4pfjets_charge[EVTSIZE];
      Int_t ak4pfjets_pdgId[EVTSIZE];
      Int_t ak4pfjets_nConstituents[EVTSIZE];
      Int_t ak4pfjets_chargedHadronMultiplicity[EVTSIZE];
      Int_t ak4pfjets_neutralHadronMultiplicity[EVTSIZE];
      Int_t ak4pfjets_neutralMultiplicity[EVTSIZE];
      Int_t ak4pfjets_chargedMultiplicity[EVTSIZE];
      Int_t ak4pfjets_HFEMMultiplicity[EVTSIZE];
      Int_t ak4pfjets_HFHadronMultiplicity[EVTSIZE];
      std::vector<Int_t> ak4pfjets_tagging;
      Double_t ak4pfjets_p[EVTSIZE];
      Double_t ak4pfjets_E[EVTSIZE];
      Double_t ak4pfjets_Et[EVTSIZE];
      Double_t ak4pfjets_mass[EVTSIZE];
      Double_t ak4pfjets_mt[EVTSIZE];
      Double_t ak4pfjets_px[EVTSIZE];
      Double_t ak4pfjets_py[EVTSIZE];
      Double_t ak4pfjets_pz[EVTSIZE];
      Double_t ak4pfjets_pt[EVTSIZE];
      Double_t ak4pfjets_phi[EVTSIZE];
      Double_t ak4pfjets_theta[EVTSIZE];
      Double_t ak4pfjets_eta[EVTSIZE];
      Double_t ak4pfjets_rapidity[EVTSIZE];
      
      /////////////////////////////
      ///-- akPFConstituents --////
      /////////////////////////////
      Int_t PFConstituents_track ;
      std::vector<Int_t> PFConstituents_charge;
      std::vector<Int_t> PFConstituents_pdgId;
      std::vector<Double_t> PFConstituents_p;
      std::vector<Double_t> PFConstituents_E;
      std::vector<Double_t> PFConstituents_Et;
      std::vector<Double_t> PFConstituents_mass;
      std::vector<Double_t> PFConstituents_mt;
      std::vector<Double_t> PFConstituents_px;
      std::vector<Double_t> PFConstituents_py;
      std::vector<Double_t> PFConstituents_pz;
      std::vector<Double_t> PFConstituents_pt;
      std::vector<Double_t> PFConstituents_phi;
      std::vector<Double_t> PFConstituents_theta;
      std::vector<Double_t> PFConstituents_eta;
      std::vector<Double_t> PFConstituents_rapidity;
      std::vector<Double_t> PFConstituents_ecalEnergy;
      std::vector<Double_t> PFConstituents_hcalEnergy;
      ////////////////////
      ///-- HGCALLAY --///
      ////////////////////
      Int_t hgcalclst_track;
      Double_t hgcalclst_x[EVTSIZE] ;
      Double_t hgcalclst_y[EVTSIZE] ;
      Double_t hgcalclst_z[EVTSIZE] ;
      Double_t hgcalclst_E[EVTSIZE] ;
      Double_t hgcalclst_eta[EVTSIZE] ;
      Double_t hgcalclst_phi[EVTSIZE] ;
      Bool_t hgcalclst_tagger[EVTSIZE];
      Double_t tagger_R[EVTSIZE];
      Int_t hgcalclst_layer[EVTSIZE];
      Int_t hgcalclst_size[EVTSIZE] ;
      std::vector<Int_t> PCLU;
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//
//
// constructors and destructor
//M
//
ntupler::ntupler(const edm::ParameterSet& iConfig):
GenParticlesToken_(consumes<vector<reco::GenParticle>>(edm::InputTag("genParticles"))),
//CaloJetsToken_(consumes<vector<reco::CaloJet>>(edm::InputTag("ak4CaloJets"))), 
AK4PFJetsToken_(consumes<vector<reco::PFJet>>(edm::InputTag("ak4PFJets"))), 
HGCalToken_( consumes<vector<reco::CaloCluster>>(edm::InputTag("hgcalLayerClusters")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
}


ntupler::~ntupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  -----------void
void
ntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<vector<reco::GenParticle>> GenParticles;
   iEvent.getByToken( GenParticlesToken_ , GenParticles );
   Handle<vector<reco::PFJet>> ak4pfjets;
   iEvent.getByToken( AK4PFJetsToken_ , ak4pfjets );
   Handle<vector<reco::CaloCluster>> hgcalclst;
   iEvent.getByToken( HGCalToken_ , hgcalclst );
   genparticle_track = 0;
   for(vector<reco::GenParticle>::const_iterator g = GenParticles->begin(); g != GenParticles->end(); ++g)
   {
      genparticle_px[genparticle_track]=g->px();
      genparticle_py[genparticle_track]=g->py();
      genparticle_pz[genparticle_track]=g->pz();
      genparticle_pt[genparticle_track]=g->pt();
      genparticle_Et[genparticle_track]=g->et();
      genparticle_E[genparticle_track]=g->energy();
      genparticle_mass[genparticle_track] = g->mass();
      genparticle_mt[genparticle_track] = g->mt();
      genparticle_phi[genparticle_track] = g->phi();
      genparticle_theta[genparticle_track] = g->theta();
      genparticle_eta[genparticle_track] = g->eta();
      genparticle_rapidity[genparticle_track] = g->rapidity();
      genparticle_charge[genparticle_track] = g->charge();
      genparticle_pdgId[genparticle_track] = g->pdgId();
      genparticle_status[ genparticle_track ]= g->status();
      genparticle_track++ ;
   }
   ak4pfjets_track=0;
   PFConstituents_track=0; 
   PFConstituents_charge.clear();
   PFConstituents_pdgId.clear();
   PFConstituents_p.clear();
   PFConstituents_px.clear();
   PFConstituents_py.clear();
   PFConstituents_pz.clear();
   PFConstituents_pt.clear();
   PFConstituents_E.clear();
   PFConstituents_Et.clear();
   PFConstituents_mass.clear();
   PFConstituents_mt.clear();
   PFConstituents_rapidity.clear();
   PFConstituents_eta.clear();
   PFConstituents_phi.clear();
   PFConstituents_theta.clear();
   PFConstituents_ecalEnergy.clear();
   PFConstituents_hcalEnergy.clear();
   ak4pfjets_tagging.clear(); 
   for(vector<reco::PFJet>::const_iterator j = ak4pfjets->begin(); j != ak4pfjets->end(); ++j)
   {
      ak4pfjets_px[ak4pfjets_track] = j->px();
      ak4pfjets_py[ak4pfjets_track] = j->py();
      ak4pfjets_pz[ak4pfjets_track] = j->pz();
      ak4pfjets_pt[ak4pfjets_track] = j->pt();
      ak4pfjets_Et[ak4pfjets_track] = j->et();
      ak4pfjets_E[ak4pfjets_track] = j->energy();
      ak4pfjets_mass[ak4pfjets_track] = j->mass();
      ak4pfjets_mt[ak4pfjets_track] = j->mt();
      ak4pfjets_phi[ak4pfjets_track] = j->phi();
      ak4pfjets_theta[ak4pfjets_track] = j->theta();
      ak4pfjets_eta[ak4pfjets_track] = j->eta();
      ak4pfjets_rapidity[ak4pfjets_track] = j->rapidity();
      ak4pfjets_charge[ak4pfjets_track] = j->charge();
      ak4pfjets_pdgId[ak4pfjets_track] = j->pdgId();
      ak4pfjets_neutralHadronMultiplicity[ak4pfjets_track] = j->neutralHadronMultiplicity();
      ak4pfjets_chargedHadronMultiplicity[ak4pfjets_track] = j->chargedHadronMultiplicity();
      ak4pfjets_chargedMultiplicity[ak4pfjets_track] = j->chargedMultiplicity();
      ak4pfjets_neutralMultiplicity[ak4pfjets_track] = j->neutralMultiplicity();
      ak4pfjets_HFEMMultiplicity[ak4pfjets_track] = j->HFEMMultiplicity();
      ak4pfjets_HFHadronMultiplicity[ak4pfjets_track] = j->HFHadronMultiplicity();
      ak4pfjets_tagging.emplace_back(0);
      int n_ak4pfj_c = j->getPFConstituents().size();
      int Index = 0;
      for(int index = PFConstituents_track; index < PFConstituents_track+n_ak4pfj_c ; index ++)
      {
          const reco::PFCandidatePtr myptr = (j->getPFConstituent(Index));
          PFConstituents_charge.push_back(myptr->charge());
          PFConstituents_pdgId.emplace_back(myptr->pdgId());
          PFConstituents_p.emplace_back(myptr->p());
          PFConstituents_px.emplace_back(myptr->px());
          PFConstituents_py.emplace_back(myptr->py());
          PFConstituents_pz.emplace_back(myptr->pz());
          PFConstituents_pt.push_back(myptr->pt());
          PFConstituents_E.emplace_back(myptr->energy());
          PFConstituents_Et.emplace_back(myptr->et());
          PFConstituents_mass.emplace_back(myptr->mass());
          PFConstituents_mt.emplace_back(myptr->mt());
          PFConstituents_rapidity.emplace_back(myptr->rapidity());
          PFConstituents_eta.emplace_back(myptr->eta());
          PFConstituents_phi.emplace_back(myptr->phi());
          PFConstituents_theta.emplace_back(myptr->theta());
          PFConstituents_hcalEnergy.emplace_back(myptr->hcalEnergy());    
          PFConstituents_ecalEnergy.emplace_back(myptr->ecalEnergy());    
          Index++;
      }
      PFConstituents_track +=n_ak4pfj_c;
      ak4pfjets_nConstituents[ak4pfjets_track] = n_ak4pfj_c;
      ak4pfjets_tagging.at(ak4pfjets_track)=PFConstituents_track;
      ak4pfjets_track++ ;
   }
   
   hgcalclst_track = 0; 
   hgcal::RecHitTools hgcal;
   for(vector<reco::CaloCluster>::const_iterator j = hgcalclst->begin() ; j != hgcalclst->end(); ++j  )
   {
       hgcalclst_x[hgcalclst_track] = j->x();
       hgcalclst_y[hgcalclst_track] = j->y();
       hgcalclst_z[hgcalclst_track] = j->z();
       hgcalclst_E[hgcalclst_track] = j->energy();
       hgcalclst_eta[hgcalclst_track] = j->eta();
       hgcalclst_phi[hgcalclst_track] = j->phi();
       hgcalclst_size[hgcalclst_track] = j->size();
       hgcalclst_tagger[hgcalclst_track] = false;
       tagger_R[hgcalclst_track]=0.;
       if(j->z() >367.699 || j->z() < -367.699)
       {
           hgcalclst_layer[hgcalclst_track] = hgcal.getLayer(j->seed().rawId())+28;
       }
       else 
           hgcalclst_layer[hgcalclst_track] = hgcal.getLayer(j->seed().rawId());
       hgcalclst_track ++;
   }
   mytree->Fill();
   

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntupler::beginJob()
{
    // -- GenParticle --//
    using namespace reco;
    mytree->Branch("genparticle_track",&genparticle_track,"genparticle_track/I");
    mytree->Branch("genparticle_charge",genparticle_charge,"genparticle_charge[genparticle_track]/I");
    mytree->Branch("genparticle_pdgId",genparticle_pdgId,"genparticle_pdgId[genparticle_track]/I");
    mytree->Branch("genparticle_status",genparticle_status,"genparticle_status[genparticle_track]/I");
    mytree->Branch("genparticle_px",genparticle_px,"genparticle_px[genparticle_track]/D");
    mytree->Branch("genparticle_py",genparticle_py,"genparticle_py[genparticle_track]/D");
    mytree->Branch("genparticle_pz",genparticle_pz,"genparticle_pz[genparticle_track]/D");
    mytree->Branch("genparticle_pt",genparticle_pt,"genparticle_pt[genparticle_track]/D");
    mytree->Branch("genparticle_E",genparticle_E,"genparticle_E[genparticle_track]/D");
    mytree->Branch("genparticle_Et",genparticle_Et,"genparticle_Et[genparticle_track]/D");
    mytree->Branch("genparticle_mass",genparticle_mass,"genparticle_mass[genparticle_track]/D");
    mytree->Branch("genparticle_mt",genparticle_mt,"genparticle_mt[genparticle_track]/D");
    mytree->Branch("genparticle_phi",genparticle_phi,"genparticle_phi[genparticle_track]/D");
    mytree->Branch("genparticle_theta",genparticle_theta,"genparticle_theta[genparticle_track]/D");
    mytree->Branch("genparticle_eta",genparticle_eta,"genparticle_eta[genparticle_track]/D");
    mytree->Branch("genparticle_rapidity",genparticle_rapidity,"genparticle_rapidity[genparticle_track]/D");
    //-- AK4PFJet --//
    mytree->Branch("ak4pfjets_track",&ak4pfjets_track,"ak4pfjets_track/I");
    mytree->Branch("ak4pfjets_charge",ak4pfjets_charge,"ak4pfjets_charge[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_pdgId",ak4pfjets_pdgId,"ak4pfjets_pdgId[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_px",ak4pfjets_px,"ak4pfjets_px[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_py",ak4pfjets_py,"ak4pfjets_py[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_pz",ak4pfjets_pz,"ak4pfjets_pz[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_pt",ak4pfjets_pt,"ak4pfjets_pt[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_E",ak4pfjets_E,"ak4pfjets_E[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_Et",ak4pfjets_Et,"ak4pfjets_Et[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_mass",ak4pfjets_mass,"ak4pfjets_mass[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_mt",ak4pfjets_mt,"ak4pfjets_mt[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_phi",ak4pfjets_phi,"ak4pfjets_phi[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_theta",ak4pfjets_theta,"ak4pfjets_theta[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_eta",ak4pfjets_eta,"ak4pfjets_eta[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_rapidity",ak4pfjets_rapidity,"ak4pfjets_rapidity[ak4pfjets_track]/D");
    mytree->Branch("ak4pfjets_nConstituents",ak4pfjets_nConstituents,"ak4pfjets_nConstituents[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_chargedHadronMultiplicity",ak4pfjets_chargedHadronMultiplicity,"ak4pfjets_chargedHadronMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_neutralHadronMultiplicity",ak4pfjets_neutralHadronMultiplicity,"ak4pfjets_neutralHadronMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_neutralMultiplicity",ak4pfjets_neutralMultiplicity,"ak4pfjets_neutralMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_chargedMultiplicity",ak4pfjets_chargedMultiplicity,"ak4pfjets_chargedMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_HFEMMultiplicity",ak4pfjets_HFEMMultiplicity,"ak4pfjets_HFEMMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_HFHadronMultiplicity",ak4pfjets_HFHadronMultiplicity,"ak4pfjets_HFHadronMultiplicity[ak4pfjets_track]/I");
    mytree->Branch("ak4pfjets_tagging",&ak4pfjets_tagging); 
    //--Ak4PFCandidate--//
    mytree->Branch("PFConstituents_charge",&PFConstituents_charge);
    mytree->Branch("PFConstituents_pdgId",&PFConstituents_pdgId);
    mytree->Branch("PFConstituents_p",&PFConstituents_p);
    mytree->Branch("PFConstituents_px",&PFConstituents_px);
    mytree->Branch("PFConstituents_py",&PFConstituents_py);
    mytree->Branch("PFConstituents_pz",&PFConstituents_pz);
    mytree->Branch("PFConstituents_pt",&PFConstituents_pt);
    mytree->Branch("PFConstituents_E",&PFConstituents_E);
    mytree->Branch("PFConstituents_Et",&PFConstituents_Et);
    mytree->Branch("PFConstituents_mass",&PFConstituents_mass);
    mytree->Branch("PFConstituents_mt",&PFConstituents_mt);
    mytree->Branch("PFConstituents_phi",&PFConstituents_phi);
    mytree->Branch("PFConstituents_theta",&PFConstituents_theta);
    mytree->Branch("PFConstituents_eta",&PFConstituents_eta);
    mytree->Branch("PFConstituents_rapidity",&PFConstituents_rapidity);
    mytree->Branch("PFConstituents_ecalEnergy",&PFConstituents_ecalEnergy);
    mytree->Branch("PFConstituents_hcalEnergy",&PFConstituents_hcalEnergy);
 //
    //--HGCal Layer --//
    mytree->Branch("hgcalclst_track",&hgcalclst_track,"hgcalclst_track/I");
    mytree->Branch("hgcalclst_x",hgcalclst_x,"hgcalclst_x[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_y",hgcalclst_y,"hgcalclst_y[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_z",hgcalclst_z,"hgcalclst_z[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_E",hgcalclst_E,"hgcalclst_E[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_tagger",hgcalclst_tagger,"hgcalclst_tagger[hgcalclst_track]/B");
    mytree->Branch("tagger_R",tagger_R,"tagger_R[hgcalclst_track]/D"); 
    mytree->Branch("hgcalclst_eta",hgcalclst_eta,"hgcalclst_eta[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_phi",hgcalclst_phi,"hgcalclst_phi[hgcalclst_track]/D");
    mytree->Branch("hgcalclst_size",hgcalclst_size,"hgcalclst_size[hgcalclst_track]/I");
    mytree->Branch("hgcalclst_layer",hgcalclst_layer,"hgcalclst_layer[hgcalclst_track]/I");
//
//
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntupler::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntupler);
