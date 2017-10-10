// -*- C++ -*-
//
// Package:    Analyzer/GMSBGenAnalyzer
// Class:      GMSBGenAnalyzer
// 
/**\class GMSBGenAnalyzer GMSBGenAnalyzer.cc Analyzer/GMSBGenAnalyzer/plugins/GMSBGenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Mon, 09 Oct 2017 11:36:42 GMT
//
//


// system include files
#include <memory>

#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/src/one/implementorsMethods.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
// ROOT
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//
using namespace std;
using namespace edm;
namespace ECAL
{

  static const Float_t etaEB = 1.4442;
  static const Float_t etaEEmin = 1.566;
  static const Float_t etaEEmax = 2.5;

  static const Float_t rEB = 129.f; // 1.29 m
  float inetaEB=tan(2.f*atan(std::exp(-etaEB)));
  static const Float_t zEB = rEB / inetaEB;
  float inetaEEmin=tan(2.f*atan(std::exp(-etaEEmin)));
  float inetaEEmax=tan(2.f*atan(std::exp(-etaEEmax)));
  static const Float_t zEE = 314.f; // 3.14 m 
  static const Float_t rEEmin = zEE * inetaEEmin;
  static const Float_t rEEmax = zEE * inetaEEmax;
};
static const Float_t sol = 2.99792458e10; // cm/s
struct Tree_struc_ {
  int run;
  int lumi;
  int event;
  float genVtxX;
  float genVtxY;
  float genVtxZ;
  float neutM;
  float neutE;
  float neutPt;
  float neutPhi;
  float neutEta;
  float neutProdVtxX;
  float neutProdVtxY;
  float neutProdVtxZ;
  float neutDecayVtxX;
  float neutDecayVtxY;
  float neutDecayVtxZ;
  float neutLabDist;
  float neutDecayR;
  float neutP;
  float neutBetaGamma;
  float neutGamma;
  float neutCTau;
  float phoE;
  float phoPt;
  float phoPhi;
  float phoEta;
  float phoTime;
};

class GMSBGenAnalyzer : public edm::EDAnalyzer {
public:
  explicit GMSBGenAnalyzer(const edm::ParameterSet&);
  ~GMSBGenAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void initTreeStructure(); 

  inline Float_t rad2 (const Float_t x, const Float_t y);
  float GetGenPhotonArrivalTime(const float r0, const float z0, const float slope);
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  //  edm::EDGetTokenT<View<reco::GenParticle> > genPartToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genPartToken_;
  // event info
  unsigned long int event;
  unsigned int run, lumi;
  float   genVtxX,genVtxY,genVtxZ;
  float neutM;
  float neutE;
  float neutPt;
  float neutPhi;
  float neutEta;
  float neutProdVtxX;
  float neutProdVtxY;
  float neutProdVtxZ;
  float neutDecayVtxX;
  float neutDecayVtxY;
  float neutDecayVtxZ;
  float neutLabDist;
  float neutDecayR;
  float neutP;
  float neutBetaGamma;
  float neutGamma;
  float neutCTau;
  float phoE;
  float phoPt;
  float phoPhi;
  float phoEta;
  float phoTime;
  // output event level ntuple
  TTree* tree; 
  Tree_struc_ tree_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GMSBGenAnalyzer::GMSBGenAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  genPartToken_ = consumes<std::vector<reco::GenParticle>> (iConfig.getUntrackedParameter<InputTag> ("GenParticlesTag", InputTag("genParticles"))) ;
  //  genpartsToken   = consumes<std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genparts"));

}


GMSBGenAnalyzer::~GMSBGenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GMSBGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //   edm::Handle<View<reco::GenParticle> genparticlesH;
   //   edm::Handle<std::vector<reco::GenParticle> > genparticlesH;
   //iEvent.getByToken(genpartsToken,   genparticlesH);

   edm::Handle<std::vector<reco::GenParticle> > genParticlesH;
   iEvent.getByToken( genPartToken_, genParticlesH );

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

   
   int   run   = iEvent.id().run();
   int lumi  = iEvent.luminosityBlock();
   int  event = iEvent.id().event();

   tree_.run=run;
   tree_.lumi=lumi;
   tree_.event=event;

   float  genVtxX=999.;
   float  genVtxY=999.;
   float  genVtxZ=999.;
   float neutM=999.;
   float neutE=999.;
   float neutPt=999.;
   float neutPhi=999.;
   float neutEta=999.;
   float  neutProdVtxX=999.;
   float neutProdVtxY=999.;
   float neutProdVtxZ=999.;
   float neutDecayVtxX=999.;
   float neutDecayVtxY=999.;
   float neutDecayVtxZ=999.;
   float neutLabDist=999.;
   float neutDecayR=999.;
   float neutP  = 999.;
   float neutBetaGamma=999.;
   float neutGamma=999.;
   float neutCTau=999.;

   float phoE=999.;
   float phoPt=999.;
   float phoPhi=999.;
   float phoEta=999.;
   float phoTime=999.;
   //   const auto & genParticles = *genparticlesH;
   //   for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
     for (const auto & genpart : *genParticlesH) // loop over gen particles
       {
	 //interaction vertex
	 if( genpart.pdgId() != 2212 || genpart.vertex().z() != 0. ) { // pdg1d=2212 is proton vtx
	   genVtxX = genpart.vertex().x();
	   genVtxY = genpart.vertex().y();
	   genVtxZ = genpart.vertex().z();
	 }

	 //neutralino production vertex
	 if (genpart.pdgId() == 1000022 && genpart.numberOfDaughters() == 2){
	 
	   neutM = genpart.mass();
	   neutE    = genpart.energy();
	   neutPt   = genpart.pt();
	   neutPhi  = genpart.phi();
	   neutEta  = genpart.eta();

	   // neutralino production vertex
	   neutProdVtxX = genpart.vx();
	   neutProdVtxY = genpart.vy();
	   neutProdVtxZ = genpart.vz();
  
	   // neutralino decay vertex (same for both daughters unless really screwed up)
	   neutDecayVtxX = genpart.daughter(0)->vx();
	   neutDecayVtxY = genpart.daughter(0)->vy();
	   neutDecayVtxZ = genpart.daughter(0)->vz();

	   //computing neutralino kinematics
	   neutLabDist=sqrt(pow((neutDecayVtxX-neutProdVtxX),2)+pow((neutDecayVtxY-neutProdVtxY),2)+pow((neutDecayVtxZ-neutProdVtxZ),2));
	   neutDecayR =sqrt(pow(neutDecayVtxX,2)+pow(neutDecayVtxY,2));
	   TLorentzVector neut_lorvec; 
	   neut_lorvec.SetPtEtaPhiE(neutPt,neutEta,neutPhi,neutE);
	   neutP  = sqrt(pow(neut_lorvec.Px(),2)+pow(neut_lorvec.Py(),2)+pow(neut_lorvec.Pz(),2));
	   neutBetaGamma=neutP/neutM;
	   neutGamma=neutE/neutM;
	   neutCTau=neutLabDist/neutBetaGamma;


	   //neutralino decay to photon +gravitino
	   if ((genpart.daughter(0)->pdgId() == 22 && genpart.daughter(1)->pdgId() == 1000039) ||
	       (genpart.daughter(1)->pdgId() == 22 && genpart.daughter(0)->pdgId() == 1000039)) {

	     const int phdaughter = (genpart.daughter(0)->pdgId() == 22)?0:1;
	         
	     phoE   = genpart.daughter(phdaughter)->energy();
	     phoPt  = genpart.daughter(phdaughter)->pt();
	     phoPhi = genpart.daughter(phdaughter)->phi();
	     phoEta = genpart.daughter(phdaughter)->eta();


	     
	     if (neutDecayR < 129. && neutDecayVtxZ < 314.) //radius within EB and z within EE 
	       {
		 float ineta=tan(2.f*atan(std::exp(-phoEta)));
		 phoTime =  GetGenPhotonArrivalTime(neutDecayR,neutDecayVtxZ,ineta);
		 phoTime = ((phoTime != -1.f) ? phoTime + neutCTau * neutGamma / sol : -1.f)*1000000000;//time in ns
	       }
	   
	   }
	 
	 }
       }
       
     tree_.genVtxX=genVtxX;
     tree_.genVtxY=genVtxY;
     tree_.genVtxZ=genVtxZ;
     tree_.neutM=neutM;
     tree_.neutE=neutE;
     tree_.neutPt=neutPt;
     tree_.neutPhi=neutPhi;
     tree_.neutEta=neutEta;
     tree_.neutProdVtxX=neutProdVtxX;
     tree_.neutProdVtxY=neutProdVtxY;
     tree_.neutProdVtxZ=neutProdVtxZ;
     tree_.neutDecayVtxX=neutDecayVtxX;
     tree_.neutDecayVtxY=neutDecayVtxY;
     tree_.neutDecayVtxZ=neutDecayVtxZ;
     tree_.neutLabDist=neutLabDist;     
     tree_.neutDecayR=neutDecayR;     
     tree_.neutP=neutP;     
     tree_.neutBetaGamma=neutBetaGamma;     
     tree_.neutGamma=neutGamma;     
     tree_.neutCTau=neutCTau;     
     tree_.phoE=phoE;
     tree_.phoPt=phoPt;
     tree_.phoPhi=phoPhi;
     tree_.phoEta=phoEta;
     tree_.phoTime=phoTime;

     
     tree->Fill();
          


}

inline Float_t GMSBGenAnalyzer::rad2 (const Float_t x, const Float_t y){return x*x + y*y;}
     
float GMSBGenAnalyzer::GetGenPhotonArrivalTime(const float r0, const float z0, const float slope)     {
  Float_t time = -1.f;
  
  const Float_t z = z0 + ( (ECAL::rEB - r0) / slope );
  if (std::abs(z) < ECAL::zEB) 
    {
      time  = std::sqrt(rad2(ECAL::rEB-r0,z-z0)) / sol;
      
      std::cout << "EB: " << time << " (z" << z << ")" ;
      
      time -= std::sqrt(rad2(ECAL::rEB,z)) / sol; // TOF correction;
	   
	   std::cout << " w TOF: " << time;
    }
  else 
    {
      const Float_t r = r0 + ( slope * (ECAL::zEE - z0) );
      if (std::abs(r) > ECAL::rEEmin && std::abs(r) < ECAL::rEEmax)
	{
	  time  = std::sqrt(rad2(r-r0,ECAL::zEE-z0)) / sol;
	  
	  std::cout << "EE: " << time << " (z" << r << ")" ;
	  
	  time -= std::sqrt(rad2(r0,ECAL::zEE)) / sol; // TOF correction;
	  
	  std::cout << " w TOF: " << time;
	}
	 }
  
  return time;
     }


// ------------ method called once each job just before starting event loop  ------------
void GMSBGenAnalyzer::beginJob()
{
  std::cout<<"begin job" <<std::endl;
  tree = fs->make<TTree>("tree","tree");

  tree->Branch("run", &tree_.run, "run/i");
  tree->Branch("lumi", &tree_.lumi, "lumi/i");
  tree->Branch("event", &tree_.event, "event/l");
  tree->Branch("genVtxX", &tree_.genVtxX, "genVtxX/f");
  tree->Branch("genVtxY", &tree_.genVtxY, "genVtxY/f");
  tree->Branch("genVtxZ", &tree_.genVtxZ, "genVtxZ/f");
  tree->Branch("neutM", &tree_.neutM, "neutM/f");
  tree->Branch("neutE", &tree_.neutE, "neutE/f");
  tree->Branch("neutPt", &tree_.neutPt, "neutPt/f");
  tree->Branch("neutPhi", &tree_.neutPhi, "neutPhi/f");
  tree->Branch("neutEta", &tree_.neutEta, "neutEta/f");
  tree->Branch("neutProdVtxX", &tree_.neutProdVtxX, "neutProdVtxX/f");
  tree->Branch("neutProdVtxY", &tree_.neutProdVtxY, "neutProdVtxY/f");
  tree->Branch("neutProdVtxZ", &tree_.neutProdVtxZ, "neutProdVtxZ/f");
  tree->Branch("neutDecayVtxX", &tree_.neutDecayVtxX, "neutDecayVtxX/f");
  tree->Branch("neutDecayVtxY", &tree_.neutDecayVtxY, "neutDecayVtxY/f");
  tree->Branch("neutDecayVtxZ", &tree_.neutDecayVtxZ, "neutDecayVtxZ/f");
  tree->Branch("neutLabDist", &tree_.neutLabDist, "neutLabDist/f");
  tree->Branch("neutP", &tree_.neutP, "neutP/f");
  tree->Branch("neutDecayR", &tree_.neutDecayR, "neutDecayR/f");
  tree->Branch("neutBetaGamma", &tree_.neutBetaGamma, "neutBetaGamma/f");
  tree->Branch("neutGamma", &tree_.neutGamma, "neutGamma/f");
  tree->Branch("neutCTau", &tree_.neutCTau, "neutCTau/f");
  tree->Branch("phoE", &tree_.phoE, "phoE/f");
  tree->Branch("phoPt", &tree_.phoPt, "phoPt/f");
  tree->Branch("phoPhi", &tree_.phoPhi, "phoPhi/f");
  tree->Branch("phoEta", &tree_.phoEta, "phoEta/f");
  tree->Branch("phoTime", &tree_.phoTime, "phoTime/f");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GMSBGenAnalyzer::endJob() 
{
}


void GMSBGenAnalyzer::initTreeStructure() {
  tree_.run = 0.;

}
// ------------ method called when starting to processes a run  ------------
/*
void 
GMSBGenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GMSBGenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GMSBGenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GMSBGenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GMSBGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GMSBGenAnalyzer);
