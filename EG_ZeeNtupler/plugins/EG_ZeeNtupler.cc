// -*- C++ -*-
//
// Package:    EgammaFabrice/EG_ZeeNtupler
// Class:      EG_ZeeNtupler
// 
/**\class EG_ZeeNtupler EG_ZeeNtupler.cc EgammaFabrice/EG_ZeeNtupler/plugins/EG_ZeeNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Fabrice Couderc
//         Created:  Fri, 15 Jan 2016 10:44:20 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "EgammaFabrice/EG_ZeeNtupler/interface/EGEnergyCorrectorSemiParam.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

#include <iostream>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class EG_ZeeNtupler : public edm::one::EDAnalyzer<>  {
public:
  explicit EG_ZeeNtupler(const edm::ParameterSet&);
  ~EG_ZeeNtupler();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  // Electron, Photon 
  edm::EDGetToken electronsToken_;
  edm::EDGetToken photonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  edm::EDGetTokenT<double> rhoToken_;  

  // ID decision objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoIdMapToken_;


  
  TTree *_zeeTree;  
  int _runId;
  int _nVtx;
  
  int _eleIx[3], _eleIy[3];
  int _tracker[3];
  int _eleQ[3], _phoId[3], _isEB[3];
  double _rho;
  float _eleCombE[3];
  float _eleE[3]  , _eleEta[3]  , _elePhi[3], _phoE[3];
  float _phoE1[3], _phoE2[3], _phoE3[3], _phoE4[3];
  
  float _eleSCE[3], _eleSCEta[3], _eleSCPhi[3];
  float _eleR9[3];
  float _genE[3];
  float _eleVtxZ[3];

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void printCutFlowResult(vid::CutFlowResult &cutflow);
  EGEnergyCorrectorSemiParam _cor1,_cor2,_cor3,_cor4;

      // ----------member data ---------------------------
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
EG_ZeeNtupler::EG_ZeeNtupler(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //  usesResource("TFileService");
  edm::Service<TFileService> fs;
  _zeeTree = fs->make<TTree> ("ZeeTree", "Zee tree selection");

  /// use the same ntuple format as EcalElf to make it compatible with IJazZ
  _zeeTree->Branch("recoFlagsEle",_tracker   , "recoFlagsEle[3]/I");
  _zeeTree->Branch("energyMC"    , _genE     , "energyMC[3]/F");
  _zeeTree->Branch("energySCEle" , _eleE     , "energySCEle[3]/F");
  _zeeTree->Branch("energySCPho" , _phoE     , "energySCPho[3]/F");
  _zeeTree->Branch("etaEle"      , _eleEta   , "etaEle[3]/F");
  _zeeTree->Branch("phiEle"      , _elePhi   , "phiEle[3]/F");
  _zeeTree->Branch("ixEle"      , _eleIx     , "ixEle[3]/I");
  _zeeTree->Branch("iyEle"      , _eleIy     , "iyEle[3]/I");

  _zeeTree->Branch("zVtxEle"     , _eleVtxZ  , "eleVtxZ[3]/F");

  _zeeTree->Branch("etaSCEle"    , _eleSCEta , "etaSCEle[3]/F");
  _zeeTree->Branch("phiSCEle"    , _eleSCPhi , "phiSCEle[3]/F");
  _zeeTree->Branch("R9Ele"       , _eleR9    , "R9Ele[3]/F");
  _zeeTree->Branch("chargeEle"   , _eleQ     , "chargeEle[3]/I");
  _zeeTree->Branch("passIdPho"   , _phoId    , "passIdPho[3]/I");
  _zeeTree->Branch("isEBEle"     , _isEB     , "isEB[3]/I");
  _zeeTree->Branch("nPV"         , &_nVtx    , "nPV/I");
  _zeeTree->Branch("rho"         , &_rho     , "rho/D");
  _zeeTree->Branch("runNumber"   , &_runId   , "runNumber/I");


  _zeeTree->Branch("energyPho1" , _phoE1     , "energyPho1[3]/F");
  _zeeTree->Branch("energyPho2" , _phoE2     , "energyPho2[3]/F");
  _zeeTree->Branch("energyPho3" , _phoE3     , "energyPho3[3]/F");
  _zeeTree->Branch("energyPho4" , _phoE4     , "energyPho4[3]/F");


  // tokens
  electronsToken_  = consumes<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));
  photonsToken_    = consumes<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photons"));

  eleIdMapToken_ = consumes<edm::ValueMap<bool> >
    (iConfig.getParameter<edm::InputTag>
     ("eleId"));

  phoIdMapToken_ = consumes<edm::ValueMap<bool> >
    (iConfig.getParameter<edm::InputTag>
     ("phoId"));

  genParticlesToken_ = consumes<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));


  vtxToken_          = consumes<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  rhoToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));



  _cor1.setGBRForestFromTooFile("EgammaFabrice/EG_ZeeNtupler/data/cms_wereg_ph_bx25_AllPositions.root" , 0);
  _cor2.setGBRForestFromTooFile("EgammaFabrice/EG_ZeeNtupler/data/cms_wereg_ph_bx25_PS.root"           , 1);
  _cor3.setGBRForestFromTooFile("EgammaFabrice/EG_ZeeNtupler/data/cms_wereg_ph_bx25_noPS.root"         , 2);
  //  _cor4.setGBRForestFromTooFile("EgammaFabrice/EG_ZeeNtupler/data/cms_wereg_ph_bx25_noPS_noHoE.root"   , 3);
  _cor4.setGBRForestFromTooFile("EgammaFabrice/EG_ZeeNtupler/data/cms_wereg_ph_bx25_db.root"   , 1);

}


EG_ZeeNtupler::~EG_ZeeNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EG_ZeeNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  /// Setup the event to zero
  for( unsigned int iobj = 0 ; iobj < 3; iobj++ ) {
    _eleCombE[iobj] = -999;
    _genE [   iobj] = -999;
    _eleE [   iobj] = _phoE[ iobj ]  = _eleSCE[  iobj] = -999;
    _eleEta[  iobj] = _elePhi[iobj]  = _eleSCEta[iobj] =_eleSCPhi[iobj] = -999;
    _eleQ  [  iobj] = _tracker[iobj] = -999;
    _eleR9[   iobj] = -999;  
    _eleVtxZ[ iobj] = -999;  
    _isEB   [ iobj] = -999;
  }

  using namespace edm;

  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  _rho = *rhoH;

  // Get Electrons / Photons
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<edm::View<reco::Photon> > photons;
  iEvent.getByToken(photonsToken_, photons);

  // Get Electron / Photon VID 
  edm::Handle<edm::ValueMap<bool> > elec_id;
  iEvent.getByToken(eleIdMapToken_ ,elec_id);
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > elec_id_cutflow;
  iEvent.getByToken(eleIdMapToken_, elec_id_cutflow);
  edm::Handle<edm::ValueMap<bool> > pho_id;
  iEvent.getByToken(phoIdMapToken_ ,pho_id);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  _runId = iEvent.getRun().run();

  /// count number of good primary vertices
  _nVtx = 0;
  for( reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx ) {
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) _nVtx++;
  }

  // Get MC particles
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( !iEvent.isRealData() ) iEvent.getByToken(genParticlesToken_,genParticles);

  int iEleSaved = 0;
  // Loop over electrons
  for (size_t iele = 0; iele < electrons->size(); ++iele){
    const auto ele = electrons->ptrAt(iele);
    
    bool passEleId = (*elec_id)[ele];
    
    if( passEleId && ele->pt() > 15 && ele->superCluster().get() != 0 ) {
      _eleE [   iEleSaved] = ele->ecalEnergy();
      _eleCombE[iEleSaved] = ele->energy();
      _eleEta[  iEleSaved] = ele->eta();
      _elePhi[  iEleSaved] = ele->phi();

      _eleQ  [  iEleSaved] = ele->charge();
      _eleSCE[  iEleSaved] = ele->superCluster()->energy();
      _eleSCEta[iEleSaved] = ele->superCluster()->position().eta();
      _eleSCPhi[iEleSaved] = ele->superCluster()->position().phi();
      
      _eleR9[   iEleSaved] = ele->full5x5_r9();      
      _tracker[ iEleSaved] = ele->ecalDrivenSeed() ? 0 : 1;      
      _isEB   [ iEleSaved] = ele->isEB() ? 1 : 0;
      const edm::Ptr<reco::CaloCluster>& theseed = ele->superCluster()->seed();
      if( ele->isEB() ) {
	EBDetId ebseedid(theseed->seed());
	_eleIx[ iEleSaved] = ebseedid.ieta();
	_eleIy[ iEleSaved] = ebseedid.iphi();
      } else {
	EEDetId eeseedid(theseed->seed());
	_eleIx[ iEleSaved] = eeseedid.ix();
	_eleIy[ iEleSaved] = eeseedid.iy();
      }

      _eleVtxZ[ iEleSaved] = ele->trackPositionAtVtx().z();
      
      /// find the corresponding photon energy 
      bool phoMatched = false;
      for( size_t ipho = 0; ipho < photons->size(); ++ipho ) {
	const auto pho = photons->ptrAt(ipho);
	double dr = reco::deltaR(ele->superCluster()->position(),pho->superCluster()->position());
	//	std::cout << " dr[ele,photon " << ipho << "] = " << dr << std::endl;

	if( pho->superCluster() == ele->superCluster() || dr < 0.02 ) {
	  if( phoMatched ) std::cout << " electron was already matched to a photon... ???? ... " << std::endl;
	  else {
	    _phoE[ iEleSaved] = pho->energy();
	    //	    _phoId[iEleSaved] = (*pho_id)[pho];

	    _phoE1[iEleSaved] = _cor1.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()).first;
	    _phoE2[iEleSaved] = _cor2.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()).first;
	    _phoE3[iEleSaved] = _cor3.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()).first;
	    _phoE4[iEleSaved] = _cor4.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()).first;

	    /*
	    std::cout << " * phoE    = " << _phoE[ iEleSaved]   << " Ep4(std) = " << pho->p4().energy() << "  p4 type: " << pho->getCandidateP4type() << std::endl
	      //		      << "                        - E[ecal] = " << pho->p4(reco::Photon::ecal_photons).energy() << std::endl
	      //		      << "                        - E[reg1] = " << pho->p4(reco::Photon::regression1).energy() << std::endl
	      //		      << "                        - E[std ] = " << pho->p4(reco::Photon::ecal_standard).energy() << std::endl
		      << "      - E1 = " << _phoE1[iEleSaved] << std::endl 
		      << "      - E2 = " << _phoE2[iEleSaved] << std::endl
	      //		                         << " " << _cor2.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()-1).first 
	      //		                         << " " << _cor2.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()+1).first 
	      //                                         << " " << _cor2.CorrectedEnergyAndResolution(*pho.get(),iSetup,_rho,vertices->size()+2).first <<std::endl 
		      << "      - E3 = " << _phoE3[iEleSaved] << std::endl 
		      << "      - E4 = " << _phoE4[iEleSaved] << std::endl ;
	    */
	    phoMatched = true;
	  }
	}

      } 

      if( !phoMatched ) std::cout << " did not find any photon matching the electron" << std::endl;
      /// VID cut flow
      vid::CutFlowResult fullCutFlowData = (*elec_id_cutflow)[ele];
      std::cout << "Will you work now?" << std::endl;
      std::cout << "DEBUG CutFlow, full info for cand with pt = " <<  ele->pt() << std::endl;
      printCutFlowResult(fullCutFlowData);
      
      if( !iEvent.isRealData() ) {
	double dR = 999;
	int iMatchedGen = -1;
	/// find the corresponding generated electrons
	for( size_t igen = 0; igen < genParticles->size(); igen++) {
	  const reco::Candidate *gen = &(*genParticles)[igen];
	  if( gen->status() != 1 || abs(gen->pdgId()) != 11 ) continue;
	  double dRloc = ROOT::Math::VectorUtil::DeltaR( ele->p4(), gen->p4() );
	  if( dRloc < 0.1 && dRloc < dR ) {
	    dR = dRloc;
	    iMatchedGen = igen;
	  }
	}
	if( iMatchedGen >= 0 ) _genE[ iEleSaved ] = (&(*genParticles)[iMatchedGen])->energy();
      }

      iEleSaved++;
      if( iEleSaved >= 3 ) break;
    }
  }
  
  if( iEleSaved >= 2 ) _zeeTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
EG_ZeeNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EG_ZeeNtupler::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EG_ZeeNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


void 
EG_ZeeNtupler::printCutFlowResult(vid::CutFlowResult &cutflow) {

  printf("    CutFlow name= %s    decision is %d\n",
         cutflow.cutFlowName().c_str(),
         (int) cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
    printf("  %2d      %50s    %d        %f          %d\n", icut,
           cutflow.getNameAtIndex(icut).c_str(),
           (int)cutflow.isCutMasked(icut),
           cutflow.getValueCutUpon(icut),
           (int)cutflow.getCutResultByIndex(icut));
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(EG_ZeeNtupler);
