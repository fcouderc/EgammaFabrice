#include "EgammaWork/EG_ZeeNtupler/interface/EGEnergyCorrectorSemiParam.h"

#include "CondFormats/DataRecord/interface/GBRDWrapperRcd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"

#include <vdt/vdtMath.h>

#include "TFile.h"
#include <iostream>
using namespace reco;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParam::EGEnergyCorrectorSemiParam() :
  foresteb_(0),
  forestee_(0),
  forestsigmaeb_(0),
  forestsigmaee_(0),
  _forestConfig(0)
{}

//--------------------------------------------------------------------------------------------------
EGEnergyCorrectorSemiParam::~EGEnergyCorrectorSemiParam()
{
  if( foresteb_ != 0 ) delete foresteb_;
  if( forestee_ != 0 ) delete forestee_;
  if( forestsigmaeb_ != 0 ) delete forestsigmaeb_;
  if( forestsigmaee_ != 0 ) delete forestsigmaee_;


}


void EGEnergyCorrectorSemiParam::setGBRForestFromTooFile( const std::string rootfilename, int config ) {
  edm::FileInPath infile(rootfilename.c_str());
  
  //load forests from file
  TFile *fgbr = TFile::Open(infile.fullPath().c_str(),"READ");    
  fgbr->GetObject("EBCorrection", foresteb_);
  fgbr->GetObject("EECorrection", forestee_);
  //  fgbr->GetObject("EBUncertainty", forestsigmaeb_);
  //  fgbr->GetObject("EEUncertainty", forestsigmaee_);
  fgbr->Close();   

  _forestConfig = config;


}


 


std::pair<double,double> EGEnergyCorrectorSemiParam::CorrectedEnergyAndResolution(const reco::Photon & pho, 
										  const edm::EventSetup& es, double rho, int nVtx) {

  std::array<float, 39> eval;
  const reco::SuperClusterRef& the_sc = pho.superCluster();
  const edm::Ptr<reco::CaloCluster>& theseed = the_sc->seed();

  const int numberOfClusters =  the_sc->clusters().size();
  const bool missing_clusters = !the_sc->clusters()[numberOfClusters-1].isAvailable();

  if( missing_clusters ) return  std::pair<double,double>(-999,-999); // do not apply corrections in case of missing info (slimmed MiniAOD electrons)

  const double raw_energy = the_sc->rawEnergy();
  const auto& ess = pho.showerShapeVariables();

  int ivar = 0;

  // SET INPUTS
  eval[ivar++]  = raw_energy;
  if( _forestConfig == 0 ) eval[ivar++]  = the_sc->position().eta();
  if( _forestConfig == 0 ) eval[ivar++]  = the_sc->position().phi();
  eval[ivar++]  = pho.r9();
  eval[ivar++]  = the_sc->etaWidth();
  eval[ivar++]  = the_sc->phiWidth();
  eval[ivar++]  = std::max(0,numberOfClusters - 1);
  if( _forestConfig != 3 ) eval[ivar++]  = pho.hadronicOverEm();
  eval[ivar++]  = rho;
  eval[ivar++]  = nVtx;
  eval[ivar++] = theseed->eta()-the_sc->position().Eta();
  eval[ivar++] = reco::deltaPhi(theseed->phi(),the_sc->position().Phi());
  eval[ivar++] = theseed->energy()/raw_energy;
  eval[ivar++] = ess.e3x3/ess.e5x5;
  eval[ivar++] = ess.sigmaIetaIeta;
  eval[ivar++] = ess.sigmaIphiIphi;
  eval[ivar++] = ess.sigmaIetaIphi/(ess.sigmaIphiIphi*ess.sigmaIetaIeta);
  eval[ivar++] = ess.maxEnergyXtal/ess.e5x5;
  eval[ivar++] = ess.e2nd/ess.e5x5;
  eval[ivar++] = ess.eTop/ess.e5x5;
  eval[ivar++] = ess.eBottom/ess.e5x5;
  eval[ivar++] = ess.eLeft/ess.e5x5;
  eval[ivar++] = ess.eRight/ess.e5x5;
  eval[ivar++] = ess.e2x5Max/ess.e5x5;
  eval[ivar++] = ess.e2x5Left/ess.e5x5;
  eval[ivar++] = ess.e2x5Right/ess.e5x5;
  eval[ivar++] = ess.e2x5Top/ess.e5x5;
  eval[ivar++] = ess.e2x5Bottom/ess.e5x5;

  const bool iseb = pho.isEB();
  EcalClusterLocal ecalLocal;
  if (iseb) {
    eval[ivar++] = pho.e5x5()/theseed->energy();
    //    EBDetId ebseedid(theseed->seed());
    //    eval[ivar++] = ebseedid.ieta();
    //    eval[ivar++] = ebseedid.iphi();
    
    float cryPhi, cryEta, thetatilt, phitilt;
    int ieta, iphi;
    ecalLocal.localCoordsEB(*(theseed),es,cryEta,cryPhi, ieta, iphi, thetatilt, phitilt);
    eval[ivar++] = ieta;
    eval[ivar++] = iphi;

    eval[ivar++] = (ieta-1*abs(ieta)/ieta)%5;
    eval[ivar++] = (iphi-1)%2;
    eval[ivar++] = (abs(ieta)<=25)*((ieta-1*abs(ieta)/ieta)%25) + (abs(ieta)>25)*((ieta-26*abs(ieta)/ieta)%20);
    eval[ivar++] = (iphi-1)%20;
    if( _forestConfig == 0 ) eval[ivar++] = cryPhi;
    if( _forestConfig == 0 ) eval[ivar++] = cryEta;
    eval[ivar++] = ieta;
    eval[ivar++] = iphi;
    
  } else {
    if( _forestConfig != 2 && _forestConfig != 3 ) {
      eval[ivar++] = the_sc->preshowerEnergy()/raw_energy;
      eval[ivar++] = the_sc->preshowerEnergyPlane1()/raw_energy;
      eval[ivar++] = the_sc->preshowerEnergyPlane2()/raw_energy;
    }
    //    EEDetId eeseedid(theseed->seed());
    //    eval[ivar++] = eeseedid.ix();
    //    eval[ivar++] = eeseedid.iy();
    
    float cryX, cryY, thetatilt, phitilt;
    int ix, iy;
    ecalLocal.localCoordsEE(*(theseed), es, cryX, cryY, ix, iy, thetatilt, phitilt);
    if( _forestConfig == 0 ) eval[ivar++] = cryX;
    if( _forestConfig == 0 ) eval[ivar++] = cryY;
    eval[ivar++] = ix;
    eval[ivar++] = iy;
   
  }
  //  std::string ecalStr[2] = {"ee", "eb"};
  //  std::cout << "  config: " << _forestConfig << " #vars: " << ivar << " Ecal: " << ecalStr[int(iseb)]  << std::endl;

  //magic numbers for MINUIT-like transformation of BDT output onto limited range
  //(These should be stored inside the conditions object in the future as well)
  const double meanlimlow  = 0.2;
  const double meanlimhigh = 2.0;
  const double meanoffset  = meanlimlow + 0.5*(meanlimhigh-meanlimlow);
  const double meanscale   = 0.5*(meanlimhigh-meanlimlow);

  const double sigmalimlow  = 0.0002;
  const double sigmalimhigh = 0.5;
  const double sigmaoffset  = sigmalimlow + 0.5*(sigmalimhigh-sigmalimlow);
  const double sigmascale   = 0.5*(sigmalimhigh-sigmalimlow);

  //these are the actual BDT responses
  const GBRForestD *cor_ph_forest = 0;
  const GBRForestD *sig_ph_forest = 0;
  if( iseb ) {
    cor_ph_forest = foresteb_;
    sig_ph_forest = forestsigmaeb_;
  } else {
    cor_ph_forest = forestee_;
    sig_ph_forest = forestsigmaee_;
  }
  //these are the actual BDT responses
  double rawcor = cor_ph_forest->GetResponse(eval.data());
  double rawsig = -1;//sig_ph_forest->GetResponse(eval.data());
  //apply transformation to limited output range (matching the training)
  double cor = meanoffset  + meanscale *vdt::fast_sin(rawcor);
  double sig = sigmaoffset + sigmascale*vdt::fast_sin(rawsig);

  double      ecor = cor*eval[0];
  if( !iseb && 
      _forestConfig != 2 && 
      _forestConfig != 3 )  ecor = cor*(eval[0]+the_sc->preshowerEnergy());

  double sigmacor = sig*ecor;
  
  return std::pair<double,double>(ecor,sigmacor);

}


