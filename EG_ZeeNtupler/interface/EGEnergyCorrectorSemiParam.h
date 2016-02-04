//--------------------------------------------------------------------------------------------------
// $Id $
//
// SCEnergyCorrectorSemiParm
//
// Helper Class for applying regression-based energy corrections with optimized BDT implementation
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef EGAMMATOOLS_SCEnergyCorrectorSemiParm_H
#define EGAMMATOOLS_SCEnergyCorrectorSemiParm_H
    
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
//#include "Geometry/CaloTopology/interface/CaloTopology.h" 
//#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
class EGEnergyCorrectorSemiParam {
 public:
  EGEnergyCorrectorSemiParam();
  ~EGEnergyCorrectorSemiParam();

  void setGBRForestFromTooFile( const std::string rootfilename, int config = 0 );
  std::pair<double,double> CorrectedEnergyAndResolution(const reco::Photon & pho, 
							const edm::EventSetup& es,
							double rho, int nVtx);

 
  
 protected:    
  const GBRForestD *foresteb_;
  const GBRForestD *forestee_;
  const GBRForestD *forestsigmaeb_;
  const GBRForestD *forestsigmaee_;

  int _forestConfig;
};
#endif
