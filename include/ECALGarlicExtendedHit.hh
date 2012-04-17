#ifndef ECALGARLICEXTENDEDHIT_HH_
#define ECALGARLICEXTENDEDHIT_HH_

#include <vector>
#include <HelixClass.h>
#include <UTIL/CellIDDecoder.h>

#include "ECALGarlicGeometryParameters.hh"

#include "EVENT/CalorimeterHit.h"

#include "ECALGarlicExtendedObjects.hh"
#include "ECALGarlicGeometryHelpers.hh"

using namespace EVENT;

using namespace ECALGarlicExtendedObjects;
using std::vector;

class ExtendedCluster2;

class PointInCalo : public ECALGarlicGeometryHelpers {

public:

  PointInCalo(const float* xyz) {
    _zone=-999;
    _pseudoLayer=-999;
    _pseudoStave=-999;
    _x0depth=-999;
    _pos = xyz;
  }    

  ~PointInCalo() {}

  int getPseudoLayer() {
    if (_pseudoLayer<-998) assignPseudoLayer();
    return _pseudoLayer;
  }

  int getZone() {
    if (_zone<-998) assignZone();
    return _zone;
  }

  float getX0depth() {
    if (_x0depth<-998) calculateX0depth();
    return _x0depth;
  }

  int getPseudoStave() {
    if (_pseudoStave<-998) assignPseudoStave();
    return _pseudoStave;
  }

  const float* getPosition() {
    return _pos;
  }

protected:
  int _pseudoLayer;
  int _zone;
  float _x0depth;
  const float* _pos;
  int _pseudoStave;

  void assignPseudoStave();
  void assignPseudoLayer();
  void assignZone();
  void calculateX0depth();
  float getAbsThickness(int layer);

  float getBarrelCorrAbsThickness(int layer, const float* pos);
  float getBarrelCorrAbsThickness() {
    if (_pseudoLayer<0) assignPseudoLayer();
    return getBarrelCorrAbsThickness(_pseudoLayer, _pos);
  }

  int getEndcapPseudoLayer(const float* pos);
  int getEndcapPseudoLayer() {return getEndcapPseudoLayer(_pos);}
  int getBarrelPseudoLayer(const float* pos);
  int getBarrelPseudoLayer() {return getBarrelPseudoLayer(_pos);}

  void movePseudoLayerAlongDirection(int dps, const float* dir, float* newpos);

};


class ExtendedHit2 : public PointInCalo {

public:

  ExtendedHit2(CalorimeterHit* hit) :  PointInCalo( hit->getPosition() ) {    
    _hit=hit;
    _cluster=NULL;
    _preShower=false;

    // int K = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (_hit)["K-1"];
    // if (this->getZone()==CALHITZONE_ENDCAP && K<5) {
    //   int zon =  this->getZone();
    //   int pl =  this->getPseudoLayer();      
    //   std::cout << "rabbits " << zon << " " << _pos[2] << " " << K << " " << pl << " " << pl-K << std::endl;
    // }

  }

  ~ExtendedHit2() {
  }

  enum ClusterLocation2 {
    CLUS2_LOCATION_UNKNOWN  = 0,
    CLUS2_LOCATION_BARREL = 1,
    CLUS2_LOCATION_ENDCAP = 2,
    CLUS2_LOCATION_OVERLAP = 3
  };

  CalorimeterHit* getCaloHit() {return _hit;}

  ExtendedCluster2* getCluster() {return _cluster;}
  bool isPreshower() {return _preShower;}
  void setCluster(ExtendedCluster2* cl) {_cluster=cl;}
  void setPreshower(bool isPS=true) {_preShower=isPS;}

  // for sorting by energy
  static bool higherEnergy(const ExtendedHit2 *a, const ExtendedHit2 *b) { 
    return a->_hit->getEnergy() > b->_hit->getEnergy();
  }

  // for sorting by pseudoLayer ID
  static bool lowerPseudoLayer(const ExtendedHit2 *a, const ExtendedHit2 *b) { 
    return a->_pseudoLayer < b->_pseudoLayer;
  }

private:

  CalorimeterHit* _hit;
  ExtendedCluster2* _cluster;
  bool _preShower;

};
    
#endif
