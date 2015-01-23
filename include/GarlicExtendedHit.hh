#ifndef ECALGARLICEXTENDEDHIT_HH_
#define ECALGARLICEXTENDEDHIT_HH_

#include <vector>
#include <HelixClass.h>
#include <UTIL/CellIDDecoder.h>

#include "GarlicGeometryParameters.hh"

#include "EVENT/CalorimeterHit.h"

#include "GarlicExtendedObjects.hh"
#include "GarlicGeometryHelpers.hh"

using namespace EVENT;

using namespace GarlicExtendedObjects;
using std::vector;

class GarlicExtendedCluster;

class PointInCalo : public GarlicGeometryHelpers {

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
    if (_pseudoLayer<-998) {
      //      cout << "getting psLay..." << endl;
      assignPseudoLayer();
    }
    return _pseudoLayer;
  }

  int getZone(bool check=false) {
    if (_zone<-998 || check) assignZone(check);
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
  void assignZone(bool check=false);
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


class GarlicExtendedHit : public PointInCalo {

public:

  GarlicExtendedHit(CalorimeterHit* hit) :  PointInCalo( hit->getPosition() ) {    
    _hit=hit;
    _cluster=NULL;
    _preShower=false;
    _energyMIP=-999;
    _layer=-999;
    _stave=-999;
    _module=-999;
    _thick0=-999;
  }

  ~GarlicExtendedHit() {
  }

  enum ClusterLocation2 {
    CLUS2_LOCATION_UNKNOWN  = 0,
    CLUS2_LOCATION_BARREL = 1,
    CLUS2_LOCATION_ENDCAP = 2,
    CLUS2_LOCATION_OVERLAP = 3
  };

  CalorimeterHit* getCaloHit() {return _hit;}

  GarlicExtendedCluster* getCluster() {return _cluster;}
  bool isPreshower() {return _preShower;}

  void setCluster(GarlicExtendedCluster* cl) {_cluster=cl;}
  void setPreshower(bool isPS=true) {_preShower=isPS;}

  // for sorting by energy
  static bool higherEnergy(const GarlicExtendedHit *a, const GarlicExtendedHit *b) { 
    return a->_hit->getEnergy() > b->_hit->getEnergy();
  }

  // for sorting by pseudoLayer ID
  static bool lowerPseudoLayer(const GarlicExtendedHit *a, const GarlicExtendedHit *b) { 
    return a->_pseudoLayer < b->_pseudoLayer;
  }

  int getLayer() { getRegion(); return _layer; }
  int getStave() { getRegion(); return _stave; }
  int getModule() { getRegion(); return _module; }

  float getEnergyMIP();

private:

  int _layer;
  int _stave;
  int _module;

  void getRegion();

  CalorimeterHit* _hit;
  GarlicExtendedCluster* _cluster;
  bool _preShower;
  float _energyMIP;
  float _thick0;

};
    
#endif
