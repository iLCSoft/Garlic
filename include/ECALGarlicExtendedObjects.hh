#ifndef ECALGARLICEXTENDEDOBJECTS_HH_
#define ECALGARLICEXTENDEDOBJECTS_HH_


#include <vector>

#include "EVENT/CalorimeterHit.h"

using namespace EVENT;

class ExtendedCluster2;

namespace ECALGarlicExtendedObjects {

  enum ClusterLocation {
    CLUS_LOCATION_UNKNOWN  = 0,
    CLUS_LOCATION_BARREL = 1,
    CLUS_LOCATION_ENDCAP = 2,
    CLUS_LOCATION_OVERLAP = 3
  };

  enum CalorimeterHitZone {
    CALHITZONE_UNKNOWN =0,
    CALHITZONE_BARREL = 1,
    CALHITZONE_ENDCAP = 2,
    CALHITZONE_OVERLAP = 3,
    CALHITZONE_RING = 4
  };

  enum GapType {
    GAP_TYPE_UNKNOWN = 0,
    GAP_TYPE_WAFER = 1,
    GAP_TYPE_ALVEOLA = 2,
    GAP_TYPE_MODULE = 3,
    GAP_TYPE_CROSS = 4
  };

  struct vec3 {
    float x;
    float y;
    float z;
  };

  struct ExtendedCluster;
  //  struct GhostCluster;
  
  //  struct GhostHit {
  //    GhostHit() : Layer(0), PseudoLayer(0), Energy(0), Type(GAP_TYPE_UNKNOWN), Count(0), Stave(0) {}
  //    ~GhostHit() {Energy=0; Type=GAP_TYPE_UNKNOWN;}
  //    vec3 Position;
  //    int Layer;
  //    int PseudoLayer;
  //    float Energy;
  //    GapType Type;
  //    int Count;
  //    int Stave;
  //  };
  //

  class GhostHit {

  public:
    GhostHit() {
      for (int i=0; i<3; i++) _Position[i]=-99999;
      _Layer=-999;
      _PseudoLayer=-999;
      _Energy=-999;
      _Type=GAP_TYPE_UNKNOWN;
      _Count=-999;
      _Stave=-999;
    }
    ~GhostHit() {};


    void set_Position   (float* pos) {for (int i=0; i<3; i++) _Position[i]=pos[i];}
    void set_Layer      (int la) {_Layer=la;}
    void set_PseudoLayer(int la) {_PseudoLayer=la;}
    void set_Energy     (float fl) {_Energy=fl;}
    void set_Type       (GapType gt) {_Type=gt;}
    void set_Count      (int la) {_Count=la;}
    void set_Stave      (int la) {_Stave=la;}

    float*  get_Position() {return _Position;}
    int     get_Layer() {return _Layer;}
    int     get_PseudoLayer() {return _PseudoLayer;}
    float   get_Energy() {return _Energy;}
    GapType get_Type() {return _Type;}
    int     get_Count() {return _Count;}
    int     get_Stave () {return _Stave;}
    
  private:
    float   _Position[3];
    int     _Layer;
    int     _PseudoLayer;
    float   _Energy;
    GapType _Type;
    int     _Count;
    int     _Stave;


  };

  class GhostCluster {

  public:
    GhostCluster() {
      _ofCluster=NULL;
    }
    ~GhostCluster() {
      if (_ghostHitVec.size()>0) {
	for (size_t i=0; i<_ghostHitVec.size(); i++) {
	  delete _ghostHitVec[i];
	}
	_ghostHitVec.clear();
      }
    }
    
    const std::vector<GhostHit* >* getGhostHits() {return &_ghostHitVec;}
    const ExtendedCluster2* getParentCluster() {return _ofCluster;}
    
    void setGhostHits(std::vector<GhostHit* >* ghv) {_ghostHitVec = *ghv;}
    void setParentCluster(ExtendedCluster2* cl) {_ofCluster=cl;}

  private:
    std::vector<GhostHit* > _ghostHitVec;
    ExtendedCluster2* _ofCluster;
    
  };

  struct Ellipsoid {
    Ellipsoid() : WidthL(0), WidthR(0), Volume(0) {}
    ~Ellipsoid() {WidthL=0; WidthR=0;Volume=0;}
    float WidthL;
    float WidthR;
    float Volume;
  };

}

#endif
