#ifndef ECALGARLICEXTENDEDOBJECTS_HH_
#define ECALGARLICEXTENDEDOBJECTS_HH_


#include <vector>
#include <ClusterParameters.hh>
#include <HelixClass.h>

#include "EVENT/CalorimeterHit.h"

using namespace EVENT;


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
  struct GhostCluster;
  
  struct GhostHit {
    GhostHit() : Layer(0), PseudoLayer(0), Energy(0), Type(GAP_TYPE_UNKNOWN), Count(0), Stave(0) {}
    ~GhostHit() {Energy=0; Type=GAP_TYPE_UNKNOWN;}
    vec3 Position;
    int Layer;
    int PseudoLayer;
    float Energy;
    GapType Type;
    int Count;
    int Stave;
  };

  struct ExtendedHit {
    ExtendedHit() : hit(0), pseudoLayer(999), cluster(0), clusterHitVec(0), preShower(0), zone(0)  {}
    ~ExtendedHit() { hit = 0; pseudoLayer=999; cluster = 0; clusterHitVec = 0; preShower=0;zone=0;}
    CalorimeterHit* hit;
    int pseudoLayer;
    ExtendedCluster *cluster;
    std::vector<ExtendedHit* > *clusterHitVec;
    bool preShower;
    int zone;
    // for sorting by energy
    static bool higherEnergy(const ExtendedHit *a, const ExtendedHit *b)
    { return a->hit->getEnergy() > b->hit->getEnergy(); }
    // for sorting by pseudoLayer ID
    static bool lowerPseudoLayer(const ExtendedHit *a, const ExtendedHit *b)
    { return a->pseudoLayer < b->pseudoLayer; }
  };

  struct GhostCluster {
    GhostCluster() : ghostHitVec(0), ofCluster(0) {}
    ~GhostCluster() {ofCluster=0;
    int NGhostHits = ghostHitVec.size();
    if(NGhostHits>0) {
      for(unsigned int i=0;i<ghostHitVec.size();i++) {
	GhostHit *gHit = ghostHitVec[i];
	delete gHit;;
      }
      ghostHitVec.clear();
    }
    }
    std::vector<GhostHit* > ghostHitVec;
    ExtendedCluster *ofCluster;
  };


  struct Ellipsoid {
    Ellipsoid() : WidthL(0), WidthR(0), Volume(0) {}
    ~Ellipsoid() {WidthL=0; WidthR=0;Volume=0;}
    float WidthL;
    float WidthR;
    float Volume;
  };


  // extended Cluster to hold ExtendedHits + information about its origin
  struct ExtendedCluster {
    ExtendedCluster() : hitVec(0), seededFrom(), PreCluster(0), Ghosts(0), Shape(0), dir(), location(CLUS_LOCATION_UNKNOWN) {}

    ~ExtendedCluster() {
      hitVec.clear(); 
      seededFrom.x=0;
      seededFrom.y=0 ;
      seededFrom.z=0 ; 
      PreCluster=0; 
      Ghosts=0; 
      Shape=0;
      dir.x=0;
      dir.y=0;
      dir.z=0;  
      location=CLUS_LOCATION_UNKNOWN;
      if (Ghosts) delete Ghosts;
    }

    std::vector<ExtendedHit* > hitVec;
    vec3 seededFrom;
    ExtendedCluster *PreCluster;
    GhostCluster *Ghosts;
    Ellipsoid* Shape;
    vec3 dir;
    vec3 seed_dir;
    ClusterLocation location;
    ClusterParameters* Parameters;
    float rawEn;
    // for sorting by #hits
    static bool moreHits(const ExtendedCluster *a, const ExtendedCluster *b)
    { return a->hitVec.size() > b->hitVec.size() ; }
    static bool higherEnergy(const ExtendedCluster *a, const ExtendedCluster *b)
    { return (a->Parameters)->Etot_g > (b->Parameters)->Etot_g ; }
    static bool higherRawEnergy(const ExtendedCluster *a, const ExtendedCluster *b)
    { return a->rawEn > b->rawEn ; }
  };
    
  struct ExtendedTrack {
    ExtendedTrack() : helix(0), ecalEntryPoint() {}
    ~ExtendedTrack() {helix=0; ecalEntryPoint.x=0; ecalEntryPoint.y=0 ;ecalEntryPoint.z=0 ;
    delete helix;
    }
    HelixClass* helix;
    vec3 ecalEntryPoint;
  };

}


#endif
