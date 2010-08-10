#ifndef ECALPRECLUSTERING_HH_
#define ECALPRECLUSTERING_HH_


#include <cmath>
#include <string>
#include <vector>

#include <lcio.h>

#include <marlin/Processor.h>

using namespace lcio;
using namespace std;

class lcio::CalorimeterHit;

class ECALPreClustering : public marlin::Processor {


 public:

  virtual marlin::Processor * newProcessor() { return new ECALPreClustering; }

  ECALPreClustering();

  virtual void init();
  virtual void processRunHeader(LCRunHeader * run);
  virtual void processEvent(LCEvent * evt);
  virtual void check(LCEvent * evt);
  virtual void end();

 private:

  // simple extended CalorimeterHit that knows if it is clustered

  struct ExtendedHit {

    ExtendedHit() : hit(0), isClustered(0) {}
    ~ExtendedHit() { hit = 0; isClustered = 0; }

    CalorimeterHit *hit;
    bool isClustered;

    static bool higherEnergy(const ExtendedHit *a, const ExtendedHit *b)
    { return ((a->hit)->getEnergy() > (b->hit)->getEnergy()); }


  };

  
  // algorithm functions
  float GetDistance(float *a,float *b);
  void PrepareHits(const LCEvent *evt, vector<ExtendedHit* > &hitVec);
  void PreCluster( vector<ExtendedHit* > &hitVec, multimap<int,ClusterImpl*>  &preClusVec);
  void ClusterNearHits(ExtendedHit *s_hit, vector<ExtendedHit* > &hits, ClusterImpl &cluster);
  void WritePreClusters(LCEvent *evt, multimap<int,ClusterImpl* > &preClusVec);


  // Collection names
  string _ecalEndcapHitCollectionName;
  string _ecalBarrelHitCollectionName;
  string _ecalOtherHitCollectionName;
  string _ecalPreClusterCollectionName;
 

  // some parameters
  int _nEvents;
  int _nPreClusters;
  int _nPreClustersInEvent;
  float _distanceCut;
  int _minHits;
  int _debug;
};

#endif
