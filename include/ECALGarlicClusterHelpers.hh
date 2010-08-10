#ifndef ECALGARLICCLUSTERHELPER_HH_
#define ECALGARLICCLUSTERHELPER_HH_

#include <vector>
#include <string>
#include <lcio.h>

#include "ECALGarlicExtendedObjects.hh"
using namespace ECALGarlicExtendedObjects;

class TH2F;
class ECALGarlicAlgorithmParameters;
class ECALGarlicGeometryParameters;

#include "ECALGarlicGeometryHelpers.hh"
#include "ECALGarlicEnergyEstimator.hh"


class IClassifierReader;

using std::vector;

class ECALGarlicClusterHelpers {

public:

  ECALGarlicClusterHelpers(ECALGarlicAlgorithmParameters* pars, ECALGarlicGeometryParameters* geopars) {
    _algoParams=pars;
    _geomParams=geopars;
    _geomHelper = new ECALGarlicGeometryHelpers(_algoParams, _geomParams);
    _energyEstimator = new ECALGarlicEnergyEstimator(_algoParams, _geomParams);
    _nn_is_setup=false;
  }

  ~ECALGarlicClusterHelpers() {}

private:

  IClassifierReader *MLPResponse0_25_B;
  IClassifierReader *MLPResponse0_35_B;
  IClassifierReader *MLPResponse0_5_B;
  IClassifierReader *MLPResponse0_75_B;
  IClassifierReader *MLPResponse1_B;
  IClassifierReader *MLPResponse1_25_B;
  IClassifierReader *MLPResponse1_5_B;
  IClassifierReader *MLPResponse1_75_B;
  IClassifierReader *MLPResponse2_25_B;
  IClassifierReader *MLPResponse3_B;
  IClassifierReader *MLPResponse5_B;
  IClassifierReader *MLPResponse10_B;
  IClassifierReader *MLPResponse20_B;

  IClassifierReader *MLPResponse0_25_EC;
  IClassifierReader *MLPResponse0_35_EC;
  IClassifierReader *MLPResponse0_5_EC;
  IClassifierReader *MLPResponse0_75_EC;
  IClassifierReader *MLPResponse1_EC;
  IClassifierReader *MLPResponse1_25_EC;
  IClassifierReader *MLPResponse1_5_EC;
  IClassifierReader *MLPResponse1_75_EC;
  IClassifierReader *MLPResponse2_25_EC;
  IClassifierReader *MLPResponse3_EC;
  IClassifierReader *MLPResponse5_EC;
  IClassifierReader *MLPResponse10_EC;
  IClassifierReader *MLPResponse20_EC;

  void setupNN();
  bool _nn_is_setup;

  ECALGarlicAlgorithmParameters* _algoParams;
  ECALGarlicGeometryParameters* _geomParams;
  ECALGarlicGeometryHelpers* _geomHelper;
  ECALGarlicEnergyEstimator* _energyEstimator;

  

public:

  void CalculatePhotonProbability(LCEvent *evt, 
				  ExtendedCluster *myCluster, 
				  std::vector<ExtendedCluster*> *clusters, 
				  int clus_ID, 
				  ClusterParameters *clusPar,
				  std::vector<ExtendedTrack*> tracks, 
				  std::map<ExtendedTrack*,std::map<int,std::vector<ExtendedHit*> > > *allRemovedHits);

  void FillClusterParameters(LCEvent *evt,
			     ExtendedCluster *myCluster,
			     std::vector<ExtendedCluster*> *clusters,
			     int clus_ID,
			     ClusterParameters *clusPar,
			     std::vector<ExtendedTrack*> tracks);

  void CalculateSeedDirection(ExtendedCluster* a_cluster);
  void CalculateClusterDirection(ExtendedCluster* a_cluster);

  void FitEllipsoid(ExtendedCluster* a_cluster);

  vec3 FollowClusterDirection(vec3 *myPos, int layersFromStart, 
			      vec3 *clusterDirection, ClusterLocation location, 
			      bool &gapHoleVeto,int &gapJumpTo);
  
  void FreeHits(ExtendedCluster *noCluster);
  void FreeHits(vector<ExtendedHit*> &noCluster);

  bool CheckClusterCriteria(ExtendedCluster *myCluster,vector<ExtendedTrack*> tracks);

  void GetFreeNeighbours(ExtendedHit *a_hit,vector<ExtendedHit* > &myNeighbours,ExtendedCluster *preCluster, int iteration, double Theta);

  void RecAddNeighbourBins(TH2F *histo, int bin_i, vector<int> *seed_c);

  void IterAddNeighbourBins(TH2F *histo, int bin_i, vector<int> *seed_c,int it);

};


#endif
