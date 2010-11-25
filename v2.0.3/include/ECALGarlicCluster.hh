#ifndef ECALGARLICCLUSTER_HH_
#define ECALGARLICCLUSTER_HH_

#include "ECALGarlicExtendedObjects.hh"

class ECALGarlicAlgorithmParameters;
class ECALGarlicGeometryParameters;

using namespace ECALGarlicExtendedObjects;

#include <vector>
#include "ECALGarlicClusterHelpers.hh"
#include "ECALGarlicGeometryHelpers.hh"

using std::vector;

class ECALGarlicCluster {


public:

  ECALGarlicCluster(ECALGarlicAlgorithmParameters* pars, ECALGarlicGeometryParameters* geompars) {
    _algoParams=pars;
    _geomParams=geompars;
    _clusterHelper = new ECALGarlicClusterHelpers(_algoParams, _geomParams);
    _geomHelper = new ECALGarlicGeometryHelpers(_algoParams, _geomParams);
  }

  ~ECALGarlicCluster() {}

  void Cluster(ExtendedCluster &preCluster, vector<vec3> &possibleSeeds, 
	       vector<ExtendedCluster* > &clusters,vector<ExtendedTrack*> tracks );

  void AddCore(vector<ExtendedHit* > &aCluster, vec3 *mySeed, ExtendedCluster &preCluster, vec3 *clusterDir);

  void BuildClusterFromNeighbours(vector<ExtendedHit* > &myCluster, vec3 *clusterDir, ExtendedCluster &preCluster);

  void MergeSatellites(vector<ExtendedCluster* > *clusters);

private:
  ECALGarlicAlgorithmParameters* _algoParams;
  ECALGarlicGeometryParameters* _geomParams;
  ECALGarlicClusterHelpers* _clusterHelper;
  ECALGarlicGeometryHelpers* _geomHelper;
};

#endif
