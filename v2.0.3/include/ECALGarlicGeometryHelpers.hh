#ifndef ECALGARLICGEOMETRYHELPER_HH_
#define ECALGARLICGEOMETRYHELPER_HH_

#include "ECALGarlicExtendedObjects.hh"
using namespace ECALGarlicExtendedObjects;

#include "ECALGarlicGeometryParameters.hh"
#include "ECALGarlicAlgorithmParameters.hh"

class ECALGarlicGeometryHelpers {

public:

  ECALGarlicGeometryHelpers(ECALGarlicAlgorithmParameters* pars, ECALGarlicGeometryParameters* geopars) {
    _algoParams=pars;
    _geomParams=geopars;
  }

  ~ECALGarlicGeometryHelpers() {}

  double Get2dProjDistance(vec3 *a,vec3 *b);
  double Get3dDistance(vec3 *a,vec3 *b);
  void   GetDistancesBetweenClusters(ExtendedCluster *a,ExtendedCluster *b, double *distances);
  void   AssignPseudoLayer(ExtendedHit* &a_hit);
  void   GetNearestHitInPseudoLayer(int ps_layer,vec3 *seed,ExtendedCluster *preCluster,ExtendedHit *&nearestHit);

private:

  ECALGarlicAlgorithmParameters* _algoParams;
  ECALGarlicGeometryParameters* _geomParams;

};


#endif
