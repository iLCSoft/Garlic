#ifndef ECALGARLICENERGYESTIMATOR_HH_
#define ECALGARLICENERGYESTIMATOR_HH_

#include "TF1.h"
#include <vector>
#include <cmath>
#include "ClusterParameters.hh"
#include <lcio.h>

#include "ECALGarlicExtendedObjects.hh"
using namespace ECALGarlicExtendedObjects;

class ECALGarlicGeometryParameters;
class ECALGarlicAlgorithmParameters;
class ECALGarlicGeometryHelpers;

class ECALGarlicEnergyEstimator {

public:

  ECALGarlicEnergyEstimator(ECALGarlicAlgorithmParameters* pars, ECALGarlicGeometryParameters* geopars);
  ~ECALGarlicEnergyEstimator() {}

private:

  ECALGarlicAlgorithmParameters* _algoParams;
  ECALGarlicGeometryParameters* _geomParams;

  ECALGarlicGeometryHelpers* _geomHelper;

  double toGeVfunction(double x, std::vector <float> * parameters);

  TF1* _corrTh;
  TF1* _corrThPar0;
  TF1* _corrThPar1;
  TF1* _corrPhi;
  TF1* _corrPhiPar0;
  TF1* _corrPhiPar1;
  TF1* _corrPhiPar2;

  TF1* _alp;
  TF1* _bet;
  TF1* _g;
  TF1* _d;
  TF1* _lam;

  TF1 *_f1_b;
  TF1 *_f2_b;
  TF1 *_f1_e;
  TF1 *_f2_e;

  TF1 *par1_f;
  TF1 *par0_f;
  TF1 *l_corr;

public:
  double getGeV(double en, bool isBarrel);
  double correctEnergyTheta();
  double correctEnergyPhi();
  void CorrectLeakage(ClusterParameters *clusPar);
  double correctEnergyTheta(double in_en, double theta);
  double correctEnergyPhi(double in_en, double phi);
  double EstimateEnergyByMix(ClusterParameters *clusPar, double in_en);
  void ApplyFullGapCorrection(LCEvent *evt,ExtendedCluster *myCluster, GhostCluster *myGhostCluster);

};

#endif
