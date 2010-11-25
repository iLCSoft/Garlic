#ifndef ECALGARLIC_HH_
#define ECALGARLIC_HH_

#include <cmath>
#include <string>
#include <vector>

#include <lcio.h>
#include <EVENT/Track.h>

#include <marlin/Processor.h>

#include "TClonesArray.h"
#include "Rtypes.h"
#include <TMVA/Reader.h>

#include <ClusterParameters.hh>
#include <MCPhoton.hh>

#include "ECALGarlicExtendedObjects.hh"

#include "ECALGarlicGeometryHelpers.hh"


using namespace lcio;
using namespace std;
using namespace ECALGarlicExtendedObjects;

class TTree;

class lcio::CalorimeterHit;
class lcio::ClusterImpl;
class lcio::Track;

class ECALGarlicClusterHelpers;
class ECALGarlicCluster;
class ECALGarlicAlgorithmParameters;
class ECALGarlicGeometryParameters;
class ECALGarlicEnergyEstimator;

class ECALGarlic : public marlin::Processor {
  
  
public:
  
  virtual marlin::Processor * newProcessor() { return new ECALGarlic; }
  
  ECALGarlic();
  
  virtual void init();
  virtual void processRunHeader(LCRunHeader * run);
  virtual void processEvent(LCEvent * evt);
  virtual void check(LCEvent * evt);
  virtual void end();
  
private:

  void printMrGarlic();
  void setUpGeometry();

  void AddDeeperHits(ExtendedHit *myHit, vec3 *clusterDir, vector<ExtendedHit* > &myCluster, ExtendedCluster &preCluster);

  bool RejectByMLPCut(ExtendedCluster *myCluster);

  // main algorithm functions
  void PreparePreClusters(LCEvent *evt, vector<ExtendedCluster* > &preClusVec);
  void PrepareTracks(const LCEvent *evt, vector<ExtendedTrack* > &trackVec);

  void RemoveHitsNearExtrapolatedTracks(LCEvent *evt,
					vector<ExtendedTrack* > &trackVec, 
					vector<ExtendedCluster* > &preClusVec, 
					map<ExtendedTrack*,
					map<int,vector<ExtendedHit*> > > *allRemovedHits);

  void ProjectForSeeding(ExtendedCluster &preClus, vector<vec3> &possibleSeeds, LCCollectionVec *seed_col);

  void WriteClusters(LCEvent *evt, map<int, vector<ExtendedCluster* > *> *clus_map);

  void cleanup(map<int,vector<ExtendedCluster* > *> *clusMap, vector<ExtendedCluster* > &preClusVec, vector<ExtendedTrack* > &trackVec);

  //----------------------------------------------------


  // Collection names
  string _mcParticleCollectionName;
  string _ecalBarrelPreShowerHitCollectionName;
  string _ecalEndcapPreShowerHitCollectionName;
  string _ecalPhotonClusterCollectionName;
  string _LDCTrackCollectionName;
  string _removedHitsCollectionName;
  string _particleCollectionName;

  // collections
  LCCollection *_mcParticleColl;

  float _thicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _thicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];

  // counters
  int _nEvents;
  int _nPreClusters;
  int _nClusters;
  int _nPhotonClusters;
  int _nOtherClusters;

  ECALGarlicAlgorithmParameters* _garPars;
  ECALGarlicGeometryParameters* _geomPars;

  ECALGarlicGeometryHelpers* _geomHelpers;

  ECALGarlicClusterHelpers* _clusterHelpers;

  ECALGarlicCluster* _clusterer;

  ECALGarlicEnergyEstimator* _energyEstimator;


  // parameters
  std::string _x_ecalPreClusterCollectionName;

  int    _x_nLayersForSeeding;
  bool   _x_correctPhi;
  bool   _x_correctTheta;
  int    _x_debug;
  float  _x_xy_gap_transit_factor;
  int    _x_nHitsMin;
  float  _x_minEnergy;
  double _x_distanceFactor;
  int    _x_minHits10X0;
  int    _x_nIterations;
  float  _x_maxSatelliteEn;
  float  _x_minSatelliteEn;
  bool   _x_applyGapCorrection;
  bool   _x_removeHitsNearTracks;
  bool   _x_rejectMLP;
  bool   _x_mergeSatellites;
  int    _x_cheatTracks;
  bool   _x_includePreShower;
  int    _x_minHitsForSeeding;
  bool   _x_correctLeakage;

  vector <float> _x_toGeVParameters;
  vector <float> _x_toGeVParameters_EC;
  vector <float> _x_mlpCuts_B;
  vector <float> _x_mlpCuts_EC;
  vector <float> _x_corrThParameters0;
  vector <float> _x_corrThParameters1;
  vector <float> _x_corrPhiParameters0;
  vector <float> _x_corrPhiParameters1;
  vector <float> _x_corrPhiParameters2;
  vector <float> _x_alp_params_B;
  vector <float> _x_alp_params_EC;
  vector <float> _x_bet_params_B;
  vector <float> _x_bet_params_EC;
  vector <float> _x_g_params_B;
  vector <float> _x_g_params_EC;
  vector <float> _x_d_params_B;
  vector <float> _x_d_params_EC;
  vector <float> _x_lam_params_B;
  vector <float> _x_lam_params_EC;
  vector <float> _x_f1_params_b;
  vector <float> _x_f2_params_b;
  vector <float> _x_f1_params_e;
  vector <float> _x_f2_params_e;
  vector <float> _x_par0_f_params;
  vector <float> _x_par1_f_params;


  double _x_rOfBarrel;
  double _x_zOfBarrel;
  double _x_rOfEndcap;
  double _x_zOfEndcap;
  float  _x_activeThickness;
  float  _x_passiveThickness;
  float  _x_firstEndcapLayerOffset;
  float  _x_firstBarrelLayerOffset;
  float  _x_bField;
  float  _x_positionBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float  _x_absThicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float  _x_absThicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];
  float  _x_padSizeEcal[MAX_NUMBER_OF_LAYERS];
  float  _x_positionEndcapLayer[MAX_NUMBER_OF_LAYERS];
  float  _x_rInnerEcalEndcap;
  float  _x_cosOfBarrel;
  float  _x_guardringSize;
  float  _x_fiberSize;
  float  _x_fiberSizeModule;
  int    _x_symmetry;
  int    _x_nPseudoLayers;
  int    _x_nCellsPerWafer;

  std::vector<vec3> _x_barrelStaveDir;

  CellIDDecoder<CalorimeterHit>* _x_defaultDecoder;



};

#endif
