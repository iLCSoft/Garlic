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

#include "GarlicExtendedObjects.hh"

#include "GarlicGeometryHelpers.hh"

using namespace lcio;
using namespace std;
using namespace GarlicExtendedObjects;

class TTree;

class GarlicClusterAlgos;

class GarlicExtendedTrack;
class GarlicExtendedCluster;
class GarlicExtendedHit;
class GarlicClusterEnergyCorrector;

class GarlicProcessor : public marlin::Processor, public GarlicGeometryHelpers {
  
public:
  
  virtual marlin::Processor * newProcessor() { return new GarlicProcessor; }
  
  GarlicProcessor();
  ~GarlicProcessor();
  
  void init();
  void processRunHeader(LCRunHeader * run);
  void processEvent(LCEvent * evt);
  void check(LCEvent * evt);
  void end();
  
private:

  void setup();
  void printMrGarlic();
  void setUpGeometry();

  // main algorithm functions
  void PreparePreClusters(LCEvent *evt, vector<GarlicExtendedCluster* > &preClusVec);
  void PrepareTracks(const LCEvent *evt, vector<GarlicExtendedTrack* > &trackVec);
  void PrepareMCTracks(vector<GarlicExtendedTrack* > &trackVec);

  vector < GarlicExtendedTrack* > selectNearbyTracks(GarlicExtendedCluster* preCluster,  vector <GarlicExtendedTrack* >* trackVec);

  void RemoveElectronHits(vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > electrons, GarlicExtendedCluster* preClus );
  //void RemoveHitsNearExtrapolatedTracks(LCEvent *evt,
  //					vector<GarlicExtendedTrack* > &trackVec, 
  //					vector<GarlicExtendedCluster* > &preClusVec);

  std::map < GarlicExtendedTrack*, GarlicExtendedCluster*> RemoveHitsNearExtrapolatedTracks(vector<GarlicExtendedTrack* > &trackVec, 
										 vector<GarlicExtendedCluster*> &preClusVec);


  void cleanup(vector<GarlicExtendedCluster* > &preClusVec);
  void cleanup(vector<GarlicExtendedTrack* > &trackVec);
  void addCollToEvent(LCCollection* col, string colname, LCEvent* evt);


  //----------------------------------------------------

  // Collection names
  string _mcParticleCollectionName;
  string _simHitCaloHitRelationCollectionName;

  string _ecalPreClusterCollectionName;
  string _TrackCollectionName;
  string _TPCTrackCollectionName;
  string _particleCollectionName;

  string _ecal_ps_col_name;
  string _conversionCollName;

  string _electronCollName;

  string _trkExtrapCollName;
  string _seedCollName;
  string _coreCollName;
  string _clusterCollName;
  string _clparsCollName;
  string _clusterParRelCollName;
  string _seedCoreRelCollName;
  string _seedClusterRelCollName;
  string _removedHitsCollectionName;

  string _electronPFOCollName;
  string _photonPFOCollName;
    
  // collections
  LCCollection *_mcParticleColl;

  LCCollectionVec* _track_extrap_col;
  LCCollectionVec* _cluster_start_col;

  float _thicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _thicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];

  // counters
  int _nEvents;
  int _nPreClusters;
  int _nClusters;
  int _nPhotonClusters;
  int _nOtherClusters;

  GarlicClusterAlgos* _clusterer;

  // parameters
  int   _x_debug;

  bool  _x_removeHitsNearTracks;
  bool  _x_cheatTracks;

  int   _x_nLayersForSeeding;
  int   _x_minHitsForSeeding;
  float _x_seedHitEnergyCut;
  float _x_seedEnergyCut;
  float _x_seedDistanceCut;

  int   _x_nlayersSection1;
  int   _x_maxHoleSection1;
  int   _x_maxHoleSection2;
  float _x_maxCoreDist;
  
  float _x_maxMergeDist;

  //  int    _x_nIterations;
  float _x_TouchingCellDistance;
  float  _x_clusterMaxDist;

  float _x_trackWindowVeto;

  float _x_stochasticTerm;
  float _x_constantTerm;
  float _x_moliereRadius;

  float _x_energy_mip_conversion;

  bool _x_debugCollections;


  float _x_MergeRatioCut             ;
  float _x_MergeEnergyDistFactor     ;
  float _x_MergeDistanceMultiplier   ;
  float _x_MergeAbsoluteLargestDist  ;
  float _x_MergePi0MassLimit         ;
  float _x_MergePi0MaxEnergyImbalance;

  float _x_ElectronTransTubeStepSize;
  int _x_ElectronTransNSteps;


  std::string _x_photonSelFile;

  float _x_forwardTrackAngle;

  CellIDDecoder<CalorimeterHit>* _x_defaultDecoder;

  bool _geomSetup;
  TFile* _fhistos;
  string _histFileName;
  int _nSaveHist;
  enum {MAXSAVEHIST=10};

  GarlicConversionFinder* _convFinder;
  GarlicClusterEnergyCorrector* _energyCorrector;

};

#endif
