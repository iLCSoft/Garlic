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

#include "ECALGarlicExtendedObjects.hh"

#include "ECALGarlicGeometryHelpers.hh"


using namespace lcio;
using namespace std;
using namespace ECALGarlicExtendedObjects;

class TTree;

class lcio::CalorimeterHit;
class lcio::ClusterImpl;
class lcio::Track;

class ECALGarlicCluster;

class ExtendedCluster2;
class ExtendedHit2;

class ECALGarlic : public marlin::Processor, public ECALGarlicGeometryHelpers {
  
public:
  
  virtual marlin::Processor * newProcessor() { return new ECALGarlic; }
  
  ECALGarlic();
  ~ECALGarlic();
  
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
  void PreparePreClusters(LCEvent *evt, vector<ExtendedCluster2* > &preClusVec);
  void PrepareTracks(const LCEvent *evt, vector<ExtendedTrack* > &trackVec);
  void PrepareMCTracks(vector<ExtendedTrack* > &trackVec);

  vector < ExtendedTrack* > selectNearbyTracks(ExtendedCluster2* preCluster,  vector <ExtendedTrack* >* trackVec);

  void RemoveElectronHits(vector < pair < ExtendedTrack*, ExtendedCluster2* > > electrons, ExtendedCluster2* preClus );
  void RemoveHitsNearExtrapolatedTracks(LCEvent *evt,
					vector<ExtendedTrack* > &trackVec, 
					vector<ExtendedCluster2* > &preClusVec);

  void cleanup(vector<ExtendedCluster2* > &preClusVec);
  void cleanup(vector<ExtendedTrack* > &trackVec);
  void addCollToEvent(LCCollection* col, string colname, LCEvent* evt);


  //----------------------------------------------------

  // Collection names
  string _mcParticleCollectionName;
  string _simHitCaloHitRelationCollectionName;

  string _ecalPreClusterCollectionName;
  string _LDCTrackCollectionName;
  string _TPCTrackCollectionName;
  string _particleCollectionName;

  string _ecal_ps_col_name;

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
    
  // collections
  LCCollection *_mcParticleColl;

  LCCollectionVec* _track_extrap_col;

  float _thicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _thicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];

  // counters
  int _nEvents;
  int _nPreClusters;
  int _nClusters;
  int _nPhotonClusters;
  int _nOtherClusters;

  ECALGarlicCluster* _clusterer;

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

  float _x_mergeTouchFrac;
  int   _x_initialLayerSeparation;


  int    _x_nIterations;

  bool   _x_rejectMLP;
  vector <float> _x_mlpCuts;

  float _x_trackWindowVeto;


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
  float _x_absorberX0;
  int    _x_symmetry;
  int    _x_nPseudoLayers;
  int    _x_nCellsPerWafer;
  int    _x_nBarEcalLayers;
  int    _x_nEndEcalLayers;

  std::vector< std::vector < float > > _x_barrelStaveDir;

  CellIDDecoder<CalorimeterHit>* _x_defaultDecoder;

  bool _geomSetup;

  TFile* _fhistos;
  string _histFileName;
  int _nSaveHist;
  enum {MAXSAVEHIST=-1};

};

#endif
