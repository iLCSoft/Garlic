#ifndef ECALGARLICALGORITHMPARAMS_HH_
#define ECALGARLICALGORITHMPARAMS_HH_

#include <cmath>
#include <vector>
#include <string>

#include "GarlicClusterSelector.hh"


using std::vector;

class GarlicAlgorithmParameters {

public:

  static GarlicAlgorithmParameters& Instance() {
    static GarlicAlgorithmParameters algopars;
    return algopars;
  }

  void SetEcalPreClusterCollectionName (std::string name) {_ecalPreClusterCollectionName=name;}

  void SetDebug (int i) {_debug=i;}

  void SetTrackCheat(int i) {_cheatTracks=i;}
  void SetTrackVetoWindow(float f) {_trackWindowVeto=f;}
  void SetTrackRemoveNearbyHits(bool b) {_removeHitsNearTracks=b;}

  void SetSeedNLayers(int n) {_nLayersForSeeding=n;}
  void SetSeedMinHits(int i) {_minHitsForSeeding=i;}
  void SetSeedHitEnergyCut(float f) {_seedHitEnergyCut=f;}
  void SetSeedEnergyCut   (float f) {_seedEnergyCut=f;}
  void SetSeedDistanceCut (float f) {_seedDistanceCut=f;}
    
  void SetCoreLayersSection1 (int i) {_coreLayerSec1 = i;}
  void SetCoreMaxHoleSection1(int i) {_coreMaxHoleSec1 = i;}
  void SetCoreMaxHoleSection2(int i) {_coreMaxHoleSec2 = i;}
  void SetCoreDistanceCut(float f) {_coreDistanceCut=f;}

  void  SetTouchingCellDistance(float f) {_touchingCellDist=f;}
  float GetTouchingCellDistance() {return _touchingCellDist;}

  void SetClusterMaxDist(float f) {_clusterMaxDist=f;}

  void SetForwardTrackAngle(float f) {_forwardTrackAngle=f;}

  void SetStochasticTerm ( float f ) {_stochasticTerm=f;}
  void SetConstantTerm   ( float f ) {_constantTerm  =f;}
  void SetMoliereRadius  ( float f ) {_moliereRadius =f;}

  void Set_MaxMergeDist (float x) {_max_merge_dist=x;}
  float Get_MaxMergeDist () {return _max_merge_dist;}
  
  // for converting hit energies to MIPs - for first layer: scale aothers according to abs thickness
  void SetEnergyMIPconversion(float x) {_en_mip_conv=x;}
  float GetEnergyMIPconversion() {return _en_mip_conv;}

  void SetMergeRatioCut              (float x ) { _merge_ratioCut                 = x; };
  void SetMergeEnergyDistFactor	     (float x ) { _merge_energy_dist_factor       = x; };
  void SetMergeAbsoluteLargestDist   (float x ) { _merge_absoluteLargestDist      = x; };
  void SetMergePi0MassLimit	     (float x ) { _merge_pi0_mass_limit           = x; };
  void SetMergePi0MaxEnergyImbalance (float x ) { _merge_pi0_max_energy_imbalance = x; };
  
  float GetMergeRatioCut              () { return _merge_ratioCut                ; };
  float GetMergeEnergyDistFactor      () { return _merge_energy_dist_factor      ; };
  float GetMergeAbsoluteLargestDist   () { return _merge_absoluteLargestDist     ; };
  float GetMergePi0MassLimit	      () { return _merge_pi0_mass_limit          ; };
  float GetMergePi0MaxEnergyImbalance () { return _merge_pi0_max_energy_imbalance; };
  
  void SetElectronTransTubeStepSize(float x) { _ElectronTransTubeStepSize=x; }
  float GetElectronTransTubeStepSize() {return _ElectronTransTubeStepSize; }

  void SetElectronTransNSteps (int x) { _ElectronTransNSteps =x; }
  int GetElectronTransNSteps () { return _ElectronTransNSteps; }


  std::string  GetEcalPreClusterCollectionName()  {return _ecalPreClusterCollectionName;}
  int   GetDebug ()         {return _debug;}

  bool  GetTrackRemoveNearbyHits() {return _removeHitsNearTracks;}
  int   GetTrackCheat()    {return _cheatTracks;}
  float GetTrackVetoWindow() {return _trackWindowVeto;}

  int   GetSeedNLayers()                 {return _nLayersForSeeding;}
  int   GetSeedMinHits() {return _minHitsForSeeding;}
  float GetSeedHitEnergyCut() {return _seedHitEnergyCut;}
  float GetSeedEnergyCut   () {return _seedEnergyCut;}
  float GetSeedDistanceCut () {return _seedDistanceCut;}

  int   GetCoreLayersSection1()  {return _coreLayerSec1;}
  int   GetCoreMaxHoleSection1() {return _coreMaxHoleSec1;}
  int   GetCoreMaxHoleSection2() {return _coreMaxHoleSec2;}
  float GetCoreDistanceCut()     {return _coreDistanceCut;}

  float GetClusterMaxDist() {return _clusterMaxDist;}

  float GetStochasticTerm () {return _stochasticTerm;}
  float GetConstantTerm   () {return _constantTerm  ;}
  float GetMoliereRadius  () {return _moliereRadius ;}

  void SetPhotonCutFile(std::string x) {_cutFileName=x;}
  std::string GetPhotonCutFile() {return _cutFileName;}

  GarlicClusterSelector* GetClusterSelector() {
    if ( !_clustersel ) {
      bool verbosesel=GetDebug()>0;
      _clustersel = new GarlicClusterSelector( _cutFileName, -999, verbosesel);
    }
    return _clustersel;
  }

  float GetForwardTrackAngle() {return _forwardTrackAngle;}


private:

  GarlicAlgorithmParameters() {}
  GarlicAlgorithmParameters(const GarlicAlgorithmParameters&);
  GarlicAlgorithmParameters& operator= (const GarlicAlgorithmParameters&);

  std::string _ecalPreClusterCollectionName;
  std::string _cutFileName;

  GarlicClusterSelector* _clustersel;

  int _debug;

  int _cheatTracks;
  bool _removeHitsNearTracks;
  float _trackWindowVeto;

  int _nLayersForSeeding;
  int _minHitsForSeeding;
  float _seedHitEnergyCut;
  float _seedEnergyCut;
  float _seedDistanceCut;

  float _coreDistanceCut;
  int _coreLayerSec1;
  int _coreMaxHoleSec1;
  int _coreMaxHoleSec2;

  float _clusterMaxDist;

  float _mergeTouchFraction;
  int _mergeInitialLayerSeparation;

  float _max_merge_dist;

  float _stochasticTerm;
  float _constantTerm;
  float _moliereRadius;

  float _en_mip_conv;

  float _forwardTrackAngle;

  float _touchingCellDist;

  float _merge_ratioCut;
  float _merge_energy_dist_factor;
  float _merge_absoluteLargestDist;
  float _merge_pi0_mass_limit;
  float _merge_pi0_max_energy_imbalance;

  float _ElectronTransTubeStepSize;
  int   _ElectronTransNSteps;

};

#endif
