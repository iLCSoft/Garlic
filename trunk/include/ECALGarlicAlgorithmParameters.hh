#ifndef ECALGARLICALGORITHMPARAMS_HH_
#define ECALGARLICALGORITHMPARAMS_HH_

#include <cmath>
#include <vector>
#include <string>

using std::vector;

class ECALGarlicAlgorithmParameters {

public:

  static ECALGarlicAlgorithmParameters& Instance() {
    static ECALGarlicAlgorithmParameters algopars;
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

  void SetClusterNIterations(int i) {_nIterations=i;}

  void SetMergeTouchFraction(float f) {_mergeTouchFraction=f;}
  void SetMergeInitalLayerSeparation(int   i) {_mergeInitialLayerSeparation=i;}

  void SetRejectMLP(bool b) {_rejectMLP=b;}
  void SetMLPCuts (vector<float> vec) {_mlpCuts=vec;}


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

  int   GetClusterNIterations()    {return _nIterations;}

  float GetMergeTouchFraction() {return _mergeTouchFraction;}
  int   GetMergeInitalLayerSeparation() {return _mergeInitialLayerSeparation;}

  bool  GetRejectMLP()      {return _rejectMLP;}
  vector <float>* GetMLPCuts() {return &_mlpCuts;}

private:

  ECALGarlicAlgorithmParameters() {}
  ECALGarlicAlgorithmParameters(const ECALGarlicAlgorithmParameters&);
  ECALGarlicAlgorithmParameters& operator= (const ECALGarlicAlgorithmParameters&);

  std::string _ecalPreClusterCollectionName;

  int _debug;

  int _cheatTracks;
  bool _removeHitsNearTracks;
  float _trackWindowVeto;

  int _nLayersForSeeding;
  int _minHitsForSeeding;
  float _seedHitEnergyCut;
  float _seedEnergyCut;
  float _seedDistanceCut;

  int _coreLayerSec1;
  int _coreMaxHoleSec1;
  int _coreMaxHoleSec2;

  int _nIterations;

  float _mergeTouchFraction;
  int _mergeInitialLayerSeparation;

  bool _rejectMLP;
  vector<float> _mlpCuts;

};

#endif
