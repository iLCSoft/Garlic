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
  void SetCoreDistanceCut(float f) {_coreDistanceCut=f;}

  //  void SetClusterNIterations(int i) {_nIterations=i;}
  void SetClusterMaxDist(float f) {_clusterMaxDist=f;}


  void SetMergeTouchFraction(float f) {_mergeTouchFraction=f;}
  void SetMergeInitalLayerSeparation(int   i) {_mergeInitialLayerSeparation=i;}

  void SetStochasticTerm ( float f ) {_stochasticTerm=f;}
  void SetConstantTerm   ( float f ) {_constantTerm  =f;}
  void SetMoliereRadius  ( float f ) {_moliereRadius =f;}

  void SetRequirePointing( bool x ) {_requirePointing=x;}

  void SetRejectMLP(bool b) {_rejectMLP=b;}
  void SetMLPCuts (vector<float> vec) {_mlpCuts=vec;}


  void  Set_photonPDF_loosecut_lowE ( float x) {_photonPDF_loosecut_lowE=x;}
  void  Set_photonPDF_loosecut_hiE  ( float x) {_photonPDF_loosecut_highE=x;}
  void  Set_photonPDF_tightcut_lowE ( float x) {_photonPDF_tightcut_lowE=x;}
  void  Set_photonPDF_tightcut_hiE  ( float x) {_photonPDF_tightcut_highE=x;}
  float Get_photonPDF_loosecut_lowE ( ) {return _photonPDF_loosecut_lowE;}
  float Get_photonPDF_loosecut_hiE  ( ) {return _photonPDF_loosecut_highE;}
  float Get_photonPDF_tightcut_lowE ( ) {return _photonPDF_tightcut_lowE;}
  float Get_photonPDF_tightcut_hiE  ( ) {return _photonPDF_tightcut_highE;}

  void  Set_electronPDF_loosecut_lowE ( float x) {_electronPDF_loosecut_lowE=x;}
  void  Set_electronPDF_loosecut_hiE  ( float x) {_electronPDF_loosecut_highE=x;}
  void  Set_electronPDF_tightcut_lowE ( float x) {_electronPDF_tightcut_lowE=x;}
  void  Set_electronPDF_tightcut_hiE  ( float x) {_electronPDF_tightcut_highE=x;}
  float Get_electronPDF_loosecut_lowE ( ) {return _electronPDF_loosecut_lowE;}
  float Get_electronPDF_loosecut_hiE  ( ) {return _electronPDF_loosecut_highE;}
  float Get_electronPDF_tightcut_lowE ( ) {return _electronPDF_tightcut_lowE;}
  float Get_electronPDF_tightcut_hiE  ( ) {return _electronPDF_tightcut_highE;}

  //  void Set_MaxTransMergeDist (float x) {_max_trans_dist=x;}
  //  float Get_MaxTransMergeDist () {return _max_trans_dist;}
  
  void Set_MaxMergeDist (float x) {_max_merge_dist=x;}
  float Get_MaxMergeDist () {return _max_merge_dist;}
  

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

  //  int   GetClusterNIterations()    {return _nIterations;}
  float GetClusterMaxDist() {return _clusterMaxDist;}

  float GetMergeTouchFraction() {return _mergeTouchFraction;}
  int   GetMergeInitalLayerSeparation() {return _mergeInitialLayerSeparation;}

  bool  GetRejectMLP()      {return _rejectMLP;}
  vector <float>* GetMLPCuts() {return &_mlpCuts;}

  bool  GetRequirePointing() {return _requirePointing;}
  float GetStochasticTerm () {return _stochasticTerm;}
  float GetConstantTerm   () {return _constantTerm  ;}
  float GetMoliereRadius  () {return _moliereRadius ;}


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

  float _coreDistanceCut;
  int _coreLayerSec1;
  int _coreMaxHoleSec1;
  int _coreMaxHoleSec2;

  //  int _nIterations;
  float _clusterMaxDist;

  float _mergeTouchFraction;
  int _mergeInitialLayerSeparation;

  bool _rejectMLP;
  vector<float> _mlpCuts;

  float _photonPDF_loosecut_lowE;
  float _photonPDF_loosecut_highE;
  float _photonPDF_tightcut_lowE;
  float _photonPDF_tightcut_highE;

  float _electronPDF_loosecut_lowE;
  float _electronPDF_loosecut_highE;
  float _electronPDF_tightcut_lowE;
  float _electronPDF_tightcut_highE;

  float _max_merge_dist;

  float _stochasticTerm;
  float _constantTerm;
  float _moliereRadius;

  bool _requirePointing;

};

#endif
