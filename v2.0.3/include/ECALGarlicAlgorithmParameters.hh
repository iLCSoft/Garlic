#ifndef ECALGARLICALGORITHMPARAMS_HH_
#define ECALGARLICALGORITHMPARAMS_HH_

#include <cmath>
#include <vector>
#include <string>

using std::vector;

class ECALGarlicAlgorithmParameters {

public:
  ECALGarlicAlgorithmParameters() {}
  ~ECALGarlicAlgorithmParameters() {}

  void SetEcalPreClusterCollectionName (std::string name) {_ecalPreClusterCollectionName=name;}
  void SetNLayersForSeeding (int n) {_nLayersForSeeding=n;}
  void SetCorrectPhi (bool b) {_correctPhi=b;}
  void SetCorrectTheta (bool b) {_correctTheta=b;}
  void SetDebug (int i) {_debug=i;}
  void SetXY_gap_transit_factor(float f) {_xy_gap_transit_factor=f;}
  void SetNHitsMin(int n) {_nHitsMin=n;}
  void SetMinEnergy(float f) {_minEnergy=f;}
  void SetDistanceFactor(float f) {_distanceFactor=f;}
  void SetMinHits10X0(int i) {_minHits10X0=i;}
  void SetNIterations(int i) {_nIterations=i;}
  void SetMaxSatEn(float f) {_maxSatelliteEn=f;}
  void SetMinSatEn(float f) {_minSatelliteEn=f;}
  void SetApplyGapCorr(bool b) {_applyGapCorrection=b;}
  void SetRemoveHitsNearTracks(bool b) {_removeHitsNearTracks=b;}
  void SetRejectMLP(bool b) {_rejectMLP=b;}
  void SetMerSat(bool b) {_mergeSatellites=b;}
  void SetCheatTracks(int i) {_cheatTracks=i;}
  void SetIncludePreshower(bool b) {_includePreShower=b;}
  void SetMinHitsForSeed(int i) {_minHitsForSeeding=i;}
  void SetCorrectLeakage(bool b) {_correctLeakage=b;}

  void SetToGeVParsBarrel(vector<float> vec) {_toGeVParameters=vec;}
  void SetToGeVParsEndcap(vector<float> vec) {_toGeVParameters_EC=vec;}
  void SetMLPCutsBarrel (vector<float> vec) {_mlpCuts_B=vec;}
  void SetMLPCutsEndcap (vector<float> vec) {_mlpCuts_EC=vec;}
  void SetCorrThetaPars0 (vector<float> vec) {_corrThParameters0=vec;}
  void SetCorrThetaPars1 (vector<float> vec) {_corrThParameters1=vec;}
  void SetCorrPhiPars0 (vector<float> vec) {_corrPhiParameters0=vec;}
  void SetCorrPhiPars1 (vector<float> vec) {_corrPhiParameters1=vec;}
  void SetCorrPhiPars2 (vector<float> vec) {_corrPhiParameters2=vec;}
  void SetAlphaParsBarrel (vector<float> vec) {alp_params_B=vec;}
  void SetAlphaParsEndcap (vector<float> vec) {alp_params_EC=vec;}
  void SetBetaParsBarrel (vector<float> vec) {bet_params_B=vec;}
  void SetBetaParsEndcap (vector<float> vec) {bet_params_EC=vec;}
  void SetGammaParsBarrel (vector<float> vec) {g_params_B=vec;}
  void SetGammaParsEndcap (vector<float> vec) {g_params_EC=vec;}
  void SetDeltaParsBarrel (vector<float> vec) {d_params_B=vec;}
  void SetDeltaParsEndcap (vector<float> vec) {d_params_EC=vec;}
  void SetLambdaParsBarrel (vector<float> vec) {lam_params_B=vec;}
  void SetLambdaParsEndcap (vector<float> vec) {lam_params_EC=vec;}
  void SetF1ParsBarrel (vector<float> vec) {_f1_params_b=vec;}
  void SetF2ParsBarrel (vector<float> vec) {_f2_params_b=vec;}
  void SetF1ParsEndcap (vector<float> vec) {_f1_params_e=vec;}
  void SetF2ParsEndcap (vector<float> vec) {_f2_params_e=vec;}
  void SetPar0FPars (vector<float> vec) {par0_f_params=vec;}
  void SetPar1FPars (vector<float> vec) {par1_f_params=vec;}

  std::string  GetEcalPreClusterCollectionName()  {return _ecalPreClusterCollectionName;}
  int   GetNLayersForSeeding ()                 {return _nLayersForSeeding;}
  bool  GetCorrectPhi ()    {return _correctPhi;}
  bool  GetCorrectTheta ()  {return _correctTheta;}
  int   GetDebug ()         {return _debug;}
  float GetXY_gap_transit_factor() {return _xy_gap_transit_factor;}
  int   GetNHitsMin()       {return _nHitsMin;}
  float GetMinEnergy()      {return _minEnergy;}
  float GetDistanceFactor() {return _distanceFactor;}
  int   GetMinHits10X0()   {return _minHits10X0;}
  int   GetNIterations()    {return _nIterations;}
  float GetMaxSatEn()       {return _maxSatelliteEn;}
  float GetMinSatEn()       {return _minSatelliteEn;}
  bool  GetApplyGapCorr()   {return _applyGapCorrection;}
  bool  GetRemoveHitsNearTracks() {return _removeHitsNearTracks;}
  bool  GetRejectMLP()      {return _rejectMLP;}
  bool  GetMerSat()         {return _mergeSatellites;}
  int   GetCheatTracks()    {return _cheatTracks;}
  bool  GetIncludePreshower() {return _includePreShower;}
  int   GetMinHitsForSeed() {return _minHitsForSeeding;}
  bool  GetCorrectLeakage() {return _correctLeakage;}
  
  vector <float> GetToGeVParsBarrel  () {return _toGeVParameters;}
  vector <float> GetToGeVParsEndcap  () {return _toGeVParameters_EC;}
  vector <float> GetMLPCutsBarrel    () {return _mlpCuts_B;}
  vector <float> GetMLPCutsEndcap    () {return _mlpCuts_EC;}
  vector <float> GetCorrThetaPars0   () {return _corrThParameters0;}
  vector <float> GetCorrThetaPars1   () {return _corrThParameters1;}
  vector <float> GetCorrPhiPars0     () {return _corrPhiParameters0;}
  vector <float> GetCorrPhiPars1     () {return _corrPhiParameters1;}
  vector <float> GetCorrPhiPars2     () {return _corrPhiParameters2;}
  vector <float> GetAlphaParsBarrel  () {return alp_params_B;}
  vector <float> GetAlphaParsEndcap  () {return alp_params_EC;}
  vector <float> GetBetaParsBarrel   () {return bet_params_B;}
  vector <float> GetBetaParsEndcap   () {return bet_params_EC;}
  vector <float> GetGammaParsBarrel  () {return g_params_B;}
  vector <float> GetGammaParsEndcap  () {return g_params_EC;}
  vector <float> GetDeltaParsBarrel  () {return d_params_B;}
  vector <float> GetDeltaParsEndcap  () {return d_params_EC;}
  vector <float> GetLambdaParsBarrel () {return lam_params_B;}
  vector <float> GetLambdaParsEndcap () {return lam_params_EC;}
  vector <float> GetF1ParsBarrel     () {return _f1_params_b;}
  vector <float> GetF2ParsBarrel     () {return _f2_params_b;}
  vector <float> GetF1ParsEndcap     () {return _f1_params_e;}
  vector <float> GetF2ParsEndcap     () {return _f2_params_e;}
  vector <float> GetPar0FPars        () {return par0_f_params;}
  vector <float> GetPar1FPars        () {return par1_f_params;}

private:
  std::string _ecalPreClusterCollectionName;

  int _nLayersForSeeding;
  bool _correctPhi;
  bool _correctTheta;
  int _debug;
  float _xy_gap_transit_factor;
  int _nHitsMin;
  float _minEnergy;
  double _distanceFactor;
  int _minHits10X0;
  int _nIterations;
  float _maxSatelliteEn;
  float _minSatelliteEn;
  bool _applyGapCorrection;
  bool _removeHitsNearTracks;
  bool _rejectMLP;
  bool _mergeSatellites;
  int _cheatTracks;
  bool _includePreShower;
  int _minHitsForSeeding;
  bool _correctLeakage;

  vector<float> _toGeVParameters;
  vector<float> _toGeVParameters_EC;
  vector<float> _mlpCuts_B;
  vector<float> _mlpCuts_EC;

  vector<float> _corrThParameters0;
  vector<float> _corrThParameters1;

  vector<float> _corrPhiParameters0;
  vector<float> _corrPhiParameters1;
  vector<float> _corrPhiParameters2;
  
  vector<float> alp_params_B;
  vector<float> alp_params_EC;

  vector<float> bet_params_B;
  vector<float> bet_params_EC;

  vector<float> g_params_B;
  vector<float> g_params_EC;

  vector<float> d_params_B;
  vector<float> d_params_EC;

  vector<float> lam_params_B;
  vector<float> lam_params_EC;

  vector <float> _f1_params_b;
  vector <float> _f2_params_b;
  vector <float> _f1_params_e;
  vector <float> _f2_params_e;
  
  vector <float> par0_f_params;
  vector <float> par1_f_params;


};

#endif
