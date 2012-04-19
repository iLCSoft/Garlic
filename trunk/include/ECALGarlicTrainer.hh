#ifndef ECALGARLICTRAINER_HH_
#define ECALGARLICTRAINER_HH_

#include <marlin/Processor.h>
#include <EVENT/MCParticle.h>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <vector>
#include <map>
#include <string>

class ECALGarlicTrainer : public marlin::Processor {

public:

  virtual marlin::Processor * newProcessor() { return new ECALGarlicTrainer; }

  ECALGarlicTrainer();
  ~ECALGarlicTrainer();

  virtual void init();
  virtual void processRunHeader(LCRunHeader * run);
  virtual void processEvent(LCEvent * evt);
  virtual void check(LCEvent * evt);
  virtual void end();

private:

  std::pair <MCParticle*, MCParticle*> findGeneratedAndECALEnterParent(MCParticle* mcp);
  void initTreeVars();

  TFile* _fout;

  enum {nintvar=4, nfloatvar=70};

  TTree* _tree;

  float _energy, _EC1, _EC4, _EC9, _start, _end, _depth, _hitDensity, _enDensity;
  float _costheta, _phi;

  float _halfLength, _fullLength;

  float _enFrac_tube5, _enFrac_tube10, _enFrac_tube15, _enFrac_tube20, _enFrac_tube30;
  float _nFrac_tube5, _nFrac_tube10, _nFrac_tube15, _nFrac_tube20, _nFrac_tube30;

  float _enFrac_x0_5, _enFrac_x0_10, _enFrac_x0_15, _enFrac_x0_20, _enFrac_x0_25;
  float _nFrac_x0_5,  _nFrac_x0_10,  _nFrac_x0_15,  _nFrac_x0_20, _nFrac_x0_25;

  float _enFracRel_x0_5, _enFracRel_x0_10, _enFracRel_x0_15, _enFracRel_x0_20, _enFracRel_x0_25;
  float _nFracRel_x0_5,  _nFracRel_x0_10,  _nFracRel_x0_15,  _nFracRel_x0_20, _nFracRel_x0_25;

  float _Es1,  _Es2,  _Es3,  _distToTrackCOG, _distToTrackPOS, _dirErr,  _Width,  _Eccentricity,  _Volume;
  float _angToTrackPOS;
  int _nhits, _Ns1, _Ns2, _Ns3;
  int _genpdg1, _genpdg2, _ecalpdg1, _ecalpdg2;
  float _genpdg1Frac, _genpdg2Frac, _ecalpdg1Frac, _ecalpdg2Frac;
  float _hitMeanEn, _hitRMSEn, _hitQ1En, _hitQ3En;
  float _fracDim2, _fracDim4, _fracDim8;
  float _transAx1, _transAx2;

  float _longFit_E, _longFit_alpha, _longFit_T, _longFit_off, _longFit_prob;

  float _transRmsMin, _transRmsMax;

  float _mol60min, _mol60max;
  float _mol80min, _mol80max;
  float _mol90min, _mol90max;
  float _mol95min, _mol95max;

  int _genP1, _genP2, _ecalP1, _ecalP2;
  float _genP1Frac, _genP2Frac, _ecalP1Frac, _ecalP2Frac;

  float _nnout;

  std::string _GARLICClusterCollectionName,
    _GARLICClusterParametersCollectionName,
    _GARLICClusterParameterRelationsCollectionName,
    _calohitRelationsCollectionName;

  std::string _outRootFile;

};

#endif
