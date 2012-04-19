#ifndef ECALGARLICEXTENDEDCLUSTER_HH_
#define ECALGARLICEXTENDEDCLUSTER_HH_

#include <vector>
#include <algorithm>
#include <HelixClass.h>

#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

#include "ECALGarlicExtendedTrack.hh"
#include "ECALGarlicExtendedObjects.hh"
#include "ECALGarlicGeometryHelpers.hh"
#include "ECALGarlicNeuralNetwork.hh"
#include "ClusterShapes.h"

#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TGraphErrors.h"

using namespace EVENT;

using namespace ECALGarlicExtendedObjects;
using std::vector;

class ExtendedHit2;

class ExtendedCluster2 : public ECALGarlicGeometryHelpers {

public:

  ExtendedCluster2( ) {init();}

  ~ExtendedCluster2() {
    if (_hEnergy) {delete _hEnergy; _hEnergy=NULL;}
    if (_hHits)   {delete _hHits;   _hHits=NULL;}
    if (_clusterShapes) {delete _clusterShapes; _clusterShapes=NULL;}
    if (_longProf) {delete _longProf; _longProf=NULL;}
  }

  void setHitRelationColl(LCCollection* coll) {_caloHitSimHitRelation=coll;}

  void setGeVtoMIPconversion(float fac) {_enConvFac=fac;}

  void setSeed ( const float* pos ) {*_seededFrom = *pos;}
  void setPreClus( ExtendedCluster2* clus ) {_PreCluster = clus;}

  void setBiggestNeighbour(ExtendedCluster2* bigNei) {_biggestNeighbourCluster=bigNei;}
  void setTracks( vector<ExtendedTrack*> * tracks ) {_tracks = tracks; changed();}
  void setID (int id ) { _clusterID = id;}

  void addHits( std::vector <ExtendedHit2*> hits ) {
    for (size_t i=0; i<hits.size(); i++) 
      _hitVec.push_back(hits[i]);
    changed();
  }

  void addHit( ExtendedHit2* hit ) {
    _hitVec.push_back(hit);    
    changed();
  }

  void removeHit ( ExtendedHit2* hit ) {
    if ( find(_hitVec.begin(), _hitVec.end(), hit)!=_hitVec.end() )
      _hitVec.erase ( find(_hitVec.begin(), _hitVec.end(), hit) );
    changed();
  }

  void freeHits();

  std::vector <ExtendedHit2*> * getHits ()  {return &_hitVec;}

  IMPL::LCGenericObjectImpl getGenericObject();

  // getters
  int getID() {return _clusterID;}

  float getEnergy() {calculateClusterShapes(); return _energy;}
  int   getNhits()  {calculateClusterShapes(); return _nhits;}

  float getTheta()  {calculateClusterShapes(); return _theta;}
  float getPhi()    {calculateClusterShapes(); return _phi;}
  int getZone()     {calculateZone(); return _zone;}
  int getLocation() {return getZone();}

  float getTrackDist_cog()   {calculateTrackDistances(); return _trackDist_cog;}
  float getTrackDist_proj()  {calculateTrackDistances(); return _trackDist_proj;}
  float getTrackDist_min()   {calculateTrackDistances(); return _trackDist_min;}
  float getTrackDist_first() {calculateTrackDistances(); return _trackDist_first;}
  float getTrackAng_proj()   {calculateTrackDistances(); return _trackAng_proj;}

  float getClusterPointAngle() {calculateClusterShapes(); return _clusterPointingAngle;}

  float getEccentricity() {calculateClusterShapes(); return _fitted_eccentricity;}
  float getWidth()        {calculateClusterShapes(); return _fitted_width;}
  float getVolume()       {calculateClusterShapes(); return _fitted_volume;}
  
  float getStart()     {calculateLongShape(); return _start;}
  float getEnd()       {calculateLongShape(); return _end;}
  float getMeanDepth() {calculateLongShape(); return _meandepth;}
  
  float* getCentreOfGravity()   {calculateClusterShapes(); return _centreOfGravity;}
  float* getProjectedPosition() {calculateClusterShapes(); return _projectedPosition;}

  float* getTubeEn() {calculateTubeFracs(); return _tube_eFracs;}
  float* getTubeN()  {calculateTubeFracs(); return _tube_nFracs;}

  float* getLongEn() {calculateLongShape(); return _long_eFracs;}
  float* getLongN()  {calculateLongShape(); return _long_nFracs;}

  float* getRelLongEn() {calculateLongShape(); return _long_rel_eFracs;}
  float* getRelLongN()  {calculateLongShape(); return _long_rel_nFracs;}

  float getHitMeanEn() {calculateHitEnMeasures(); return _hitMeanEn;}
  float getHitRMSEn()  {calculateHitEnMeasures(); return _hitRMSEn;}
  float getHitQ1En()   {calculateHitEnMeasures(); return _hitQ1En;}
  float getHitQ3En()   {calculateHitEnMeasures(); return _hitQ3En;}


  std::pair<float, float> get1dMol60() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol60, _maxMol60); };
  std::pair<float, float> get1dMol80() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol80, _maxMol80); };
  std::pair<float, float> get1dMol90() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol90, _maxMol90); };
  std::pair<float, float> get1dMol95() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol95, _maxMol95); };




  float getNNval() {calculateNNResults(); return _NNval;}
  bool getNNsel()  {calculateNNResults(); return _NNsel;}

  float* getFractalDimension() {calculateFractalDimension(); return _fractalDimension;}

  //  float* getEMFitParameters() {calculateEMFit(); return _EMShowerParameters;}

  float* getLongFitPars() {makeLongitudinalProjectionGraph(); return _longFitPars;}

  std::pair <float, float> getTransverseRMS() {calculateTransverseProjectionShapes(); return std::pair <float, float> (_transRmsMin, _transRmsMax);}

  void saveTransverseProjection(TFile* saveFile, TString hname) {
    calculateTransverseProjectionShapes();
    TH2F* hh = this->GetProjectionHistos().first;
    if (saveFile && hh->GetEntries()>0) {
      saveFile->cd();
      TString hname2 = "transverseProjection_";
      if (hname=="") hname2+=_thisCounter;
      else hname2+=hname;
      hh->SetNameTitle(hname2, hname2);
      // hh->Write();
    }
    return;
  }

  std::pair <float, float> getTransverseAxisLengths() {calculateTransverseProjectionShapes(); 
    return std::pair <float, float> (_transProjAxisLengths[0], _transProjAxisLengths[1]);};

  void saveLongProf(TFile* histsave=NULL, TString name="") {
    makeLongitudinalProjectionGraph();
    if (histsave) {
      histsave->cd();
      if (name=="") name+=counter;
      name = "longProf_"+name;
      _longProf->SetNameTitle(name, name);
      _longProf->Write();
    }
    return;
  }

  ClusterShapes* getClusterShape() {calculateClusterShapes(); return _clusterShapes;};

  enum {
    CLUS2_LOCATION_UNKNOWN  = 0,
    CLUS2_LOCATION_BARREL = 1,
    CLUS2_LOCATION_ENDCAP = 2,
    CLUS2_LOCATION_OVERLAP = 3
  };

  // for sorting by #hits
  static bool moreHits(const ExtendedCluster2 *a, const ExtendedCluster2 *b) { 
    return a->_hitVec.size() > b->_hitVec.size(); 
  }

  static bool higherEnergy(ExtendedCluster2 *a, ExtendedCluster2 *b) { 
    return a->getEnergy() > b->getEnergy(); 
  }
  
  static bool higherRawEnergy(ExtendedCluster2 *a, ExtendedCluster2 *b) { 
    return a->getEnergy() > b->getEnergy(); 
  }

  void GetDistancesToCluster(ExtendedCluster2 *b, double *distances);

  std::pair < TH2F*, TH2I* > GetProjectionHistos(float encut=0, int playercut=1000) {
    MakeProjectionHistos(encut, playercut); 
    return std::pair < TH2F*, TH2I* > (_hEnergy, _hHits);
  }

  void GetGlobalPositionFromLocal(float localx, float localy, float* pos);

  float* getOrigin() {return _origin;}
  float* getLocalXaxis() {return _xprimeaxis;}
  float* getLocalYaxis() {return _yprimeaxis;}

  int* getMCCalPdg() {calculateMCnature(); return _MCCalPdg;}
  float* getMCCalFrac() {calculateMCnature(); return _MCCalFrac;}

  int* getMCGenPdg() {calculateMCnature(); return _MCGenPdg;}
  float* getMCGenFrac() {calculateMCnature(); return _MCGenFrac;}

  IMPL::ClusterImpl getClusterImpl();

private:

  void init();

  void changed() {
    _validClusterShapes=false;
    _validLongProj=false;
    _validEMFit=false;
    _validTransClusterShapes=false;
    _validHitEnMeasures=false;
    _validTrackDistances=false;
    _validTubeFracs=false;
    _validLongShape=false;
    _validZone=false;
    _validNNResults=false;
    //    _validHistos=false;
    _validFractal=false;
    _validMC=false;
    return;
  };

  static int counter;

  void MakeProjectionHistos(float encut, int pseudolayerCut);

  void makeLongitudinalProjectionGraph();


  std::vector<ExtendedHit2* > _hitVec;

  float _seededFrom[3];
  float _dir[3];
  float _seed_dir[3];

  ExtendedCluster2* _PreCluster;
  GhostCluster* _Ghosts;
  float _rawEn;

  bool _validClusterShapes;
  bool _validLongProj;
  bool _validEMFit;
  bool _validTransClusterShapes;
  bool _validHitEnMeasures;
  bool _validTrackDistances;
  bool _validTubeFracs;
  bool _validLongShape;
  bool _validZone;
  bool _validNNResults;
  //  bool _validHistos;
  bool _validFractal;
  bool _validMC;

  vector<ExtendedTrack*> * _tracks;

  int getStack (int layer) {
    int st = layer<10 ? 0 : (layer<25 ? 1 : 2);
    return st;
  }

  ExtendedCluster2* _biggestNeighbourCluster;

  ClusterShapes* _clusterShapes;

  float getDistToClusterLine(ExtendedCluster2* cl);

  void calculateClusterShapes();
  void calculateEMFit();
  void calculateHitEnMeasures();
  void calculateTrackDistances();
  void calculateTubeFracs();
  void calculateLongShape();
  void calculateZone();
  void calculateNNResults();
  void calculateFractalDimension();
  void calculateMCnature();
  void calculateTransverseProjectionShapes();

  float get1dmoliere(TH1F* h, float fraction);


  MCParticleVec generatedDecInCaloMCAncestors( MCParticle* mcp );
  MCParticleVec allMCAncestors( MCParticle* mcp );

  MCParticle* lastGeneratedMCAncestor( MCParticle* mcp );
  MCParticle* firstDecInCaloMCAncestor( MCParticle* mcp );

  enum {MAXHITS=40000};

  vector < float > _radii;
  vector <float> _X0s;


  TH2F* _hEnergy;
  TH2I* _hHits;

  TGraphErrors* _longProf;

  float _origin[3];
  float _xprimeaxis[3];
  float _yprimeaxis[3];

  int _thisCounter;

  float _enConvFac;


  // the cluster paramters
  enum {
    ntubes=5,
    nlong=5,
    nfrac=3,
    nlongfitpars=5
  };

  LCCollection* _caloHitSimHitRelation;

  int _clusterID;

  float _energy;
  int _nhits;

  float _theta, _phi;
  int _zone;

  float _trackDist_cog;
  float _trackDist_proj;
  float _trackDist_min;
  float _trackDist_first;
  float _trackAng_proj;

  float _clusterPointingAngle;

  float _fitted_eccentricity;
  float _fitted_width;
  float _fitted_volume;
  
  float _start, _end, _meandepth;

  float _centreOfGravity[3];
  float _projectedPosition[3];

  float _tube_eFracs[ntubes];
  float _tube_nFracs[ntubes];

  float _long_eFracs[nlong];
  float _long_nFracs[nlong];

  float _long_rel_eFracs[nlong];
  float _long_rel_nFracs[nlong];

  float _hitMeanEn, _hitRMSEn, _hitQ1En, _hitQ3En;

  float _fractalDimension[nfrac];

  float _projectionAxMax;
  float _projectionAxMin;

  float _transRmsMin;
  float _transRmsMax;

  float _NNval;
  bool _NNsel;

  int _MCCalPdg[2];
  float _MCCalFrac[2];

  int _MCGenPdg[2];
  float _MCGenFrac[2];

  float _transProjAxisLengths[2];

  //  float _EMShowerParameters[6];

  float _longFitPars[nlongfitpars];

  float _minMol60, _maxMol60;
  float _minMol80, _maxMol80;
  float _minMol90, _maxMol90;
  float _minMol95, _maxMol95;
  
};
    
#endif
