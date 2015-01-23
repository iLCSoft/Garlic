#ifndef GARLICEXTENDEDCLUSTER_HH_
#define GARLICEXTENDEDCLUSTER_HH_

#include <vector>
#include <algorithm>
#include <HelixClass.h>

#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

#include "GarlicExtendedTrack.hh"
#include "GarlicExtendedHit.hh"
#include "GarlicExtendedObjects.hh"
#include "GarlicGeometryHelpers.hh"
//#include "ECALGarlicNeuralNetwork.hh"
#include "ClusterShapes.h"

#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TGraphErrors.h"

using namespace EVENT;

using namespace GarlicExtendedObjects;
using std::vector;

class GarlicExtendedCluster : public GarlicGeometryHelpers {

public:

  GarlicExtendedCluster( ) {init();}

  ~GarlicExtendedCluster() {
    cleanUpHistograms();
    if (_clusterShapes) {delete _clusterShapes; _clusterShapes=NULL;}
    if (_longProf) {delete _longProf; _longProf=NULL;}
  }

  enum {TRANS, LONGT, MISCL, POINT, E_ONP, TOTAL};
  enum {
    CLUS2_LOCATION_UNKNOWN  = 0,
    CLUS2_LOCATION_BARREL = 1,
    CLUS2_LOCATION_ENDCAP = 2,
    CLUS2_LOCATION_OVERLAP = 3
  };



  void cleanUpHistograms() {
    if (_hEnergy) {
      delete _hEnergy; _hEnergy=NULL;}
    if (_hHits)   {
      delete _hHits;   _hHits=NULL;}
  }

  void setReferenceDirection( const float* dir ) {
    float tt(0);
    for (int i=0; i<3; i++) {
      _refDirection[i] = dir[i];
      tt+=pow(dir[i], 2);
    }
    if (tt>0.001) _validDirection=true;
  }

  void setHitRelationColl(LCCollection* coll) {_caloHitSimHitRelation=coll;}

  void setSeed ( const float* pos ) {for (int i=0; i<3; i++) _seededFrom[i] = pos[i];}
  void setPreClus( GarlicExtendedCluster* clus ) {_PreCluster = clus;}

  void setBiggestNeighbour(GarlicExtendedCluster* bigNei) {_biggestNeighbourCluster=bigNei;}
  void setTracks( vector<GarlicExtendedTrack*> * tracks ) {_tracks = tracks; changed();}
  void setID (int id ) { _clusterID = id;}

  void setAssociatedTrack( GarlicExtendedTrack* trk ) { _associatedTrack = trk; changed(); }
  GarlicExtendedTrack* getAssociatedTrack () { return _associatedTrack; }

  void addHits( std::vector <GarlicExtendedHit*> hits ) {
    for (size_t i=0; i<hits.size(); i++) {
      // check that it's not already there
      if ( find ( _hitVec.begin(), _hitVec.end(), hits[i] ) != _hitVec.end() ) {
	cout << "WARNING, asking to add already existing hit to cluster! refusing. probably an ERROR ?" << i << " " << hits[i] << " " << hits[i]->getCaloHit()->getEnergy() << endl;
      } else {
	_hitVec.push_back(hits[i]);
      }
    }
    changed();
  }

  void addHit( GarlicExtendedHit* hit ) {
    if ( find ( _hitVec.begin(), _hitVec.end(), hit ) != _hitVec.end() ) {
      cout << "WARNING, asking to add already existing hit to cluster! refusing. probably an ERROR ?" << endl;
    } else {
      _hitVec.push_back(hit);    
      changed();
    }
  }

  void removeHit ( GarlicExtendedHit* hit ) {
    if ( find(_hitVec.begin(), _hitVec.end(), hit)!=_hitVec.end() ) {
      _hitVec.erase ( find(_hitVec.begin(), _hitVec.end(), hit) );
      changed();
    }
  }

  //  float getDistToClusterAxis ( const float* pos, bool assumeIpPointing=false );
  // assumption==0: line from cog to IP
  //             1: measured cluster axis
  //             2: reference direction
  float getDistToClusterAxis ( const float* pos, int directionAssumption );

  void Print();

  void freeHits();

  std::vector <GarlicExtendedHit*> * getHits ()  {return &_hitVec;}

  IMPL::LCGenericObjectImpl getGenericObject();

  // getters
  int getID() {return _clusterID;}

  float getEnergy() {calculateClusterShapes(); return _energy;}
  int   getNhits()  {calculateClusterShapes(); return _nhits;}

  float getTheta()  {calculateClusterShapes(); return _theta;}
  float getPhi()    {calculateClusterShapes(); return _phi;}
  int getZone()     {calculateZone(); return _zone;}
  int getLocation() {return getZone();}


  std::pair < float, float > getCombinedInvariantMass( GarlicExtendedCluster* cl2 );

  float getTrackDist_cog()   {calculateTrackDistances(); return _trackDist_cog;}
  float getTrackDist_proj()  {calculateTrackDistances(); return _trackDist_proj;}
  float getTrackDist_min()   {calculateTrackDistances(); return _trackDist_min;}
  float getTrackDist_first() {calculateTrackDistances(); return _trackDist_first;}
  float getTrackAng_proj()   {calculateTrackDistances(); return _trackAng_proj;}

  GarlicExtendedTrack* getClosestTrack() { calculateTrackDistances(); return _closestTrack; }

  float getClusterPointAngle() {calculateClusterShapes(); return _clusterPointingAngle;}

  float getEccentricity() {calculateClusterShapes(); return _fitted_eccentricity;}
  float getWidth()        {calculateClusterShapes(); return _fitted_width;}
  float getVolume()       {calculateClusterShapes(); return _fitted_volume;}
  
  float getStart()     {calculateLongShape(); return _start;}
  float getEnd()       {calculateLongShape(); return _end;}
  float getMeanDepth() {calculateLongShape(); return _meandepth;}
  float getRelMeanDepth() {calculateLongShape(); return _relmeandepth;}
  
  float* getCentreOfGravity()   {calculateClusterShapes(); return _centreOfGravity;}
  float* getProjectedPosition() {calculateClusterShapes(); return _projectedPosition;}

  float* getTubeEn() {calculateTubeFracs(); return _tube_eFracs;}
  float* getTubeN()  {calculateTubeFracs(); return _tube_nFracs;}

  float* getLongEn() {calculateLongShape(); return _long_eFracs;}
  float* getLongN()  {calculateLongShape(); return _long_nFracs;}

  float* getRelLongEn() {calculateLongShape(); return _long_rel_eFracs;}
  float* getRelLongN()  {calculateLongShape(); return _long_rel_nFracs;}

  float* getRelRelLongEn() {calculateLongShape(); return _long_rel_rel_eFracs;}

  float getHitMeanEn() {calculateHitEnMeasures(); return _hitMeanEn;}
  float getHitRMSEn()  {calculateHitEnMeasures(); return _hitRMSEn;}
  float getHitQ1En()   {calculateHitEnMeasures(); return _hitQ1En;}
  float getHitQ2En()   {calculateHitEnMeasures(); return _hitQ2En;}
  float getHitQ3En()   {calculateHitEnMeasures(); return _hitQ3En;}

  float getClusterMass() {calculateClusterMass(); return _clusterMass;}

  float getFracPseudoLayers() { calculateClusterShapes(); return _fracPLayers; }
  int getNPseudoLayers() { calculateClusterShapes(); return _nPLayers; }
  int getNLayers() { calculateClusterShapes(); return _nLayers; }

  std::pair<float, float> get1dMol60() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol60, _maxMol60); };
  std::pair<float, float> get1dMol80() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol80, _maxMol80); };
  std::pair<float, float> get1dMol90() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol90, _maxMol90); };
  std::pair<float, float> get1dMol95() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_minMol95, _maxMol95); };

  std::pair<float, float> getEarly1dMol90() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_earlyMinMol90, _earlyMaxMol90); };
  std::pair<float, float> getEarly1dMol95() { calculateTransverseProjectionShapes(); return std::pair< float, float > (_earlyMinMol95, _earlyMaxMol95); };

  bool  getIsMipLike()        {calculateTrackLikeVars(); return _isMipLike;}
  int   getInteractionPlayer() {calculateTrackLikeVars(); return _interactionPlayer;}
  float* getInteractionPosition() {calculateTrackLikeVars(); return _interactionPosition;}


  float getNNval() {calculateNNResults(); return _NNval;}
  bool  getNNsel()  {calculateNNResults(); return _NNsel;}

  float getClusterProperty( std::string varname );


  int getPhotonCutSel() {return int( getPhotonCutScore( TOTAL ) ); }

  float getPhotonCutScore(int ivars=TOTAL) {
    calculatePhotonCutSelResults(); 
    switch(ivars) {
    case TRANS:
      return _photonCutSelT;
      break;
    case LONGT:
      return _photonCutSelL;
      break;
    case MISCL:
      return _photonCutSelM;
      break;
    case POINT:
      return _photonCutSelP;
      break;
    case E_ONP:
      cout << "ERROR: ASKING FOR E/P for photon!" << endl;
      return -999;
      break;
    default:
      return _photonCutSel;
    }
  }

  int getElectronCutSel() {return int(getElectronCutScore(TOTAL)); }

  float getElectronCutScore(int ivars=TOTAL) {
    calculateElectronCutSelResults(); 
    switch(ivars) {
    case TRANS:
      return _electronCutSelT;
      break;
    case LONGT:
      return _electronCutSelL;
      break;
    case MISCL:
      return _electronCutSelM;
      break;
    case POINT:
      return _electronCutSelP;
      break;
    case E_ONP:
      return _electronCutSelEP;
      break;
    default:
      return _electronCutSel;
    }
  }

  float* getFractalDimension() {calculateFractalDimension(); return _fractalDimension;}

  //  float* getEMFitParameters() {calculateEMFit(); return _EMShowerParameters;}

  float* getLongFitPars() {makeLongitudinalProjectionGraph(); return _longFitPars;}

  std::pair <float, float> getTransverseRMS() {calculateTransverseProjectionShapes(); return std::pair <float, float> (_transRmsMin, _transRmsMax);}

  void saveTransverseProjection(TFile* saveFile, TString hname) {
    calculateTransverseProjectionShapes();
    //    TH2F* hh = this->GetProjectionHistos().first;
    TH2F* hh = _hEnergy; // DANIEL CHANGED HERE
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

  float getBiggestPlayHole() {calculatePlayerHoles(); return _biggestPlayHole;}
  float getBiggestPlayHole_eFracBefore() {calculatePlayerHoles(); return _biggestPlayHole_eFracBefore;}
  float getBiggestPlayHole_eFracAfter () {calculatePlayerHoles(); return _biggestPlayHole_eFracAfter;}

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


  float* getClusterStartPosition() {calculateClusterStartPosition(); return _clusterStartPos;}
  int getClusterStartPlayer() {calculateClusterStartPosition(); return _clusterStartPlayer;}

  // for sorting by #hits
  static bool moreHits(const GarlicExtendedCluster *a, const GarlicExtendedCluster *b) { 
    return a->_hitVec.size() > b->_hitVec.size(); 
  }

  static bool higherEnergy(GarlicExtendedCluster *a, GarlicExtendedCluster *b) { 
    return a->getEnergy() > b->getEnergy(); 
  }
  
  static bool higherRawEnergy(GarlicExtendedCluster *a, GarlicExtendedCluster *b) { 
    return a->getEnergy() > b->getEnergy(); 
  }

  void GetDistancesToCluster(GarlicExtendedCluster *b, double *distances);

  //  std::pair < TH2F*, TH2I* > GetProjectionHistos(float encut=0, int playercut=1000) {
  std::pair < TH2F*, TH2I* > GetProjectionHistos(float encut=0, float maxX0=1000) {
    MakeProjectionHistos(encut, maxX0); 
    TH2F* he = _hEnergy ? (TH2F*) _hEnergy->Clone(_hEnergy->GetName() + TString("gotten") ) : 0;
    TH2I* hh = _hHits ? (TH2I*) _hHits->Clone(_hHits->GetName() + TString("gotten") ) : 0 ;
    return std::pair < TH2F*, TH2I* > ( he, hh );
  }

  void GetGlobalPositionFromLocal(float localx, float localy, float* pos);

  float* getReferencePoint() {
    return getProjectedPosition();
    //  return _refPoint;
  }
  float* getReferenceDirection() {
    if ( _validDirection ) 
      return _refDirection;
    else 
      return getProjectedPosition(); // assume pointing
  }

  //  float* getOrigin() {return _origin;}
  float* getLocalXaxis() {return _xprimeaxis;}
  float* getLocalYaxis() {return _yprimeaxis;}

  int* getMCCalPdg() {calculateMCnature(); return _MCCalPdg;}
  float* getMCCalFrac() {calculateMCnature(); return _MCCalFrac;}

  int* getMCGenPdg() {calculateMCnature(); return _MCGenPdg;}
  float* getMCGenFrac() {calculateMCnature(); return _MCGenFrac;}

  IMPL::ClusterImpl getClusterImpl();

private:

  void init();

  void changed();

  static int counter;

  //  void MakeProjectionHistos(float encut, int pseudolayerCut);
  void MakeProjectionHistos(float encut, float maxX0);

  void makeLongitudinalProjectionGraph();

  GarlicExtendedTrack* _associatedTrack;

  //  photonPDFs* _photonPDFs;

  std::vector<GarlicExtendedHit* > _hitVec;

  float _seededFrom[3];
  float _dir[3];
  float _seed_dir[3];

  GarlicExtendedCluster* _PreCluster;
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
  bool _validPhotonCutSelResults;
  bool _validElectronCutSelResults;
  bool _validPlayerHoles;
  bool _validFractal;
  bool _validMC;
  bool _validTrackLike;
  bool _validClusterMass;
  bool _validClusterStart;

  int _nPLayers;
  int _nLayers;
  float _fracPLayers;

  bool _isMipLike;
  int _interactionPlayer;
  float _interactionPosition[3];

  vector<GarlicExtendedTrack*> * _tracks;

  GarlicExtendedTrack* _closestTrack;


  int getStack (int layer) {
    int st = layer<10 ? 0 : (layer<25 ? 1 : 2);
    return st;
  }

  GarlicExtendedCluster* _biggestNeighbourCluster;

  ClusterShapes* _clusterShapes;

  float getDistToClusterLine(GarlicExtendedCluster* cl);

  void calculateClusterShapes();
  void calculateEMFit();
  void calculateHitEnMeasures();
  void calculateTrackDistances();
  void calculateTubeFracs();
  void calculateLongShape();
  void calculateZone();
  void calculateNNResults();
  void calculatePhotonCutSelResults();
  void calculateElectronCutSelResults();
  void calculateFractalDimension();
  void calculateMCnature();
  void calculateTransverseProjectionShapes();
  void calculatePlayerHoles(float min_energy_fraction=0.05);
  void calculateTrackLikeVars();
  void calculateClusterMass();
  void calculateClusterStartPosition();

  float get1dmoliere(TH1F* h, float fraction);
  std::vector < std::pair < float, float > > getTransRadii(TH2F* hh);


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

  float _refPoint[3];
  float _refDirection[3];
  bool _validDirection;

  float _clusterStartPos[3];
  int _clusterStartPlayer;

  //  float _origin[3];
  float _xprimeaxis[3];
  float _yprimeaxis[3];

  int _thisCounter;

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
  
  float _start, _end, _relmeandepth, _meandepth;

  float _centreOfGravity[3];
  float _projectedPosition[3];

  float _tube_eFracs[ntubes];
  float _tube_nFracs[ntubes];

  float _long_eFracs[nlong];
  float _long_nFracs[nlong];

  float _long_rel_eFracs[nlong];
  float _long_rel_nFracs[nlong];

  float _long_rel_rel_eFracs[nlong];

  float _hitMeanEn, _hitRMSEn, _hitQ1En, _hitQ2En, _hitQ3En;

  float _fractalDimension[nfrac];

  float _projectionAxMax;
  float _projectionAxMin;

  float _transRmsMin;
  float _transRmsMax;

  float _NNval;
  bool _NNsel;

  float _photonCutSel;   // total
  float _photonCutSelT;  // transverse
  float _photonCutSelL;  // longit
  float _photonCutSelM;  // miscellaneous
  float _photonCutSelP;  // pointing

  float _electronCutSel;
  float _electronCutSelT;
  float _electronCutSelL;
  float _electronCutSelM;
  float _electronCutSelEP;
  float _electronCutSelP;

  int _MCCalPdg[2];
  float _MCCalFrac[2];

  int _MCGenPdg[2];
  float _MCGenFrac[2];

  float _transProjAxisLengths[2];

  //  float _EMShowerParameters[6];

  float _longFitPars[nlongfitpars];

  float _clusterMass;

  float _minMol60, _maxMol60;
  float _minMol80, _maxMol80;
  float _minMol90, _maxMol90;
  float _minMol95, _maxMol95;
  
  float _earlyMinMol95, _earlyMaxMol95;
  float _earlyMinMol90, _earlyMaxMol90;

  float _biggestPlayHole;
  float _biggestPlayHole_eFracBefore;
  float _biggestPlayHole_eFracAfter;



};
    
#endif
