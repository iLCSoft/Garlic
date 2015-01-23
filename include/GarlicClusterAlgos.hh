#ifndef GARLICCLUSTERALG_HH_
#define GARLICCLUSTERALG_HH_

#include "GarlicExtendedObjects.hh"
using namespace GarlicExtendedObjects;

#include <vector>
#include "GarlicGeometryHelpers.hh"
#include "GarlicExtendedTrack.hh"

#include "TFile.h"

class GarlicExtendedCluster;
class GarlicExtendedHit;

using std::vector;

class GarlicClusterAlgos : public GarlicGeometryHelpers {

public:

  GarlicClusterAlgos() {
    _nSaveHist=0;
    _fhistos=0;
    _saveHists=false;
    _hnn="blah";
    _verbose=false;
  }

  ~GarlicClusterAlgos() {}

  //  std::vector <CalorimeterHit*> getSeeds(GarlicExtendedCluster* preClus);
  std::map <CalorimeterHit*, bool> getSeeds(GarlicExtendedCluster* preClus);
  std::map < CalorimeterHit*, GarlicExtendedCluster* > getClusters( GarlicExtendedCluster* preCluster, std::map < CalorimeterHit*, GarlicExtendedCluster* > cores );
  std::map < CalorimeterHit*, GarlicExtendedCluster* > getCores( GarlicExtendedCluster* preCluster, std::vector <CalorimeterHit*> seeds );

  bool mergeCandidate( GarlicExtendedCluster* primary, GarlicExtendedCluster* secondary );

  //  void mergeSatellites( std::map < CalorimeterHit*, GarlicExtendedCluster* > & clusters );

  //  void mergeSatellitesAndElectrons(std::map < CalorimeterHit*, GarlicExtendedCluster* > & clusters ,
  //				   vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > & allElectrons );

  void mergeAllSatellitesAndElectrons(std::map < CalorimeterHit*, GarlicExtendedCluster* > & clusters ,
  				      vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > & allElectrons );

  void saveHistos(TFile* hfile=NULL, TString hname="blah") {
    _saveHists=(hfile!=0); 
    _fhistos=hfile; 
    _hnn=hname;
  }

  vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > getElectrons( GarlicExtendedCluster* preCluster, vector <GarlicExtendedTrack* >, bool force_a_cluster=false );

  void setVerbose(bool x=true) {_verbose=x;}

private:

  GarlicExtendedCluster* BuildCore(CalorimeterHit* mySeed, GarlicExtendedCluster* preCluster, const float *clusterDir);

  TFile* _fhistos;
  bool _saveHists;
  int _nSaveHist;
  TString _hnn;
  bool _verbose;

};

#endif
