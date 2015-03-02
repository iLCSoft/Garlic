#ifndef ECALGARLICCLUSTER_HH_
#define ECALGARLICCLUSTER_HH_

#include "ECALGarlicExtendedObjects.hh"
using namespace ECALGarlicExtendedObjects;

#include <vector>
//#include "ECALGarlicClusterHelpers.hh"
#include "ECALGarlicGeometryHelpers.hh"
#include "ECALGarlicExtendedTrack.hh"

#include "TFile.h"

class ExtendedCluster2;
class ExtendedHit2;

using std::vector;

class ECALGarlicCluster : public ECALGarlicGeometryHelpers {

public:

  ECALGarlicCluster() {
    _nSaveHist=0;
    _fhistos=0;
    _saveHists=false;
    _hnn="blah";
    _verbose=false;
  }

  ~ECALGarlicCluster() {}

  //  std::vector <CalorimeterHit*> getSeeds(ExtendedCluster2* preClus);
  std::map <CalorimeterHit*, bool> getSeeds(ExtendedCluster2* preClus);
  std::map < CalorimeterHit*, ExtendedCluster2* > getClusters( ExtendedCluster2* preCluster, std::map < CalorimeterHit*, ExtendedCluster2* > cores );
  std::map < CalorimeterHit*, ExtendedCluster2* > getCores( ExtendedCluster2* preCluster, std::vector <CalorimeterHit*> seeds );

  bool mergeCandidate( ExtendedCluster2* primary, ExtendedCluster2* secondary );

  void mergeSatellites( std::map < CalorimeterHit*, ExtendedCluster2* > & clusters );
  void mergeSatellites_byPDF( std::map < CalorimeterHit*, ExtendedCluster2* > & clusters, bool usepdforcuts );
  void mergeSatellitesAndElectrons_byPDF(std::map < CalorimeterHit*, ExtendedCluster2* > & clusters ,
					 vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > & allElectrons );

  void mergeUnclusteredAndElectronsPhotons_byPDF(std::vector < ExtendedHit2* > & unclusteredHits, 
						 vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > & allElectrons,
						 std::map < CalorimeterHit*, ExtendedCluster2* > & clusters );


  void saveHistos(TFile* hfile=NULL, TString hname="blah") {
    _saveHists=(hfile!=0); 
    _fhistos=hfile; 
    _hnn=hname;
  }

  vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > getElectrons( ExtendedCluster2* preCluster, vector <ExtendedTrack* > );

  void setVerbose(bool x=true) {_verbose=x;}

private:

  ExtendedCluster2* BuildCore(CalorimeterHit* mySeed, ExtendedCluster2* preCluster, const float *clusterDir);

  TFile* _fhistos;
  bool _saveHists;
  int _nSaveHist;
  TString _hnn;
  bool _verbose;

};

#endif
