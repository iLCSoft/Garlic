/*********************************************************************
 * ECALGarlic
 *
 * GAmma Reconstruction for the LInear Collider:
 * A Marlin processor that clusters Calorimeter Hits and photon ID.
 *
 * 26 November 2007  Marcel Reinhard (LLR)
 *********************************************************************/

/*

V 1-1
-----

- corrected bug in Get2dProjDistance
- new method of seed finding
- introduced angled projection in the endcaps
- changed Core finding :
    1.) no search for first hit in cluster
    2.) take up to 4 cells per layer ( d < cellSize )
- don't add deeper hits when clustering: loss in speed but easier to control number of iterations needed.
- neighbour definition extended to 3 following layers
- distance cut increased in first 1/3 of calo
- removed bug in theta projection
- increased distance for additional hits to < sqrt(8)
- removed 2 bugs in overstepping the gap
- introduced function to obtain Energy in GeV from hitcount or energy measurement


V 1-2
-----

- gap correction revised
- using cell_index to be more precise in neighbour criterion
- Included EndcapRing completely
- exclude preShower hits from energy measurement and gap correction
- Array of cluster Parameters now filled with empty parameters to avoid event overlaps


V 1-3
-----

- dedicated to single particle studies (gamma, pi, tau)
- satellite merging if incosistent with pi0 mass  -- REMOVED
- ClusterParameters now part of ExtendedCluster
- Extend to more than one ROI but Histogram only biggest
- Gap Correction for alveola only deleted
- Using full chain of Digi and Tracking
- Moved to XML steering file format
- Take events with tracks
- Include MC Event info
- include PreShower hits seperately, and count them as cluster hits
- since in chain with LDCCaloDigi, Etot is already weighted and in GeV, thus omitting Etot_w, introducing Etot_g and E_GeV
- can read tracks via file, specify CheatTracks 2 in steering

V 1-4
-----
Di-Jet analysis

V 1-5
-----
Modified Gap coroection
Theta-dependent clustering distance

V 2-0 
-----
tidy up, rewrite...(Daniel Jeans)

V 2-1
-----
further tidy up & rewrite...(Daniel Jeans)
new NN variables and trainings

*/


#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCRelation.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>

#include <marlin/Global.h>
#include <marlin/ProcessorMgr.h>

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <gear/CalorimeterParameters.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/LayerLayout.h>

//Root stuff
#include <TROOT.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TString.h>

#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryParameters.hh"

#include "ECALGarlicCluster.hh"

#include "ECALGarlicExtendedCluster.hh"
#include "ECALGarlicExtendedHit.hh"

#include <HelixClass.h>
#include <ECALGarlic.hh>

using namespace marlin;
using namespace lcio;
using namespace std;

ECALGarlic anECALGarlic;

ECALGarlic::ECALGarlic() : Processor("ECALGarlic") {

  streamlog_out ( DEBUG1 ) << "hello from ECALGarlic constructor " << this << std::endl;

  // processor description
  _description = "Clustering and photon recognition";


  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollection" , 
                           "Name of the MCParticle input collection"  ,
			   _mcParticleCollectionName,
                           std::string("MCParticle") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "simHitCaloHitRelations",
			   "name of sim->calo hit relation collection",
			   _simHitCaloHitRelationCollectionName,
			   std::string("RelationCaloHit") );
			   
  // input pre-clusters
  registerInputCollection( LCIO::CLUSTER,
                           "EcalPreClusterCollection" , 
			   "Name of PreCluster collection",
			   _ecalPreClusterCollectionName,
			   std::string("ECAL_PreClusters") );

  // input track collection
  registerInputCollection( LCIO::TRACK,
			   "LDCTrackCollection",
			   "LDC track collection name",
			   _LDCTrackCollectionName,
			   std::string("LDCTracks") );

  registerInputCollection( LCIO::TRACK,
			   "TPCTrackCollection",
			   "TPC track collection name",
			   _TPCTrackCollectionName,
			   std::string("TPCTracks") );


  // output collections

  registerOutputCollection( LCIO::CALORIMETERHIT,
			    "SeedCollName",
			    "Name of garlic seed collection",
			    _seedCollName,
			    std::string("GARLICSeeds") );

  registerOutputCollection( LCIO::CALORIMETERHIT,
			    "TrackExtrapolationCollName",
			    "Name of track extrapolation collection",
			    _trkExtrapCollName,
			    std::string("GARLICTrackExtrapolations") );

  registerOutputCollection( LCIO::CLUSTER,
			    "ElectronClusterCollName",
			    "Name of electron cluster collection",
			    _electronCollName,
			    std::string("GARLICElectronClusters") );

  registerOutputCollection( LCIO::CLUSTER,
			    "CoreCollName",
			    "Name of garlic core collection",
			    _coreCollName,
			    std::string("GARLICCores") );

  registerOutputCollection( LCIO::CLUSTER,
			    "ClusterCollName",
			    "Name of garlic cluster collection",
			    _clusterCollName,
			    std::string("GARLICClusters") );

  registerOutputCollection( LCIO::LCGENERICOBJECT,
			    "ClusterParametersCollName",
			    "Name of garlic cluster parameters collection",
			    _clparsCollName,
			    std::string("GARLICClusterParameters") );

  registerOutputCollection( LCIO::LCRELATION,
			    "ClusterParLinksCollName",
			    "Name of cluster to parameter relations",
			    _clusterParRelCollName,
			    std::string("GARLICClusterParameterLinks") );

  registerOutputCollection( LCIO::LCRELATION,
			    "SeedCoreLinksCollName",
			    "Name of seed to core relations",
			    _seedCoreRelCollName,
			    std::string("GARLICSeedCoreLinks") );

  registerOutputCollection( LCIO::LCRELATION,
			    "SeedClusterLinksCollName",
			    "Name of seed to cluster relations",
			    _seedClusterRelCollName,
			    std::string("GARLICSeedClusterLinks") );

  registerOutputCollection( LCIO::CALORIMETERHIT,
			    "RemovedHitsCollection",
			    "collection name of removed hits (near tracks)",
			    _removedHitsCollectionName,
			    std::string("GARLICRemovedHits") );

  // parameters

  registerProcessorParameter("TrackCheat",
			     "take MC info for tracks?",
			     _x_cheatTracks,
			     false);

  registerProcessorParameter("TrackRemoveNearbyHits",
			     "Should remove Hits near extrapolated tracks to reject pions?",
			     _x_removeHitsNearTracks,
			     true);

  registerProcessorParameter("TrackVetoWindow",
			     "window around track in which to remove hits",
			     _x_trackWindowVeto, 
			     float (10.) );


  // seeding
  registerProcessorParameter("SeedNLayers",
			     "Number of ECAL pseudo layers used for projecting to obtain a seed, typically equivalent to 5 X0.",
			     _x_nLayersForSeeding,
			     int(8));

  registerProcessorParameter("SeedMinHits",
			     "Minimum number of hits to accept a seed.",
			     _x_minHitsForSeeding,
			     int(4));

  registerProcessorParameter("SeedMinHitEnergy",
			     "consider only hits above this threshold (in MIPs) in the seeding",
			     _x_seedHitEnergyCut,
			     float(3.5) );

  registerProcessorParameter("SeedMinEnergy",
			     "reject seeds below this energy (in MIPs)",
			     _x_seedEnergyCut,
			     float (10.) );

  registerProcessorParameter("SeedMaxDistance",
			     "radius of seed cylinder (in mm)",
			     _x_seedDistanceCut,
			     float (15.) );

  // core building
  registerProcessorParameter("CoreNlayersSection1",
			     "first calo section definition for core building",
			     _x_nlayersSection1,
			     int (16) );

  registerProcessorParameter("CoreMaxHoleSection1",
			     "max hole in core building for section 1",
			     _x_maxHoleSection1,
			     int (3) );

  registerProcessorParameter("CoreMaxHoleSection2",
			     "max hole in core building for section 2",
			     _x_maxHoleSection2,
			     int (2) );

  // clustering
  registerProcessorParameter("ClusterNIterations",
			     "Number of Iterations to apply the neighbouring criterion.",
			     _x_nIterations,
			     int(3));

  // cluster merging
  registerProcessorParameter("MergeTouchFrac",
			     "what fraction of layers must have adjacent hits in order to merge 2 clusters",
			     _x_mergeTouchFrac,
			     float (0.7) );
			     
  registerProcessorParameter("MergeInitSepVeto",
			     "veto merge of 2 clusters if first N layers are separated",
			     _x_initialLayerSeparation,
			     int (3) );

  // photon ID
  registerProcessorParameter("RejectClustersByMLPCut",
			     "Should remove clusters that do not pass MLP cuts?",
			     _x_rejectMLP,
			     true);


  _x_mlpCuts.clear();

  _x_mlpCuts.push_back(0.5);
  _x_mlpCuts.push_back(0.5);
  _x_mlpCuts.push_back(0.5); 
  _x_mlpCuts.push_back(0.5);
  _x_mlpCuts.push_back(0.5);
  registerProcessorParameter("MLPCuts",
			     "Cuts applied on the NN decision",
			     _x_mlpCuts,
			     _x_mlpCuts,
			     _x_mlpCuts.size() 
			     );

  registerProcessorParameter("DebugMode",
			     "Talk a lot? (0-3)",
			     _x_debug,
			     int(0));

  // detector parameters
  // should this not be taken from the gear file???

  registerProcessorParameter("SiliconThickness",
			     "Thickness of Silicon defined as layerThickness - absorberThickness - ThicknessOfPassiveMaterial(PCB etc.)",
			     _x_activeThickness,
			     float(1.0));

  registerProcessorParameter("PCBThickness",
			     "Thickness of PassiveMaterial(PCB etc.) in a layer that is not the absorber",
			     _x_passiveThickness,
			     float(2.15));

  registerProcessorParameter("FirstEndcapLayerOffset",
			     "Offset ECAL Layers to match reconstructed hits",
			     _x_firstEndcapLayerOffset,
			     float(-0.4));

  registerProcessorParameter("FirstBarrelLayerOffset",
			     "Offset ECAL Layers to match reconstructed hits",
			     _x_firstBarrelLayerOffset,
			     float(-5.14));

  registerProcessorParameter("clusterCheckHistoFile", "name of file in which to save clustering histograms",
                             _histFileName, std::string(""));

  _clusterer=NULL;
  _fhistos=NULL;
  _track_extrap_col=NULL;
  _mcParticleColl=NULL;

}


ECALGarlic::~ECALGarlic() {
  streamlog_out ( DEBUG1 ) << "hello from ECALGarlic destructor " << this << std::endl;

  if (_clusterer) {delete _clusterer; _clusterer=NULL;}
  if (_fhistos) {delete _fhistos; _fhistos=NULL;}

}

void ECALGarlic::init() {
  streamlog_out ( DEBUG1 ) << "hello from ECALGarlic init() " << this << std::endl;

}

void ECALGarlic::setup()
{

  streamlog_out ( DEBUG1 ) << "hello from ECALGarlic::setup() " << std::endl;

  printMrGarlic();
  printParameters();

  
  ECALGarlicAlgorithmParameters::Instance().SetEcalPreClusterCollectionName (_ecalPreClusterCollectionName);

  ECALGarlicAlgorithmParameters::Instance().SetDebug (_x_debug);

  ECALGarlicAlgorithmParameters::Instance().SetTrackCheat     (_x_cheatTracks);
  ECALGarlicAlgorithmParameters::Instance().SetTrackRemoveNearbyHits(_x_removeHitsNearTracks);
  ECALGarlicAlgorithmParameters::Instance().SetTrackVetoWindow (_x_trackWindowVeto);

  ECALGarlicAlgorithmParameters::Instance().SetSeedNLayers (_x_nLayersForSeeding);
  ECALGarlicAlgorithmParameters::Instance().SetSeedMinHits(_x_minHitsForSeeding);
  ECALGarlicAlgorithmParameters::Instance().SetSeedHitEnergyCut(_x_seedHitEnergyCut);
  ECALGarlicAlgorithmParameters::Instance().SetSeedEnergyCut   (_x_seedEnergyCut);
  ECALGarlicAlgorithmParameters::Instance().SetSeedDistanceCut (_x_seedDistanceCut);

  ECALGarlicAlgorithmParameters::Instance().SetClusterNIterations(_x_nIterations);

  ECALGarlicAlgorithmParameters::Instance().SetRejectMLP(_x_rejectMLP);
  ECALGarlicAlgorithmParameters::Instance().SetMLPCuts(_x_mlpCuts);

  ECALGarlicAlgorithmParameters::Instance().SetCoreLayersSection1(_x_nlayersSection1);
  ECALGarlicAlgorithmParameters::Instance().SetCoreMaxHoleSection1(_x_maxHoleSection1);
  ECALGarlicAlgorithmParameters::Instance().SetCoreMaxHoleSection2(_x_maxHoleSection2);

  ECALGarlicAlgorithmParameters::Instance().SetMergeTouchFraction(_x_mergeTouchFrac);
  ECALGarlicAlgorithmParameters::Instance().SetMergeInitalLayerSeparation(_x_initialLayerSeparation);


  _clusterer = new ECALGarlicCluster();

  _nEvents=0;
  _nPhotonClusters=0;


  if (_histFileName!="")
    _fhistos = new TFile(_histFileName.c_str(),"recreate");


  _nSaveHist=0;
  _geomSetup = false;

  return;
}



void ECALGarlic::processRunHeader(LCRunHeader * run)
{
  return;
}


void ECALGarlic::processEvent(LCEvent * evt)   // main !
{
  streamlog_out ( DEBUG1 ) << endl << "Event: " << evt->getEventNumber() << endl;
  if(evt->getEventNumber()%100==0)
    streamlog_out ( DEBUG1 ) << endl << "Event: " << evt->getEventNumber() << endl;

  // set up the detector geometry
  if (!_geomSetup) {
    setUpGeometry();
    _geomSetup=true;
  }

  _nEvents++;
  
  _mcParticleColl = 0;
  try {
    _mcParticleColl = evt->getCollection(_mcParticleCollectionName);
  }
  catch(DataNotAvailableException err) {};

  LCCollection* _caloHitRelationColl = 0;
  try {
    _caloHitRelationColl = evt->getCollection(_simHitCaloHitRelationCollectionName);
  } 
  catch(DataNotAvailableException err) {};

  LCCollectionVec* coreClusterColl=0;
  coreClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  coreClusterColl->setFlag( coreClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  LCCollectionVec* finalClusterColl=0;
  LCCollectionVec* finalSelectedClusterColl=0;
  LCCollectionVec* finalRejectedClusterColl=0;
  LCCollectionVec* finalClusterCollParameters=0;
  LCCollectionVec* finalCluParLinkColl=0;
  LCCollectionVec* electronClusterColl=0;

  LCCollectionVec* MCphotonClusterColl=0;
  LCCollectionVec* MCelectronClusterColl=0;
  LCCollectionVec* MCchHadClusterColl=0;
  LCCollectionVec* MCneuHadClusterColl=0;
  LCCollectionVec* MCotherClusterColl=0;

  LCCollectionVec* wronglySelectedClusterColl=0;
  LCCollectionVec* wronglyRejectedClusterColl=0;

  finalClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalClusterColl->setFlag( finalClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  finalSelectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalSelectedClusterColl->setFlag( finalSelectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalSelectedClusterColl->setSubset();

  finalRejectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalRejectedClusterColl->setFlag( finalRejectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalRejectedClusterColl->setSubset();

  finalClusterCollParameters = new LCCollectionVec(LCIO::LCGENERICOBJECT);
  finalCluParLinkColl = new LCCollectionVec(LCIO::LCRELATION);

  electronClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  electronClusterColl->setFlag( electronClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  MCphotonClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  MCphotonClusterColl->setFlag( MCphotonClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  MCphotonClusterColl->setSubset();

  MCelectronClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  MCelectronClusterColl->setFlag( MCelectronClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  MCelectronClusterColl->setSubset();

  MCchHadClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  MCchHadClusterColl->setFlag( MCchHadClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  MCchHadClusterColl->setSubset();
  
  MCneuHadClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  MCneuHadClusterColl->setFlag( MCneuHadClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  MCneuHadClusterColl->setSubset();

  MCotherClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  MCotherClusterColl->setFlag( MCotherClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  MCotherClusterColl->setSubset();

  wronglySelectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  wronglySelectedClusterColl->setFlag(wronglySelectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  wronglySelectedClusterColl->setSubset();
  
  wronglyRejectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  wronglyRejectedClusterColl->setFlag(wronglyRejectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  wronglyRejectedClusterColl->setSubset();

  IMPL::LCCollectionVec *seed_col = 0;
  _track_extrap_col = 0;
  seed_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  seed_col->setFlag(seed_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

  IMPL::LCCollectionVec *rej_seed_col = 0;
  _track_extrap_col = 0;
  rej_seed_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  rej_seed_col->setFlag(rej_seed_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

  _track_extrap_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  _track_extrap_col->setFlag(_track_extrap_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));


  IMPL::LCCollectionVec *seedCoreLink=0;
  seedCoreLink = new LCCollectionVec(LCIO::LCRELATION);

  IMPL::LCCollectionVec *seedFinalClusLink=0;
  seedFinalClusLink = new LCCollectionVec(LCIO::LCRELATION);

  // prepare the Track collection
  vector <ExtendedTrack* > trackVec;
  vector <ExtendedCluster2* > preClusVec;

  if ( ECALGarlicAlgorithmParameters::Instance().GetTrackCheat() == 0 )
    PrepareTracks(evt,trackVec);
  else 
    PrepareMCTracks(trackVec);

  // 1.) transform PreClusters to ExtendedClusters (pseudoIDs etc)
  streamlog_out ( DEBUG2 ) << "Preparing PreClusters... " << endl;
  PreparePreClusters(evt, preClusVec);
  streamlog_out ( DEBUG2 ) << preClusVec.size() << " PreClusters prepared" << endl;

  if(preClusVec.size()==0) {
    streamlog_out ( MESSAGE ) << "no preclusters found!" << endl;
    return;
  }

  // find the tracks near each precluster
  std::map < ExtendedCluster2*, vector <ExtendedTrack* > > nearbyTracks;
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {
    nearbyTracks[ preClusVec[roi_i] ] = selectNearbyTracks( preClusVec[roi_i], &trackVec);
  }    

  // look for electrons (seeded by tracks)
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {
    vector < pair < ExtendedTrack*, ExtendedCluster2* > > electrons = _clusterer->getElectrons(preClusVec[roi_i], nearbyTracks[ preClusVec[roi_i] ]);    
    // add to event
    if (electronClusterColl) {
      for (size_t ic=0; ic<electrons.size(); ic++) {
	ClusterImpl* climp = new ClusterImpl(electrons[ic].second->getClusterImpl());
	electronClusterColl->addElement(climp);
      }
    }
    RemoveElectronHits( electrons, preClusVec[roi_i] );

    for (size_t i=0; i<electrons.size(); i++) {
      if (electrons[i].second) {
	delete electrons[i].second;
	electrons[i].second=NULL;
      }
    }

  }


  // 1b.) clear hits near to extrapolated tracks
  if(ECALGarlicAlgorithmParameters::Instance().GetTrackRemoveNearbyHits()) {
    streamlog_out ( DEBUG2 ) << "Removing hits near tracks..." << endl;
    RemoveHitsNearExtrapolatedTracks(evt,trackVec,preClusVec);
  }

  // the next steps are done per ROI (preCluster)
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {


    ExtendedCluster2 *preCluster=preClusVec[roi_i];

    streamlog_out ( DEBUG2 ) 
      << "evt" << evt->getEventNumber() << " roi " << roi_i << " nhits=" << preCluster->getHits()->size() << endl
      << "Building histogram for seeding for PreCluster " << roi_i << endl;

    if (preCluster->getHits()->size()==0) {
      streamlog_out ( DEBUG2 ) << "precluster/RoI with no hits, ignoring!" << endl;
      continue;
    }

    // get tracks close to this precluster
    //    vector <ExtendedTrack* > nearbyTracks = selectNearbyTracks(preCluster, &trackVec);

    // get the seeds
    if (_fhistos && _nSaveHist++<MAXSAVEHIST) {
      TString hn = "evt"; hn+=evt->getEventNumber();
      hn+="_roi"; hn+=roi_i;
      _clusterer->saveHistos(_fhistos, hn);
    } else _clusterer->saveHistos(NULL);

    //    std::vector <CalorimeterHit*> seeds = _clusterer->getSeeds(preCluster);
    std::map <CalorimeterHit*, bool> allseeds = _clusterer->getSeeds(preCluster);
    // add seeds to collection
    std::vector <CalorimeterHit*> goodseeds;
    for ( std::map <CalorimeterHit*, bool>::iterator ijs=allseeds.begin(); ijs!=allseeds.end(); ijs++ ) {
      if (ijs->second) {
	goodseeds.push_back(ijs->first);
	seed_col->addElement(ijs->first);
      } else {
	rej_seed_col->addElement(ijs->first);
      }
    }

    // for (size_t is=0; is<seeds.size(); is++) seed_col->addElement(seeds[is]);
    // streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " nseeds = " << seeds.size() << endl;
    // for (size_t is=0; is<seeds.size(); is++) 
    //   streamlog_out ( DEBUG2 ) << "    " << seeds[is]->getPosition()[0] << " " << seeds[is]->getPosition()[1] << " " << seeds[is]->getPosition()[2] << endl;

    if (goodseeds.size()==0) continue;


    // build the cores
    std::map < CalorimeterHit*, ExtendedCluster2* > cores =  _clusterer->getCores(preCluster, goodseeds);
    streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " ncores = " << cores.size() << endl;

    for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt=cores.begin(); itt!=cores.end(); itt++) {
      ClusterImpl* climp = new ClusterImpl(itt->second->getClusterImpl());
      coreClusterColl->addElement(climp);
      streamlog_out ( DEBUG2 ) << "   saving core with nhits = " << climp->getCalorimeterHits().size() << endl;
      LCRelationImpl* rel = new LCRelationImpl(itt->first, climp);
      seedCoreLink->addElement(rel);
    }

    // build the clusters
    std::map < CalorimeterHit*, ExtendedCluster2* > Cclusters = _clusterer->getClusters(preCluster, cores);
    streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " nclusters = " << Cclusters.size() << endl;

    // merge clusters
    _clusterer->mergeSatellites(Cclusters);
    streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " nclusters (after merge) = " << Cclusters.size() << endl;

    // add tracks to clusters
    for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt=Cclusters.begin(); itt!=Cclusters.end(); itt++) {
      itt->second->setTracks( & nearbyTracks[ preClusVec[roi_i] ] );
    }

    for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt=Cclusters.begin(); itt!=Cclusters.end(); itt++) {

      ClusterImpl* climp = new ClusterImpl(itt->second->getClusterImpl());
      finalClusterColl->addElement(climp);

      bool selected = itt->second->getNNsel();

      if (selected)
	finalSelectedClusterColl->addElement(climp);
      else 
	finalRejectedClusterColl->addElement(climp);

      IMPL::LCGenericObjectImpl* genobj = new IMPL::LCGenericObjectImpl(itt->second->getGenericObject());
      finalClusterCollParameters->addElement(genobj);

      LCRelationImpl* rel = new LCRelationImpl(climp, genobj);
      finalCluParLinkColl->addElement(rel);

      rel = new LCRelationImpl(itt->first, climp);
      seedFinalClusLink->addElement(rel);

      itt->second->setHitRelationColl(_caloHitRelationColl);
      int pdg = itt->second->getMCCalPdg()[0]; // dominant pdg
      //	float frac = itt->second->getMCCalFrac()[0]; // the fraction of energy from this pdg

      TString PDG;

      switch (pdg) {
      case (22):
	MCphotonClusterColl->addElement(climp);
	PDG="photon";
	break;
      case (11):
	MCelectronClusterColl->addElement(climp);
	PDG="electron";
	break;
      case (211):
      case (321):
      case (2212):
	PDG="chHad";
	MCchHadClusterColl->addElement(climp);
	break;
      case (130):
      case (310):
      case (2112):
	PDG="neuHad";
	MCneuHadClusterColl->addElement(climp);
	break;
      default:
	PDG="other";
	MCotherClusterColl->addElement(climp);
	if (pdg!=-999) 
	  streamlog_out( DEBUG ) << "undefined pdg " << pdg << endl;
	break;
      }
      
      if (selected && pdg!=22) wronglySelectedClusterColl->addElement(climp);
      else if (!selected && pdg==22) wronglyRejectedClusterColl->addElement(climp);

      // the seed has been added to a collection, as has the cluster.
      // however, we need to delete the extendedcluster
      delete itt->second; itt->second=NULL;

    }
  }

  addCollToEvent(seed_col, _seedCollName, evt);
  addCollToEvent(rej_seed_col, _seedCollName+"Rejected", evt);
  addCollToEvent(_track_extrap_col, _trkExtrapCollName, evt);
  addCollToEvent(coreClusterColl, _coreCollName, evt);
  addCollToEvent(finalClusterColl, _clusterCollName, evt);
  addCollToEvent(finalSelectedClusterColl, _clusterCollName+"Selected", evt);
  addCollToEvent(finalRejectedClusterColl, _clusterCollName+"Rejected", evt);
  addCollToEvent(finalClusterCollParameters, _clparsCollName, evt);
  addCollToEvent(finalCluParLinkColl, _clusterParRelCollName, evt);
  addCollToEvent(electronClusterColl, _electronCollName, evt);
  addCollToEvent(MCphotonClusterColl, _clusterCollName+"MCphoton", evt);
  addCollToEvent(MCelectronClusterColl, _clusterCollName+"MCelectron", evt);
  addCollToEvent(MCchHadClusterColl, _clusterCollName+"MCchHad", evt);
  addCollToEvent(MCneuHadClusterColl, _clusterCollName+"MCneuHad", evt);
  addCollToEvent(MCotherClusterColl, _clusterCollName+"MCother", evt);
  addCollToEvent(wronglyRejectedClusterColl, _clusterCollName+"WronglyRejected", evt);
  addCollToEvent(wronglySelectedClusterColl, _clusterCollName+"WronglySelected", evt);
  addCollToEvent(seedCoreLink, _seedCoreRelCollName, evt);
  addCollToEvent(seedFinalClusLink, _seedClusterRelCollName, evt);

  cleanup(preClusVec);
  cleanup(trackVec);

  streamlog_out ( DEBUG1 ) << "Finished Event " << evt->getEventNumber() << endl;

  return;  
}


void ECALGarlic::addCollToEvent(LCCollection* col, string colname, LCEvent* evt) {
  if (col) {
    if (evt && colname.find("IGNORE")==string::npos) {
      evt->addCollection(col, colname);
    } else {
      delete col;
    }
  }
  return;
}

void ECALGarlic::cleanup(vector<ExtendedCluster2* > &preClusVec) {
  if(preClusVec.size()>0) {
    for (size_t i = 0; i < preClusVec.size(); i++) {
      ExtendedCluster2 *a_ext_preCluster = preClusVec[i];
      if (a_ext_preCluster) {
	vector<ExtendedHit2* > ext_preHit_vec = *(a_ext_preCluster->getHits());
	for (size_t pHit_i=0; pHit_i<ext_preHit_vec.size(); pHit_i++) {
	  ExtendedHit2 *a_ext_pHit = dynamic_cast<ExtendedHit2*>(ext_preHit_vec[pHit_i]);
	  if (a_ext_pHit) {
	    delete a_ext_pHit;
	    a_ext_pHit=NULL;
	  }
	}
	delete a_ext_preCluster;
	a_ext_preCluster=NULL;
      }
    }
  }
  return;
}

void ECALGarlic::cleanup(vector<ExtendedTrack* > &trackVec) {
  for (size_t i=0; i<trackVec.size(); i++) 
    if (trackVec[i]) {delete trackVec[i]; trackVec[i]=NULL;}
  trackVec.clear();
  streamlog_out ( DEBUG2 ) << "Cleanup finished" << endl;
}



void ECALGarlic::check(LCEvent * evt)
{
}



void ECALGarlic::end()
{
  streamlog_out ( MESSAGE ) << endl
			    << "ECALGarlic Report: " << endl
			    << _nEvents << " events processed" << endl
			    << "Found " << _nClusters << " Cluster" << endl;
  if(_nEvents>0)
    streamlog_out ( MESSAGE ) << "= " << ((float)_nClusters)/_nEvents  << " /event" << endl;

  if (_fhistos) {
    _fhistos->Write(0);
    _fhistos->Close();
  }

}



void ECALGarlic::PrepareMCTracks(vector<ExtendedTrack* > &trackVec)
{
  streamlog_out ( DEBUG2 ) << "Getting MC tracks..." << endl
			   << "BField: " << ECALGarlicGeometryParameters::Instance().Get_bField() << " T" << endl;

  if (!_mcParticleColl) {
    streamlog_out ( DEBUG1 ) << "no MC collection found..." << endl;
    return;
  }

  for(int jmcp=0; jmcp<_mcParticleColl->getNumberOfElements(); jmcp++) {
    MCParticle* mcp = dynamic_cast <MCParticle*> (_mcParticleColl->getElementAt(jmcp));
    if (mcp && mcp->getCharge()!=0 && !mcp->isDecayedInTracker()) {
      ExtendedTrack* trk = new ExtendedTrack(mcp);
      trackVec.push_back(trk);
    }
  }
  streamlog_out ( DEBUG2 ) << "Have prepared " << trackVec.size() << " MC tracks" << endl;
}


void ECALGarlic::PrepareTracks(const LCEvent *evt, vector<ExtendedTrack* > &trackVec)
{
  streamlog_out ( DEBUG2 ) << "Getting tracks..." << endl
			   << "BField: " << ECALGarlicGeometryParameters::Instance().Get_bField() << " T" << endl;

  TrackVec allSubTracks;

  try {
    LCCollection* trackColl = evt->getCollection(_LDCTrackCollectionName);
    streamlog_out ( DEBUG2 )  << trackColl->getNumberOfElements() << " tracks in LDC track collection" << endl;
    for(int track_i=0; track_i<trackColl->getNumberOfElements(); track_i++) {
      streamlog_out ( DEBUG2 ) << "Track " << track_i << endl;
      Track *a_track = dynamic_cast<Track*>(trackColl->getElementAt(track_i));

      TrackVec subtracks = a_track->getTracks();
      for (size_t f=0; f<subtracks.size(); f++) {
	allSubTracks.push_back(subtracks[f]);
      }

      ExtendedTrack *a_ext_track = new ExtendedTrack(a_track);
      trackVec.push_back(a_ext_track);
      if (_track_extrap_col) {
	CalorimeterHitImpl *extrap=new CalorimeterHitImpl();
	extrap->setEnergy(1000+track_i);
	extrap->setPosition( a_ext_track->getEcalEntryPos() );
	_track_extrap_col->addElement(extrap);
      }
    }
  } catch (DataNotAvailableException err) {};

  // also try TPC tracks
  try {
    LCCollection* trackColl = evt->getCollection(_TPCTrackCollectionName);
    streamlog_out ( DEBUG2 )  << trackColl->getNumberOfElements() << " tracks in TPC track collection" << endl;
    for(int track_i=0; track_i<trackColl->getNumberOfElements(); track_i++) {
      streamlog_out ( DEBUG2 ) << "Track " << track_i << endl;
      Track *a_track = dynamic_cast<Track*>(trackColl->getElementAt(track_i));

      // This TPC track already considered as part of LDC track, skip it
      if (find ( allSubTracks.begin(), allSubTracks.end(), a_track )!=allSubTracks.end()) continue;

      ExtendedTrack *a_ext_track = new ExtendedTrack(a_track);
      trackVec.push_back(a_ext_track);
      if (_track_extrap_col) {
	CalorimeterHitImpl *extrap=new CalorimeterHitImpl();
	extrap->setEnergy(1000+track_i);
	extrap->setPosition( a_ext_track->getEcalEntryPos() );
	_track_extrap_col->addElement(extrap);
      }
    }
  } catch (DataNotAvailableException err) {};



  streamlog_out ( DEBUG2 ) << "Have " << trackVec.size() << " tracks" << endl;
}



void ECALGarlic::PreparePreClusters(LCEvent *evt, vector<ExtendedCluster2* > &preClusVec) {

  LCCollection *preClusColl = 0;
  try {
    preClusColl = evt->getCollection(ECALGarlicAlgorithmParameters::Instance().GetEcalPreClusterCollectionName());
  } catch(DataNotAvailableException &exc) {
    if (ECALGarlicAlgorithmParameters::Instance().GetDebug()>1)
      std::cerr << "In ECALGarlic::preparePreClusters "	<< exc.what() << std::endl;
    return;
  }

  int NPreCluster = preClusColl->getNumberOfElements();
  if (!ECALGarlicGeometryParameters::Instance().Get_defaultDecoder() && NPreCluster>0) {
    CellIDDecoder<CalorimeterHit>* dec = new CellIDDecoder<CalorimeterHit> (preClusColl);
    ECALGarlicGeometryParameters::Instance().Set_defaultDecoder(dec);
  }
  streamlog_out ( DEBUG2 ) << "Number of PreClusters in original collection: " << preClusColl->getNumberOfElements() << endl;

  for (int preC_i = 0; preC_i < NPreCluster; preC_i++) {
    ClusterImpl *a_cluster = dynamic_cast<ClusterImpl* >( preClusColl->getElementAt(preC_i) );

    ExtendedCluster2 *a_ext_cluster = new ExtendedCluster2( );
    const CalorimeterHitVec &preHitVec=a_cluster->getCalorimeterHits();
    streamlog_out ( DEBUG2 ) << "Reading PreCluster " << preC_i << " " << preHitVec.size() << endl;
    vector<ExtendedHit2 *> extHitVec;

    for (EVENT::CalorimeterHitVec::const_iterator hit_it=preHitVec.begin(); hit_it!=preHitVec.end(); hit_it++) {
      CalorimeterHitImpl *a_hit = dynamic_cast<CalorimeterHitImpl*>( *hit_it );
      ExtendedHit2 *a_ext_hit = new ExtendedHit2(a_hit);
      a_ext_hit->setCluster(NULL);
      a_ext_hit->setPreshower(0);
      extHitVec.push_back(a_ext_hit);      
    }

    // sort hits by pseudo layer 
    sort(extHitVec.begin(),extHitVec.end(),ExtendedHit2::lowerPseudoLayer);
    a_ext_cluster->addHits(extHitVec);
    preClusVec.push_back(a_ext_cluster);
  }

  return;
}

void ECALGarlic::RemoveElectronHits(vector < pair < ExtendedTrack*, ExtendedCluster2* > > electrons, ExtendedCluster2* preClus ) {

  for (size_t iel=0; iel<electrons.size(); iel++) {
    std::vector <ExtendedHit2*> * electronHits = electrons[iel].second->getHits();
    
    for (size_t ih=0; ih<electronHits->size(); ih++) {
      preClus->removeHit( electronHits->at(ih) );
    }

  }
  return;
}


void ECALGarlic::RemoveHitsNearExtrapolatedTracks(LCEvent *evt,
						  vector<ExtendedTrack* > &trackVec, 
						  vector<ExtendedCluster2*> &preClusVec) {

  IMPL::LCCollectionVec *removedHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  removedHits->setFlag(removedHits->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  removedHits->setSubset();

  float min_dist_cut = ECALGarlicAlgorithmParameters::Instance().GetTrackVetoWindow();
  float dist_cut = max( ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1], min_dist_cut);

  for (size_t track_i=0; track_i<trackVec.size(); track_i++) {

    streamlog_out ( DEBUG2 ) << " Checking track " << track_i << " : " ;

    ExtendedTrack *a_ext_track = dynamic_cast<ExtendedTrack* >(trackVec[track_i]);

    // check that extrapolation to ecal looks reasonably well controlled

    // int nhits = a_ext_track->getTrack()->getTrackerHits().size();
    // TrackerHit* lastHit =  a_ext_track->getTrack()->getTrackerHits()[nhits-1];

    if (a_ext_track->getMinDist_HitEcalEntry() > 400) {
      continue;
    }

    for(size_t pclus_i=0; pclus_i<preClusVec.size(); pclus_i++) {

      streamlog_out ( DEBUG2 ) << "Removing hits from PreCluster " << pclus_i << endl;

      vector<ExtendedHit2* > *hit_vec = preClusVec[pclus_i]->getHits();

      int nhits = hit_vec->size();

      for (int ih=nhits-1; ih>=0; ih--) {
	ExtendedHit2 *a_ext_hit = dynamic_cast<ExtendedHit2* > (hit_vec->at(ih));
	float dist = a_ext_track->getDistanceToPoint( a_ext_hit->getCaloHit()->getPosition() );
	if(dist<dist_cut) {       
	  hit_vec->erase(hit_vec->begin()+ih);
	  removedHits->addElement( dynamic_cast <CalorimeterHitImpl*> (a_ext_hit->getCaloHit()) );
	} else {
	  streamlog_out ( DEBUG2 ) << "Hit not removed: " << a_ext_hit->getCaloHit()->getCellID0() << 
	    ", with dist = " << dist << " while asked: "<< dist_cut <<  endl;
	}
      }
      streamlog_out ( DEBUG2 ) << "Done with PreCluster " << pclus_i << endl;
    }
  }

  evt->addCollection(removedHits,_removedHitsCollectionName);
}


vector < ExtendedTrack* > ECALGarlic::selectNearbyTracks(ExtendedCluster2* preCluster,  vector <ExtendedTrack* >* trackVec) {
  float distCut = 10*preCluster->getWidth();
  vector < ExtendedTrack* > nearbytracks;
  for (size_t it=0; it<trackVec->size(); it++) {
    ExtendedTrack* trk = trackVec->at(it);
    float dist(0);
    for (int j=0; j<3; j++) {
      dist += pow ( trk->getEcalEntryPos()[j] - preCluster->getProjectedPosition()[j] , 2);
    }
    if (dist<distCut) nearbytracks.push_back(trk);
    else {
      dist=0;
      for (int j=0; j<3; j++)
	dist += pow ( trk->getEcalEntryPos()[j] - preCluster->getCentreOfGravity()[j] , 2);
      dist = sqrt(dist);
      if (dist<distCut) nearbytracks.push_back(trk);
    }
  }
  return nearbytracks;
}



void ECALGarlic::printMrGarlic() {

  streamlog_out ( MESSAGE ) 
    <<"                       `.+MY'WMa," << endl
    <<"                       .Mt.4&JnHTN," << endl
    <<"                       d@,o..Z,'`4WN," << endl
    <<"                       .WMHHM#m,  JyMx" << endl
    <<"     `     `    `       `    ,Mp4. j,M,  `        `   `    `   `    `" << endl
    <<"                             `J#.t ,|`N.     `" << endl
    <<"                    ` `  `` .d#!J`  S ZN. `" << endl
    <<"                `     ...gMY'`,=    .h.7N,      `" << endl
    <<"     `     `   ` ..gMB'7`  `         `?&,THNaJ,         `    `   `" << endl
    <<"          ` `..H#'!..v^                  ?4..`'WNa, ` `" << endl
    <<"       `` .JM9^  .J^ ..J^`         `  ?7i.. `G, ..TMm,  `             `" << endl
    <<"     `  .H#^ ` .J' .Y` ..v''3J.  ..v7=i.  .S. Jl``n.`WN," << endl
    <<"     `.H@`    .Z  .^  J'     ``nJ!      7,  ^  .L` T, .TN.`    `" << endl
    <<"     .H       J`     .bJ.......dG.......JR      ,|  ?,  ?N.         `" << endl
    <<"     d#      .%       4.dN#NJbCTYOd##bJSJt       4   4   M|7" << endl
    <<"     db `    ,!       `J4#W'  `  `  ?#WY! ..     `.  J   M       `" << endl
    <<"     JN. j   ,|         JJ%          ,.$  .Y``   ,:  Z  .#`" << endl
    <<"      Wb .h`  4.       ` ?S,     `  .Y= `.M:   ` J  `! .M      `" << endl
    <<"     `.Wh. 4.  S.`   `.J,   73J...v'`  .JMF    `.% ?``.Mt" << endl
    <<"      ` 7MJ ?o  4,       ?=&.......Jg#N#M^    `.% `  .M^" << endl
    <<"        `.TMa.7+.?o.`        ` ?7''YBB'!    ` .= `..M'  `    `      `" << endl
    <<"            ?WNadG,?a.`   `     `       ``  .Y...HB^" << endl
    <<"               .7'HMMMNJ,.` ``  `  `  ` ..gMHMB'!" << endl
    <<"     `               `.dMMMMHW+....JWMHMZMN,                     `" << endl
    <<"                  ` .d#TdNWWpWWNWWMkWHWHZ1ITN,`" << endl
    <<"               `  .d#OydNWt?YC!?I!?1?7dHR+zzvWb       `    `" << endl
    <<"             `   `J#Jd4HW#.zz+.+l:+z1:dWKvNOzvUb  `            `" << endl
    <<"     `       `.Mm.MtdDdHW%?!v!?1v:+v?1HW$,Jy.+.M|                   `" << endl
    <<"               dNdBJdxdpWr11z.++z+J1udHH1v.D1I:Jb" << endl
    <<"           `  `.dBnOVHdHWb..+.`?!`vjHpWtjaJO++.d%  `    `    `" << endl
    <<"     `       `.M^. ?7M##HHJvz:1z+x1HpWDzJRTmo:JN," << endl
    <<"              dh,  .Gd@;;??7TYY9VYTYYWH'! 4JzGM'                 `" << endl
    <<"               .M,.,d#;?;;;;;;;;;;;;+SJ.`  .B=`" << endl
    <<"       `   `     ?!dM1;;?;;;;;;;;;;;?+b..,l,r         `    `          `" << endl
    <<"              `   .M@;;?;??;??;++?;;?;+vYWM= `" << endl
    <<"     `       ..MMMMMmx;;;;;;;?jM#;?;;;;;;d@                    `" << endl
    <<"          `  dM9VOVTHdHJ?;?;jd#4Nx;?;;;;;?Hm.                       `" << endl
    <<"            .MZlltltltWHp;;jMF` WN+;?;??;;?JMe`       `    `" << endl
    <<"          ` ,MZllllltltd@c?d#    TNx;;;;?;jHM#    `              `" << endl
    <<"     `     ` WNyttlltwAdHIjM'    `.TMNkmWHUuXML `" << endl
    <<"              ?WMNUUZlltdMMF       JMSuuzuuuzM#                `" << endl
    <<"                .MNOlttOtd#        .WMNkQkWHNHD`        `" << endl
    <<"                  TMNgggg#^ `     `   ?''''!                 `" << endl
    <<"     `               ??;                                            `" << endl
    <<"                                                  `   `" << endl
    <<"           `                                               `" << endl
    <<"     `                                                           `" << endl
    <<"                 `                `" << endl
    <<"" << endl;

  return;
}


void ECALGarlic::setUpGeometry() {

  setup();

  // size of dead zones, ncells, si thickness, (not in GEAR file, specified in steering file...)
  ECALGarlicGeometryParameters::Instance().Set_activeThickness (_x_activeThickness);
  ECALGarlicGeometryParameters::Instance().Set_defaultDecoder  (NULL);

  // b field
  ECALGarlicGeometryParameters::Instance().Set_bField( Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z() );

  // Calorimeter geometry from GEAR
  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
  const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();

  ECALGarlicGeometryParameters::Instance().Set_symmetry ( pEcalBarrel.getSymmetryOrder() );

  std::vector< vector < float > > barrelStaveDir;
  if(pEcalBarrel.getSymmetryOrder() > 1){
    float phi0 = pEcalBarrel.getPhi0();
    for(int i=0; i<pEcalBarrel.getSymmetryOrder(); i++){
      float phi  = phi0 + i*2*acos(-1.)/pEcalBarrel.getSymmetryOrder();
      vector < float > staveDir;
      staveDir.push_back( -sin(phi) );
      staveDir.push_back( cos(phi) );
      staveDir.push_back( 0. );
      barrelStaveDir.push_back(staveDir);
    }
  }
  ECALGarlicGeometryParameters::Instance().Set_barrelStaveDir(barrelStaveDir);

  // radius of barrel front face
  ECALGarlicGeometryParameters::Instance().Set_rOfBarrel ( pEcalBarrel.getExtent()[0] );

  // z of barrel end
  ECALGarlicGeometryParameters::Instance().Set_zOfBarrel ( pEcalBarrel.getExtent()[3] );

  // inner r of endcap
  ECALGarlicGeometryParameters::Instance().Set_rInnerEcalEndcap ( pEcalEndcap.getExtent()[0] );

  // outer r of endcap
  ECALGarlicGeometryParameters::Instance().Set_rOfEndcap ( pEcalEndcap.getExtent()[1] );

  // z of endcap front face
  ECALGarlicGeometryParameters::Instance().Set_zOfEndcap ( pEcalEndcap.getExtent()[2] );

  // number of layers (include the preshower layer)
  ECALGarlicGeometryParameters::Instance().Set_nBarrelEcalLayers( ecalBarrelLayout.getNLayers()+1 );
  ECALGarlicGeometryParameters::Instance().Set_nEndcapEcalLayers( ecalEndcapLayout.getNLayers()+1 );
  ECALGarlicGeometryParameters::Instance().Set_nPseudoLayers( max(ecalBarrelLayout.getNLayers(), ecalEndcapLayout.getNLayers()) + 1 );

  // thickness of passive material (not absorber or silicon) (assumed constant)
  float totalThick = ecalBarrelLayout.getThickness(2);
  float absorbThick = ecalBarrelLayout.getAbsorberThickness(2);
  float passiveThick = totalThick-absorbThick-_x_activeThickness; // all that is not absorber or silicon
  ECALGarlicGeometryParameters::Instance().Set_passiveThickness(passiveThick);

  // layer0 is "preshower" layer
  // layer1 is "first" ECAL layer

  // absorber thickness arrays in barrel and endcap

  float absThick[MAX_NUMBER_OF_LAYERS]={0};
  float cellSize[MAX_NUMBER_OF_LAYERS]={0};
  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    absThick[i+1] = ecalBarrelLayout.getAbsorberThickness(i);
    cellSize[i+1] = ecalBarrelLayout.getCellSize0(i);
  }
  absThick[0]=0;
  cellSize[0]=cellSize[1];
  ECALGarlicGeometryParameters::Instance().Set_absThicknessBarrelLayer(absThick);
  ECALGarlicGeometryParameters::Instance().Set_padSizeEcal(cellSize);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) absThick[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    absThick[i+1] = ecalEndcapLayout.getAbsorberThickness(i);
  }
  absThick[0]=0;
  ECALGarlicGeometryParameters::Instance().Set_absThicknessEndcapLayer(absThick);

  // layer positions: this should be approx position of centre of silicon layer

  float positions[MAX_NUMBER_OF_LAYERS]={0};
  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    positions[i] = ecalBarrelLayout.getDistance(i) + ecalBarrelLayout.getThickness(i)/2;
  }
  positions[0]=positions[1]-ecalBarrelLayout.getThickness(0);
  positions[ecalBarrelLayout.getNLayers()] = positions[ecalBarrelLayout.getNLayers()-1] + ecalBarrelLayout.getThickness(ecalBarrelLayout.getNLayers()-1);
  ECALGarlicGeometryParameters::Instance().Set_positionBarrelLayer(positions);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) positions[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    positions[i] = ecalEndcapLayout.getDistance(i) + ecalEndcapLayout.getThickness(i)/2;
  }
  positions[ecalEndcapLayout.getNLayers()] = positions[ecalEndcapLayout.getNLayers()-1] + ecalEndcapLayout.getThickness(ecalEndcapLayout.getNLayers()-1);

  ECALGarlicGeometryParameters::Instance().Set_positionEndcapLayer(positions);


  return;
}

