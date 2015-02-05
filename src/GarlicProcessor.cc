/*********************************************************************
 * GarlicProcessor
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

v 3-0-pre
---------
Jan 2015, Daniel Jeans
significant redesign: 
- many changes
- no neural networks, simple cuts

optimising for separation of nearby photons

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

#include "GarlicAlgorithmParameters.hh"
#include "GarlicGeometryParameters.hh"

#include "GarlicClusterAlgos.hh"

#include "GarlicExtendedCluster.hh"
#include "GarlicExtendedHit.hh"
#include "GarlicConversionFinder.hh"
#include "GarlicClusterEnergyCorrector.hh"

#include <HelixClass.h>
#include <GarlicProcessor.hh>

using namespace marlin;
using namespace lcio;
using namespace std;

GarlicProcessor anGarlicProcessor;

GarlicProcessor::GarlicProcessor() : Processor("GarlicProcessor") {

  //  streamlog_out ( DEBUG1 ) << "hello from GarlicProcessor constructor " << this << std::endl;

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
                           "TrackCollection",
                           "track collection name",
                           _TrackCollectionName,
                           std::string("MarlinTrkTracks") );

  registerInputCollection( LCIO::TRACK,
                           "TPCTrackCollection",
                           "TPC track collection name",
                           _TPCTrackCollectionName,
                           std::string("ClupatraTracks") );


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
                            std::string("GARLICPhotonClusters") );

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

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "GarlicConversionsPFOColName",
                            "collection of identified conversions",
                            _conversionCollName,
                            std::string("GARLICConversionPFOs") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "GarlicElectronsPFOColName",
                            "collection of identified electrons",
                            _electronPFOCollName,
                            std::string("GARLICElectronPFOs") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "GarlicPhotonsPFOColName",
                            "collection of photon PFOs",
                            _photonPFOCollName,
                            std::string("GARLICPhotonPFOs") );

  // parameters

  registerProcessorParameter("GarlicDebugCollections",
			     "write out debug info collections to lcio file?",
			     _x_debugCollections,
			     false);


  registerProcessorParameter("TrackCheat",
                             "take MC info for tracks?",
                             _x_cheatTracks,
                             false);

  registerProcessorParameter("TrackRemoveNearbyHits",
                             "Should remove Hits near extrapolated tracks to reject pions?",
                             _x_removeHitsNearTracks,
                             true);

  registerProcessorParameter("TrackVetoWindow",
                             "window around track in which to remove hits (mm)",
                             _x_trackWindowVeto,
                             float (10.) );


  // seeding
  registerProcessorParameter("SeedNLayers",
                             "Number of ECAL pseudo layers used for projecting to obtain a seed, typically equivalent to 5 X0.",
                             _x_nLayersForSeeding,
                             int(12));

  registerProcessorParameter("SeedMinHits",
                             "Minimum number of hits to accept a seed.",
                             _x_minHitsForSeeding,
                             int(2));

  registerProcessorParameter("SeedMinHitEnergy",
                             "consider only hits above this threshold (in MIPs) in the seeding",
                             _x_seedHitEnergyCut,
                             float(2.5) );

  registerProcessorParameter("SeedMinEnergy",
                             "reject seeds below this energy (in MIPs)",
                             _x_seedEnergyCut,
                             float (5.) );

  registerProcessorParameter("SeedMaxDistance",
                             "radius of seed cylinder (in terms of cell size)",
                             _x_seedDistanceCut,
                             float (1.5) );

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

  registerProcessorParameter("CoreMaxDistance",
                             "max transverse distance for adding hits to core (in terms of cell size)",
                             _x_maxCoreDist,
                             float (1.5) );


  // electron specific
  registerProcessorParameter("ElectronTransTubeStepSize",
			     "fraction of moliere radius with which to increment electron hit tube at each iteration",
			     _x_ElectronTransTubeStepSize,
			     float (0.75) );

  registerProcessorParameter("ElectronTransTubeNSteps",
			     "number of iterations for electron tube making",
			     _x_ElectronTransNSteps,
			     int (3) );


  // clustering
  registerProcessorParameter("TouchingCellDistance",
			     "definition of distance between cells to be considered touching (multiplier of cell size)",
			     _x_TouchingCellDistance,
			     float ( 2.4 ) );

  registerProcessorParameter("ClusterMaxDist",
                             "Maximum distance from core to added hits (in Moliere Radii)",
                             _x_clusterMaxDist,
                             float(2.));


  // cluster merging


  registerProcessorParameter("MaxMergeDistLowEn",
                             "maximum distance between clusters to consider merging: max(MaxMergeDistLowEn, MaxMergeDistLowEn + MaxMergeDistEnDep*log10(clusEn))",
                             _x_maxMergeDistLowEn,
                             float(2.) );


  registerProcessorParameter("MaxMergeDistEnDep",
                             "maximum distance between clusters to consider merging: max(MaxMergeDistLowEn, MaxMergeDistLowEn + MaxMergeDistEnDep*log10(clusEn))",
                             _x_maxMergeDistEnDep,
                             float(1.5) );


  registerProcessorParameter("stochasticTerm",
                             "assumed stochastic term of energy resolution",
                             _x_stochasticTerm,
                             float(0.17) );

  registerProcessorParameter("constantTerm",
                             "assumed constant term of energy resolution",
                             _x_constantTerm,
                             float(0.01) );

  registerProcessorParameter("MoliereRadius",
                             "assumed Moliere radius of ECAL (mm)",
                             _x_moliereRadius,
                             float( 20. ) );

  registerProcessorParameter("MergeRatioCut",
			     "don't consider merging two clusters if E_min/E_max > MergeRatioCut",
			     _x_MergeRatioCut,
			     float ( 0.25 ) );

  registerProcessorParameter("MergeEnergyDistFactor",
			     "don't consider merging two clusters if E_min/E_max > MergeEnergyDistFactor/(cl-cl dist **2) [cl-cl dist = transverse dist between clusters in mm]",
			     _x_MergeEnergyDistFactor,
			     float(120.) );

  registerProcessorParameter("MergeAbsoluteLargestDist",
			     "don't consider merging two clusters if distance between cluster centres > MergeAbsoluteLargestDist",
			     _x_MergeAbsoluteLargestDist,
			     float (500.) );

  registerProcessorParameter("MergePi0MassLimit",
			     "to veto merging 2 clusters, ask for their invariant mass to be consistent with a value > MergePi0MassLimit @ 2 sigma, and that E_max/E_min < MergePi0MaxEnergyImbalance",
			     _x_MergePi0MassLimit,
			     float(0.1) );

  registerProcessorParameter("MergePi0MaxEnergyImbalance",
			     "to veto merging 2 clusters, ask for their invariant mass to be consistent with a value > MergePi0MassLimit @ 2 sigma, and that E_max/E_min < MergePi0MaxEnergyImbalance",
			     _x_MergePi0MaxEnergyImbalance,
			     float(50.) );


  // photon ID
  registerProcessorParameter("PhotonCutFile",
                             "file containing the photon selection cuts",
                             _x_photonSelFile,
                             std::string("") );

  registerProcessorParameter("EcalEnergyMipConversion",
                             "factor to convert between ECAL hit energy and MIP scale, for first layer",
                             _x_energy_mip_conversion,
                             float (140.) );

  registerProcessorParameter("ForwardTrackAngle",
			     "Tracks with polar angle below this have less stringent cuts for electron finding (rad)",
			     _x_forwardTrackAngle,
			     float (0.15) );


  registerProcessorParameter("DebugMode",
                             "Talk a lot? (0-3)",
                             _x_debug,
                             int(0));

  registerProcessorParameter("clusterCheckHistoFile", "name of file in which to save clustering histograms",
                             _histFileName, std::string(""));

  _clusterer=NULL;
  _fhistos=NULL;
  _track_extrap_col=NULL;
  _mcParticleColl=NULL;
  _cluster_start_col=NULL;

}


GarlicProcessor::~GarlicProcessor() {
  //  streamlog_out ( DEBUG1 ) << "hello from GarlicProcessor destructor " << this << std::endl;

  if (_clusterer) {delete _clusterer; _clusterer=NULL;}
  if (_fhistos) {delete _fhistos; _fhistos=NULL;}

}

void GarlicProcessor::init() {
  streamlog_out ( DEBUG1 ) << "hello from GarlicProcessor init() " << this << std::endl;
  printParameters();
}

void GarlicProcessor::setup()
{

  streamlog_out ( DEBUG1 ) << "hello from GarlicProcessor::setup() " << std::endl;

  printMrGarlic();
  printParameters();


  GarlicAlgorithmParameters::Instance().SetEcalPreClusterCollectionName (_ecalPreClusterCollectionName);

  GarlicAlgorithmParameters::Instance().SetDebug (_x_debug);

  GarlicAlgorithmParameters::Instance().SetTrackCheat     (_x_cheatTracks);
  GarlicAlgorithmParameters::Instance().SetTrackRemoveNearbyHits(_x_removeHitsNearTracks);
  GarlicAlgorithmParameters::Instance().SetTrackVetoWindow (_x_trackWindowVeto);

  GarlicAlgorithmParameters::Instance().SetSeedNLayers (_x_nLayersForSeeding);
  GarlicAlgorithmParameters::Instance().SetSeedMinHits(_x_minHitsForSeeding);
  GarlicAlgorithmParameters::Instance().SetSeedHitEnergyCut(_x_seedHitEnergyCut);
  GarlicAlgorithmParameters::Instance().SetSeedEnergyCut   (_x_seedEnergyCut);
  GarlicAlgorithmParameters::Instance().SetSeedDistanceCut (_x_seedDistanceCut);

  GarlicAlgorithmParameters::Instance().SetClusterMaxDist(_x_clusterMaxDist);

  GarlicAlgorithmParameters::Instance().SetCoreLayersSection1(_x_nlayersSection1);
  GarlicAlgorithmParameters::Instance().SetCoreMaxHoleSection1(_x_maxHoleSection1);
  GarlicAlgorithmParameters::Instance().SetCoreMaxHoleSection2(_x_maxHoleSection2);
  GarlicAlgorithmParameters::Instance().SetCoreDistanceCut(_x_maxCoreDist);

  GarlicAlgorithmParameters::Instance().SetTouchingCellDistance(_x_TouchingCellDistance);

  GarlicAlgorithmParameters::Instance().SetStochasticTerm ( _x_stochasticTerm );
  GarlicAlgorithmParameters::Instance().SetConstantTerm   ( _x_constantTerm   );
  GarlicAlgorithmParameters::Instance().SetMoliereRadius  ( _x_moliereRadius  );

  GarlicAlgorithmParameters::Instance().SetEnergyMIPconversion( _x_energy_mip_conversion );

  GarlicAlgorithmParameters::Instance().SetPhotonCutFile(_x_photonSelFile);

  GarlicAlgorithmParameters::Instance().SetForwardTrackAngle(_x_forwardTrackAngle);


  GarlicAlgorithmParameters::Instance().SetMergeMaxDistAtLowEn( _x_maxMergeDistLowEn );
  GarlicAlgorithmParameters::Instance().SetMergeMaxDistEnDep( _x_maxMergeDistEnDep );

  GarlicAlgorithmParameters::Instance().SetMergeRatioCut              (_x_MergeRatioCut             );
  GarlicAlgorithmParameters::Instance().SetMergeEnergyDistFactor      (_x_MergeEnergyDistFactor     );
  GarlicAlgorithmParameters::Instance().SetMergeAbsoluteLargestDist   (_x_MergeAbsoluteLargestDist  );
  GarlicAlgorithmParameters::Instance().SetMergePi0MassLimit          (_x_MergePi0MassLimit         );
  GarlicAlgorithmParameters::Instance().SetMergePi0MaxEnergyImbalance (_x_MergePi0MaxEnergyImbalance);


  GarlicAlgorithmParameters::Instance().SetElectronTransTubeStepSize ( _x_ElectronTransTubeStepSize );
  GarlicAlgorithmParameters::Instance().SetElectronTransNSteps ( _x_ElectronTransNSteps );


  _clusterer = new GarlicClusterAlgos();

  if ( GarlicAlgorithmParameters::Instance().GetDebug () > 2 )
    _clusterer->setVerbose(true);

  _energyCorrector = new GarlicClusterEnergyCorrector();


  _nEvents=0;
  _nPhotonClusters=0;

  if (_histFileName!="")
    _fhistos = new TFile(_histFileName.c_str(),"recreate");

  _convFinder=NULL;

  _nSaveHist=0;
  _geomSetup = false;

  return;
}



void GarlicProcessor::processRunHeader(LCRunHeader * run)
{
  return;
}


void GarlicProcessor::processEvent(LCEvent * evt)   // main !
{
  streamlog_out ( DEBUG1 ) << endl
                           << endl << "Event: " << evt->getEventNumber() << endl;

  // set up some utilities for the first event:
  // ---- the detector geometry
  if (!_geomSetup) {
    setUpGeometry();
    _geomSetup=true;
  }
  // --- the conversion finder
  if ( ! _convFinder ) {
    _convFinder = new GarlicConversionFinder(GarlicGeometryParameters::Instance().Get_bField(), true);
  }


  _nEvents++;

  // define the input and output collections

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

  LCCollectionVec* finalClusterColl=0;
  LCCollectionVec* finalTightSelectedClusterColl=0;
  LCCollectionVec* finalLooseSelectedClusterColl=0;
  LCCollectionVec* finalVeryLooseSelectedClusterColl=0;
  LCCollectionVec* finalRejectedClusterColl=0;

  LCCollectionVec* finalClusterCollParameters=0;
  LCCollectionVec* finalCluParLinkColl=0;

  LCCollectionVec* electronTightClusterColl=0;
  LCCollectionVec* electronLooseClusterColl=0;
  LCCollectionVec* electronVeryLooseClusterColl=0;

  LCCollectionVec* tightConversionColl=0;
  LCCollectionVec* looseConversionColl=0;
  LCCollectionVec* trackClusterColl=0;

  LCCollectionVec* MCphotonClusterColl=0;
  LCCollectionVec* MCelectronClusterColl=0;
  LCCollectionVec* MCchHadClusterColl=0;
  LCCollectionVec* MCneuHadClusterColl=0;
  LCCollectionVec* MCotherClusterColl=0;

  LCCollectionVec* wronglySelectedClusterColl=0;
  LCCollectionVec* wronglyRejectedClusterColl=0;

  LCCollectionVec* electronTightPFOColl=0;
  LCCollectionVec* electronLoosePFOColl=0;
  LCCollectionVec* electronVeryLoosePFOColl=0;

  LCCollectionVec* photonTightPFOColl=0;
  LCCollectionVec* photonLoosePFOColl=0;
  LCCollectionVec* photonVeryLoosePFOColl=0;

  IMPL::LCCollectionVec *seedCoreLink=0;
  IMPL::LCCollectionVec *seedFinalClusLink=0;

  IMPL::LCCollectionVec *seed_col = 0;
  IMPL::LCCollectionVec *rej_seed_col = 0;
  _track_extrap_col = 0;

  finalClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalClusterColl->setFlag( finalClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  finalTightSelectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalTightSelectedClusterColl->setFlag( finalTightSelectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalTightSelectedClusterColl->setSubset();

  finalLooseSelectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalLooseSelectedClusterColl->setFlag( finalLooseSelectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalLooseSelectedClusterColl->setSubset();

  finalVeryLooseSelectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalVeryLooseSelectedClusterColl->setFlag( finalVeryLooseSelectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalVeryLooseSelectedClusterColl->setSubset();

  finalRejectedClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  finalRejectedClusterColl->setFlag( finalRejectedClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  finalRejectedClusterColl->setSubset();

  photonTightPFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  photonLoosePFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  photonVeryLoosePFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  electronTightClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  electronTightClusterColl->setFlag( electronTightClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  electronLooseClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  electronLooseClusterColl->setFlag( electronLooseClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  electronVeryLooseClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  electronVeryLooseClusterColl->setFlag( electronVeryLooseClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

  electronTightPFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  electronLoosePFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  electronVeryLoosePFOColl = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  looseConversionColl     = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  tightConversionColl     = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  if ( _x_debugCollections ) { // these are collections of objects potentially useful for (expert) debugging

    finalClusterCollParameters = new LCCollectionVec(LCIO::LCGENERICOBJECT);

    finalCluParLinkColl = new LCCollectionVec(LCIO::LCRELATION);

    coreClusterColl = new LCCollectionVec(LCIO::CLUSTER);
    coreClusterColl->setFlag( coreClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

    trackClusterColl = new LCCollectionVec(LCIO::CLUSTER);
    trackClusterColl->setFlag( trackClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );

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

    seed_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    seed_col->setFlag(seed_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

    rej_seed_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    rej_seed_col->setFlag(rej_seed_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

    _track_extrap_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    _track_extrap_col->setFlag(_track_extrap_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

    _cluster_start_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    _cluster_start_col->setFlag(_cluster_start_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

    seedFinalClusLink = new LCCollectionVec(LCIO::LCRELATION);
    seedCoreLink = new LCCollectionVec(LCIO::LCRELATION);

  }




  //-------------------------------------------
  // prepare the Track collection
  //-------------------------------------------
  vector <GarlicExtendedTrack* > trackVec;
  vector <GarlicExtendedCluster* > preClusVec;

  if ( GarlicAlgorithmParameters::Instance().GetTrackCheat() == 0 )
    PrepareTracks(evt,trackVec);
  else
    PrepareMCTracks(trackVec);


  //-------------------------------------------
  // 1.) transform PreClusters to ExtendedClusters (pseudoIDs etc)
  //-------------------------------------------
  streamlog_out ( DEBUG2 ) << "Preparing PreClusters... " << endl;
  PreparePreClusters(evt, preClusVec);
  streamlog_out ( DEBUG2 ) << preClusVec.size() << " PreClusters prepared" << endl;
  if(preClusVec.size()==0) {
    streamlog_out ( MESSAGE ) << "no preclusters found!" << endl;
    return;
  }

  //-------------------------------------------
  // find the tracks near each precluster
  //-------------------------------------------
  std::map < GarlicExtendedCluster*, vector <GarlicExtendedTrack* > > nearbyTracks;
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {
    nearbyTracks[ preClusVec[roi_i] ] = selectNearbyTracks( preClusVec[roi_i], &trackVec);
  }
  streamlog_out ( DEBUG2 ) << "associated tracks to preclusters...now look for electrons" << endl;


  //-------------------------------------------
  // look for electrons (seeded by tracks)
  //-------------------------------------------
  vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > allElectrons;
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {
    streamlog_out ( DEBUG2 )
      << "looking for electrons in roi " << roi_i << endl
      << " -- preclus energy, ntrks: " << preClusVec[roi_i]->getEnergy() << " " << nearbyTracks[ preClusVec[roi_i] ].size() << endl;
    vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > electrons = _clusterer->getElectrons(preClusVec[roi_i], nearbyTracks[ preClusVec[roi_i] ]);
    streamlog_out ( DEBUG2 ) << "  found " << electrons.size() << "electrons" << endl;
    if (electrons.size()>0) {
      for (size_t i=0; i<electrons.size(); i++) {
        streamlog_out ( DEBUG2 ) << i << " " << electrons[i].first->getTotalMomentum() << " " << electrons[i].second->getEnergy() << " " << electrons[i].first << endl;
        allElectrons.push_back(electrons[i]);
      }
    }
    RemoveElectronHits( electrons, preClusVec[roi_i] );

  }
  streamlog_out ( DEBUG2 ) << "got electrons, now look for conversions" << endl;

  //-------------------------------------------
  // look for conversions using all tracks
  //-------------------------------------------
  _convFinder->setTracks( trackVec );
  std::vector < GarlicConversionInfo > identifiedConversions = _convFinder->getConversions();
  streamlog_out ( DEBUG2 ) << " found " << identifiedConversions.size() << " conversion candidates" << endl;
  streamlog_out ( DEBUG2 ) << "got electrons and conversions: force all conv tracks as electrons" << endl;

  // for conversion tracks not ID'd as electrons, force them to become electrons
  for (size_t j=0; j<identifiedConversions.size(); j++) {
    for (int k=0; k<2; k++) {
      if ( identifiedConversions[j].electronID[k]==0 ) {

        streamlog_out ( DEBUG2 ) << "  conv. " << j << " trk " << k << endl;

        GarlicExtendedTrack* etrk = identifiedConversions[j].etrks[k]; // this etrk from a conversion has not been identified as an electron

        // find the precluster near this track
        GarlicExtendedCluster* thisPreClus(0);
        for ( std::map < GarlicExtendedCluster*, vector <GarlicExtendedTrack* > >::iterator itt=nearbyTracks.begin(); itt!=nearbyTracks.end(); itt++) {
          vector <GarlicExtendedTrack* > tt = itt->second;
          for (size_t l=0; l<tt.size(); l++) {
            if ( tt[l]==etrk ) {
              thisPreClus=itt->first; // got it!
              break;
            }
          }
        }

        streamlog_out ( DEBUG2 ) << "   trk & cluster: " << etrk << " " << thisPreClus << endl;

        if ( etrk && thisPreClus ) {
          std::vector < GarlicExtendedTrack* > ttrk;
          ttrk.push_back( etrk );
          bool forceElectron=true;
          std::vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > electrons = _clusterer->getElectrons( thisPreClus, ttrk, forceElectron);

          streamlog_out ( DEBUG2 ) << "forced electron from conversion track! " << electrons.size() << endl;

          if (electrons.size()>0) {
            for (size_t i=0; i<electrons.size(); i++) {
              streamlog_out ( DEBUG2 ) << i << " " << electrons[i].first->getTotalMomentum() << " " << electrons[i].second->getEnergy() << endl;
              allElectrons.push_back(electrons[i]);
            }
          }
          RemoveElectronHits( electrons, thisPreClus );
        }

      }
    }
  }

  //-------------------------------------------
  // 1b.) clear hits near to extrapolated tracks
  //-------------------------------------------
  streamlog_out ( DEBUG2 ) << "now veto hits near tracks" << endl;

  std::map < GarlicExtendedTrack*, GarlicExtendedCluster*> track_clusters;
  if(GarlicAlgorithmParameters::Instance().GetTrackRemoveNearbyHits()) {
    streamlog_out ( DEBUG2 ) << "Removing hits near tracks..." << endl;
    track_clusters = RemoveHitsNearExtrapolatedTracks(trackVec,preClusVec);
    if ( trackClusterColl ) {
      for ( std::map < GarlicExtendedTrack*, GarlicExtendedCluster*>::iterator itt=track_clusters.begin(); itt!=track_clusters.end(); itt++) {
	ClusterImpl* climp = new ClusterImpl(itt->second->getClusterImpl());
	trackClusterColl->addElement( climp );
	// cout << "does this track cluster look like MIP? " << itt->second->getIsMipLike() << endl;
      }
    }
  }

  streamlog_out ( DEBUG2 ) << "now loop over preclusters, look for photon clusters" << endl;

  //-------------------------------------------
  // the next steps are done per ROI (preCluster)
  //-------------------------------------------
  for (size_t roi_i=0; roi_i<preClusVec.size(); roi_i++) {

    GarlicExtendedCluster *preCluster=preClusVec[roi_i];

    streamlog_out ( DEBUG2 )
      << "evt" << evt->getEventNumber() << " roi " << roi_i << " nhits=" << preCluster->getHits()->size() << endl
      << "Building histogram for seeding for PreCluster " << roi_i << endl;

    if (preCluster->getHits()->size()==0) {
      streamlog_out ( DEBUG2 ) << "precluster/RoI with no hits, ignoring!" << endl;
      continue;
    }

    // get tracks close to this precluster
    //    vector <GarlicExtendedTrack* > nearbyTracks = selectNearbyTracks(preCluster, &trackVec);

    //-------------------------------------------
    // get the seeds
    //-------------------------------------------
    if (_fhistos && _nSaveHist++<MAXSAVEHIST) {
      _fhistos->cd();
      TString hn = "evt"; hn+=evt->getEventNumber();
      hn+="_roi"; hn+=roi_i;
      _clusterer->saveHistos(_fhistos, hn);
    } else _clusterer->saveHistos(NULL);

    std::map <CalorimeterHit*, bool> allseeds = _clusterer->getSeeds(preCluster);
    // add seeds to collection
    std::vector <CalorimeterHit*> goodseeds;
    for ( std::map <CalorimeterHit*, bool>::iterator ijs=allseeds.begin(); ijs!=allseeds.end(); ijs++ ) {
      if (ijs->second) {
        goodseeds.push_back(ijs->first);
        if ( seed_col ) seed_col->addElement(ijs->first);
      } else {
        if ( rej_seed_col ) rej_seed_col->addElement(ijs->first);
      }
    }

    if (goodseeds.size()==0) continue;

    streamlog_out ( DEBUG2 ) << "GOOD SEEDS = " << goodseeds.size() << endl;


    //-------------------------------------------
    // build the cores
    //-------------------------------------------
    std::map < CalorimeterHit*, GarlicExtendedCluster* > cores =  _clusterer->getCores(preCluster, goodseeds);
    streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " ncores = " << cores.size() << endl;

    if ( coreClusterColl ) {
      for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator itt=cores.begin(); itt!=cores.end(); itt++) {
	ClusterImpl* climp = new ClusterImpl(itt->second->getClusterImpl());
	coreClusterColl->addElement(climp);
	streamlog_out ( DEBUG2 ) << "   saving core with nhits = " << climp->getCalorimeterHits().size() << endl;

	if ( seedCoreLink ) {
	  LCRelationImpl* rel = new LCRelationImpl(itt->first, climp);
	  seedCoreLink->addElement(rel);
	}
      }
    }

    //-------------------------------------------
    // build the clusters
    //-------------------------------------------
    std::map < CalorimeterHit*, GarlicExtendedCluster* > Cclusters = _clusterer->getClusters(preCluster, cores);
    streamlog_out ( DEBUG2 ) << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " nclusters = " << Cclusters.size() << endl;

    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator gg=Cclusters.begin(); gg!=Cclusters.end(); gg++) {
      streamlog_out ( DEBUG2 ) << "cluster: " << gg->second->getEnergy() << endl;
    }


    //-------------------------------------------
    // merge photon clusters and electron clusters with each other
    //-------------------------------------------
    _clusterer->mergeAllSatellitesAndElectrons(Cclusters, allElectrons);

    RemoveElectronHits( allElectrons, preClusVec[roi_i] ); // after merging...not sure if we need to do this...

    streamlog_out ( DEBUG2 )
      << "evt" << evt->getEventNumber() << " --- roi " << roi_i << " nclusters (after merge) = " << Cclusters.size() << endl;

    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator gg=Cclusters.begin(); gg!=Cclusters.end(); gg++) {
      streamlog_out ( DEBUG2 ) << "cluster: " << gg->second->getEnergy() << endl;
    }

    //-------------------------------------------
    // look at remaining unclustered hits
    //-------------------------------------------
    std::vector < GarlicExtendedHit* > unclusteredHits;
    for (size_t ihh=0; ihh<preCluster->getHits()->size(); ihh++) {
      GarlicExtendedHit* ehit = (*(preCluster->getHits()))[ihh];
      bool foundit=false;

      // first check the photon clusters
      for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator icc=Cclusters.begin(); icc!=Cclusters.end(); icc++) {
        std::vector <GarlicExtendedHit*> clhits = *(icc->second->getHits());
        if ( find ( clhits.begin(), clhits.end(), ehit ) != clhits.end() ) {
          foundit=true;
          break;
        }
      }
      if (foundit) {
        //      cout << "WEIRD, found this hit in a photon cluster..." << foundit << endl;
        continue;
      }

      // and the electrons
      // vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > allElectrons;
      for (size_t iel=0; iel<allElectrons.size(); iel++) {
        std::vector <GarlicExtendedHit*> clhits = *(allElectrons[iel].second->getHits());
        if ( find ( clhits.begin(), clhits.end(), ehit ) != clhits.end() ) {
          foundit=true;
          break;
        }
      }
      if (foundit) {
        streamlog_out( WARNING ) << "found this hit in an electron... " << foundit << endl;
        continue;
      }

      unclusteredHits.push_back(ehit);
    }

    //-------------------------------------------
    // add tracks to clusters
    //-------------------------------------------
    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator itt=Cclusters.begin(); itt!=Cclusters.end(); itt++) {
      itt->second->setTracks( & nearbyTracks[ preClusVec[roi_i] ] );
    }

    //    cout << "final number of photon clusters (before final selection): " << Cclusters.size() << endl;

    //-------------------------------------------
    // apply cuts to identified photon clusters
    //  these are cuts which typically depend on environment, not only on the cluster itself
    //   this is to reject hadron fragments
    //-------------------------------------------

    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator itt=Cclusters.begin(); itt!=Cclusters.end(); itt++) {
      ClusterImpl* climp = new ClusterImpl(itt->second->getClusterImpl());
      finalClusterColl->addElement(climp);

      streamlog_out ( DEBUG2 ) << "applying basic cuts..." << itt->second->getEnergy() << endl;

      int photonsel = itt->second->getPhotonCutSel();

      bool cutSelected_T  = photonsel>=GarlicClusterSelector::TIGHT;
      bool cutSelected_L  = photonsel>=GarlicClusterSelector::LOOSE;
      bool cutSelected_VL = photonsel>=GarlicClusterSelector::VERYLOOSE;

      streamlog_out ( DEBUG2 ) << "photon selection, with pointing requirement: " << itt->second->getPhotonCutSel() << endl;

      streamlog_out ( DEBUG2 )  << "done applying cuts...VL, L, T = " << cutSelected_VL << " " << cutSelected_L << " " << cutSelected_T << endl
                                << " p-layer frac = " << itt->second->getFracPseudoLayers() << endl
                                << " shower length (X0) = " << itt->second->getEnd() - itt->second->getStart() << endl
                                << " shower start (X0) = " << itt->second->getStart() << endl;


      // --------------------------------------------------
      // distance of cluster from nearest charged track
      // -------------------------------------------------

      // is the nearest track an electron?
      GarlicExtendedTrack* clTrk = itt->second->getClosestTrack();
      bool closestIsElectron(false);
      for ( size_t kk=0; kk<allElectrons.size(); kk++) {
        if ( allElectrons[kk].first == clTrk ) {
          closestIsElectron=true;
          break;
        }
      }

      streamlog_out ( DEBUG2 ) << "is the closest track identified as an electron? " << closestIsElectron << " ( " << clTrk << " " << endl;
      float minTrackDist = 10; // hard coded!!!
      float dd[4];
      dd[0] = itt->second->getTrackDist_cog();  //
      dd[1] = itt->second->getTrackDist_proj(); // different measures of distance to track
      dd[2] = itt->second->getTrackDist_min();  //
      dd[3] = itt->second->getTrackDist_first();//

      streamlog_out ( DEBUG2 ) << "distances " << dd[0] << " " << dd[1] << " " << dd[2] << " " << dd[3] << " , " << minTrackDist << endl;

      int nclose(0);
      int nclose3(0);
      int nclose10(0);

      for (int i=0; i<4; i++) {
        if (dd[i]<minTrackDist) nclose++;
        if (dd[i]<3*minTrackDist) nclose3++;
        if (dd[i]<10*minTrackDist) nclose10++;
      }

      // ----------------------------------------------------------
      // check how many track removed hits are within 1 and 2 moliere radii
      // ----------------------------------------------------------
      int nRemHits_1rm(0);
      int nRemHits_2rm(0);
      float enRemHits_1rm(0);
      float enRemHits_2rm(0);

      for ( std::map < GarlicExtendedTrack*, GarlicExtendedCluster*>::iterator itc = track_clusters.begin(); itc!=track_clusters.end(); itc++) {
        GarlicExtendedCluster* ecl = itc->second;
        for ( size_t j=0; j<ecl->getHits()->size(); j++) {
          CalorimeterHit* hit = ecl->getHits()->at(j)->getCaloHit();
          float dist = min( itt->second->getDistToClusterAxis( hit->getPosition(), 0 ), itt->second->getDistToClusterAxis( hit->getPosition(), 1 ) );
          if ( dist < GarlicAlgorithmParameters::Instance().GetMoliereRadius() ) {
            nRemHits_1rm++;
            enRemHits_1rm+=hit->getEnergy();
          }
          if ( dist < 2*GarlicAlgorithmParameters::Instance().GetMoliereRadius() ) {
            nRemHits_2rm++;
            enRemHits_2rm+=hit->getEnergy();
          }
        }
      }

      streamlog_out ( DEBUG2 ) << "removed hits, energy within 1, 2 rM = " <<
        nRemHits_1rm << " " << enRemHits_1rm << " , " <<
        nRemHits_2rm << " " << enRemHits_2rm << endl;

      // ----------------------------------------------------------
      // is there a bigger photon within a couple of RM?
      // ----------------------------------------------------------
      bool nearbyLargeCluster(false);
      float* thipos = itt->second->getCentreOfGravity();
      for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator jtt=Cclusters.begin(); jtt!=Cclusters.end(); jtt++) {
        if (jtt==itt) continue;
        if ( jtt->second->getEnergy() < 5.*itt->second->getEnergy() ) continue; // must be larger in energy
        float* clpos = jtt->second->getCentreOfGravity();
        float dist(0);
        for (int i=0; i<3; i++) dist+=pow( thipos[i]-clpos[i],2);
        dist=sqrt(dist);
        if ( dist < GarlicAlgorithmParameters::Instance().GetMoliereRadius() ) {
          nearbyLargeCluster=true;
          break;
        }
      }

      streamlog_out ( DEBUG2 ) << " has nearby large cluster? " << nearbyLargeCluster << endl;

      // does this cluster point to a nearby identified interaction position on a track cluster?
      float mindistToInt(99999);
      for (std::map < GarlicExtendedTrack*, GarlicExtendedCluster*>::iterator jtt=track_clusters.begin(); jtt!=track_clusters.end(); jtt++) {
        mindistToInt = min(mindistToInt, itt->second->getDistToClusterAxis( jtt->second->getInteractionPosition(), 1 ) );
      }
      streamlog_out ( DEBUG2 ) << "cluster distance to identified calo-track interaction " << mindistToInt << endl;

      // in the case of low energy or late starting showers, apply some more stringent cuts
      bool isLateShower = itt->second->getStart() > 5; // x0
      bool isLowEn      = itt->second->getEnergy() < 0.8; // GeV

      bool mipLike =  itt->second->getIsMipLike();


      streamlog_out ( DEBUG2 ) << " late shower, lowEn, miplike ? " << isLateShower << " " << isLowEn << " " << mipLike << endl;


      // ------------------------------------------------
      // this is a bunch of hard-coded cuts whose main role is to remove hadronic fragments
      // this code should be better organised in future
      // ------------------------------------------------

      std::vector < int > rejectReasons;
      if      ( nclose>=2 && !closestIsElectron )
        rejectReasons.push_back(1);

      if ( isLowEn && nclose3>=2 && !closestIsElectron)
        rejectReasons.push_back(2);

      if ( itt->second->getNLayers() < 3 || itt->second->getNPseudoLayers() < 3 )
        rejectReasons.push_back(3);

      if ( ( !isLowEn && itt->second->getFracPseudoLayers() < 0.75 ) || ( isLowEn && itt->second->getFracPseudoLayers() < 0.50 )  )
        rejectReasons.push_back(4);

      if ( enRemHits_1rm > 3*itt->second->getEnergy() )
        rejectReasons.push_back(5);

      if ( (isLateShower && enRemHits_2rm > 2*itt->second->getEnergy()) ||
           (isLowEn      && enRemHits_2rm > 2*itt->second->getEnergy()) ||
           (                enRemHits_2rm > 6*itt->second->getEnergy()) )
        rejectReasons.push_back(6);

      if ( isLateShower && nclose10>=2 )
        rejectReasons.push_back(7);

      if ( isLateShower && nearbyLargeCluster )
        rejectReasons.push_back(8);

      if ( mipLike )
        rejectReasons.push_back(9);

      if ( mindistToInt < GarlicAlgorithmParameters::Instance().GetMoliereRadius() )
        rejectReasons.push_back(10);

      if ( rejectReasons.size()>0 ) {
        streamlog_out ( DEBUG2 ) << "rejection reasons: ";
        for (size_t i=0; i<rejectReasons.size(); i++)
          streamlog_out ( DEBUG2 ) << rejectReasons[i] << " ";
        streamlog_out ( DEBUG2 ) << endl;
      }


      int selected(0);
      if      ( cutSelected_T && rejectReasons.size()==0 ) selected=GarlicClusterSelector::TIGHT;
      else if ( cutSelected_L && rejectReasons.size()==0 ) selected=GarlicClusterSelector::LOOSE;
      else if ( cutSelected_VL && rejectReasons.size()==0 ) selected=GarlicClusterSelector::VERYLOOSE;


      // -------------------------------------------------
      // now make the output PFO objects for selected photons
      // ------------------------------------------------

      if ( selected>0 ) {

        // cout << "selected photon cluster with energy " << climp->getEnergy() << endl;

        ReconstructedParticleImpl* gammaRP = new ReconstructedParticleImpl();
        gammaRP->addCluster( climp );
        float mom[3];
        float en = _energyCorrector->getCorrectedEnergy(itt->second);
        float tt(0);
        for (int i=0; i<3; i++) {
          mom[i]=itt->second->getCentreOfGravity()[i];
          tt+=pow(mom[i],2);
        }
        tt=sqrt(tt);
        for (int i=0; i<3; i++) {
          mom[i]*=en/tt;
        }
        gammaRP->setMomentum( mom );
        gammaRP->setEnergy( en );

        if (selected==GarlicClusterSelector::TIGHT) {
          streamlog_out ( DEBUG2 ) << "tight selected cluster! " << climp->getEnergy() << endl;
          finalTightSelectedClusterColl->addElement(climp);
          photonTightPFOColl->addElement( gammaRP );
        } else if (selected==GarlicClusterSelector::LOOSE) {
          streamlog_out ( DEBUG2 )  << "loose selected cluster! " << climp->getEnergy() << endl;
          finalLooseSelectedClusterColl->addElement(climp);
          photonLoosePFOColl->addElement( gammaRP );
        } else if (selected==GarlicClusterSelector::VERYLOOSE) {
          streamlog_out ( DEBUG2 )  << "very loose selected cluster! " << climp->getEnergy() << endl;
          finalVeryLooseSelectedClusterColl->addElement(climp);
          photonVeryLoosePFOColl->addElement( gammaRP );
        }
      } else { // just store the cluster
        streamlog_out ( DEBUG2 ) << "rejected cluster! " << climp->getEnergy() << endl;
        finalRejectedClusterColl->addElement(climp);
      }

      CalorimeterHitImpl *clstart=new CalorimeterHitImpl();
      clstart->setEnergy(500);
      clstart->setPosition( itt->second->getClusterStartPosition() );
      if ( _cluster_start_col ) _cluster_start_col->addElement(clstart);


      // -------------------------------------
      // links between objects
      // -------------------------------------

      if (finalClusterCollParameters && finalCluParLinkColl ) {
	IMPL::LCGenericObjectImpl* genobj = new IMPL::LCGenericObjectImpl(itt->second->getGenericObject());
	finalClusterCollParameters->addElement(genobj);

	LCRelationImpl* rel = new LCRelationImpl(climp, genobj);
	finalCluParLinkColl->addElement(rel);
      }

      if ( seedFinalClusLink ) {
	LCRelationImpl* rel = new LCRelationImpl(itt->first, climp);
	seedFinalClusLink->addElement(rel);
      }

      itt->second->setHitRelationColl(_caloHitRelationColl);
      int pdg = itt->second->getMCCalPdg()[0]; // dominant pdg
      //        float frac = itt->second->getMCCalFrac()[0]; // the fraction of energy from this pdg


      TString PDG;

      switch (pdg) {
      case (22):
        if ( MCphotonClusterColl ) MCphotonClusterColl->addElement(climp);
        PDG="photon";
        break;
      case (11):
        if ( MCelectronClusterColl ) MCelectronClusterColl->addElement(climp);
        PDG="electron";
        break;
      case (211):
      case (321):
      case (2212):
        PDG="chHad";
        if ( MCchHadClusterColl ) MCchHadClusterColl->addElement(climp);
        break;
      case (130):
      case (310):
      case (2112):
        PDG="neuHad";
        if ( MCneuHadClusterColl ) MCneuHadClusterColl->addElement(climp);
        break;
      default:
        PDG="other";
        if ( MCotherClusterColl ) MCotherClusterColl->addElement(climp);
        if (pdg!=-999)
          streamlog_out( DEBUG ) << "undefined pdg " << pdg << endl;
        break;
      }

      if ( wronglySelectedClusterColl ) {
	if (selected && pdg!=22) {
	  wronglySelectedClusterColl->addElement(climp);
	  streamlog_out ( DEBUG2 ) << "making a mistake! selecting a " << pdg << " " << climp->getEnergy() << " event " << evt->getEventNumber() << endl;
	}
      }

      if ( wronglyRejectedClusterColl ) {
	if (!selected && pdg==22) {
	  wronglyRejectedClusterColl->addElement(climp);
	  streamlog_out ( DEBUG2 )  << "making a mistake! rejecting a " << pdg << " " << climp->getEnergy() << " event " << evt->getEventNumber() << endl;
	}
      }

    }


    // the seed has been added to a collection, as has the cluster.
    // however, we need to delete the extendedcluster
    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator itt=Cclusters.begin(); itt!=Cclusters.end(); itt++) {
      delete itt->second; itt->second=NULL;
    }

  } // RoI loop


  // -------------------------------------
  // make the electron clusters
  // -------------------------------------
  std::vector < ClusterImpl* > elecClusters;
  if (allElectrons.size()>0) {
    // add to event
    for (size_t ic=0; ic<allElectrons.size(); ic++) {
      ClusterImpl* climp = new ClusterImpl(allElectrons[ic].second->getClusterImpl());

      int sel = allElectrons[ic].second->getElectronCutSel();

      if      ( sel==GarlicClusterSelector::TIGHT )     electronTightClusterColl->addElement(climp);
      else if ( sel==GarlicClusterSelector::LOOSE )     electronLooseClusterColl->addElement(climp);
      else if ( sel==GarlicClusterSelector::VERYLOOSE ) electronVeryLooseClusterColl->addElement(climp);
      
      //        electronClusterColl->addElement(climp);
      elecClusters.push_back(climp); // for later adding to PFOs
    }
  }



  // -------------------------------------
  // make the conversion PFOs
  // -------------------------------------
  streamlog_out ( DEBUG2 ) << "making conversion PFOs! " << identifiedConversions.size() << endl;
  for (size_t jj=0; jj<identifiedConversions.size(); jj++) {

    ReconstructedParticleImpl* convRP = new ReconstructedParticleImpl();
    convRP->addTrack( identifiedConversions[jj].trks[0] );
    convRP->addTrack( identifiedConversions[jj].trks[1] );
    // check if we have an electron cluster for either of these tracks
    // if so, add it to the conversion PFO
    for (size_t ic=0; ic<allElectrons.size(); ic++) {
      if ( allElectrons[ic].first == identifiedConversions[jj].etrks[0] ||
           allElectrons[ic].first == identifiedConversions[jj].etrks[1] ) {
        convRP->addCluster( elecClusters[ic] );
      }
    }
    float mom[3];
    float en(0);
    for (int i=0; i<3; i++) {
      mom[i]=identifiedConversions[jj].momentum[i];
      en+=pow(mom[i],2);
    }

    en=sqrt(en);
    convRP->setMomentum( mom );
    convRP->setEnergy( en );

    // check how many tracks have electron ID
    // also if they are curlers: the tracks can break and cause problems in the extrapoaltion to the CALOs

    int possibleElectrons(0);
    for (int k=0; k<2; k++) {
      streamlog_out( DEBUG ) << "conversion track " << identifiedConversions[jj].etrks[k]->getTotalMomentum() << endl;
      streamlog_out( DEBUG ) << "  electron id = " << identifiedConversions[jj].electronID[k] << endl;
      streamlog_out( DEBUG ) << "ecal entry position: ";
      for (int kk=0; kk<3; kk++)
        streamlog_out( DEBUG ) << identifiedConversions[jj].etrks[k]->getEcalEntryPos()[kk] << " " ;
      streamlog_out( DEBUG ) << endl;
      // is this track a curler?
      float nCurls(0);
      if ( fabs( identifiedConversions[jj].etrks[k]->getEcalEntryPos()[2] ) > GarlicGeometryParameters::Instance().Get_zOfBarrel() ) {
        float endcap_s = GarlicGeometryParameters::Instance().Get_zOfBarrel() /
          identifiedConversions[jj].trks[k]->getTanLambda();
        nCurls = fabs( endcap_s / ( 2.*acos(-1)/identifiedConversions[jj].trks[k]->getOmega() ) );
        streamlog_out( DEBUG ) << " ; ncurls = " << nCurls << endl;
      }

      if ( track_clusters.find( identifiedConversions[jj].etrks[k] )==track_clusters.end() ) {
        streamlog_out( DEBUG ) << "this track has NO track cluster" << endl;
      }

      if ( identifiedConversions[jj].electronID[k] || nCurls > 0.5 ) possibleElectrons++;

    }

    // require at least one electron ID, or 2 curlers
    if ( possibleElectrons>0 ) {
      if      ( identifiedConversions[jj].selTight )
        tightConversionColl->addElement( convRP );
      else if ( identifiedConversions[jj].selLoose )
        looseConversionColl->addElement( convRP );
    } else {
      delete convRP; convRP=NULL;
    }

  }


  // -------------------------------------
  // make the electron PFOs
  // -------------------------------------
  streamlog_out ( DEBUG2 ) << "making electron collection: " << allElectrons.size() << endl;
  // make collection of electron objects
  //  vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > allElectrons;
  if (allElectrons.size()>0) {
    // add to event
    for (size_t ic=0; ic<allElectrons.size(); ic++) {
      // was this electron identified as a conversion?
      bool isConv=false;
      for (size_t jj=0; jj<identifiedConversions.size(); jj++) {
	if ( allElectrons[ic].first == identifiedConversions[jj].etrks[0] ||
	     allElectrons[ic].first == identifiedConversions[jj].etrks[1] ) {
	  isConv=true;
	}
      }
      if ( isConv ) 
	continue; // already added cluster to converion PFO
      else {
	ReconstructedParticleImpl* eleRP = new ReconstructedParticleImpl();
	eleRP->addTrack( allElectrons[ic].first->getTrack() );
	eleRP->addCluster( elecClusters[ic] );
	float mom[3];
	float en(0);
	for (int i=0; i<3; i++) {
	  mom[i]=allElectrons[ic].first->getMomentum()[i];
	  en+=pow(mom[i],2);
	}
	en=sqrt(en);
	eleRP->setMomentum( mom );
	eleRP->setEnergy( en );
	int sel = allElectrons[ic].second->getElectronCutSel();

	streamlog_out( DEBUG ) << "adding electron PFO, with selection " << sel << endl;

	if      ( sel==GarlicClusterSelector::TIGHT )     electronTightPFOColl->addElement(eleRP);
	else if ( sel==GarlicClusterSelector::LOOSE )     electronLoosePFOColl->addElement(eleRP);
	else if ( sel==GarlicClusterSelector::VERYLOOSE ) electronVeryLoosePFOColl->addElement(eleRP);
      }
    }

    streamlog_out ( DEBUG2 ) << "done electrons" << endl;


    for (size_t i=0; i<allElectrons.size(); i++) {
      if (allElectrons[i].second) {
        streamlog_out ( DEBUG2 )  << "deleting electron extended cluster! " << allElectrons[i].second << endl;
        delete allElectrons[i].second;
        allElectrons[i].second=NULL;
      }
    }
  }

  addCollToEvent(seed_col, _seedCollName, evt);
  addCollToEvent(rej_seed_col, _seedCollName+"Rejected", evt);

  addCollToEvent(_track_extrap_col, _trkExtrapCollName, evt);

  addCollToEvent(coreClusterColl, _coreCollName, evt);

  addCollToEvent(finalClusterColl, _clusterCollName+"All", evt);
  addCollToEvent(finalTightSelectedClusterColl, _clusterCollName+"TightSel", evt);
  addCollToEvent(finalLooseSelectedClusterColl, _clusterCollName+"LooseSel", evt);
  addCollToEvent(finalVeryLooseSelectedClusterColl, _clusterCollName+"VeryLooseSel", evt);
  addCollToEvent(finalRejectedClusterColl, _clusterCollName+"Rejected", evt);

  addCollToEvent(finalClusterCollParameters, _clparsCollName, evt);
  addCollToEvent(finalCluParLinkColl, _clusterParRelCollName, evt);

  addCollToEvent(electronTightClusterColl, _electronCollName+"TightSel", evt);
  addCollToEvent(electronLooseClusterColl, _electronCollName+"LooseSel", evt);
  addCollToEvent(electronVeryLooseClusterColl, _electronCollName+"VeryLooseSel", evt);

  addCollToEvent(_cluster_start_col, "GARLICClusterStartPositions", evt);

  addCollToEvent(trackClusterColl, "GARLICTrackClusters", evt);

  addCollToEvent(MCphotonClusterColl, _clusterCollName+"MCphoton", evt);
  addCollToEvent(MCelectronClusterColl, _clusterCollName+"MCelectron", evt);
  addCollToEvent(MCchHadClusterColl, _clusterCollName+"MCchHad", evt);
  addCollToEvent(MCneuHadClusterColl, _clusterCollName+"MCneuHad", evt);
  addCollToEvent(MCotherClusterColl, _clusterCollName+"MCother", evt);
  addCollToEvent(wronglyRejectedClusterColl, _clusterCollName+"WronglyRejected", evt);
  addCollToEvent(wronglySelectedClusterColl, _clusterCollName+"WronglySelected", evt);
  addCollToEvent(seedCoreLink, _seedCoreRelCollName, evt);
  addCollToEvent(seedFinalClusLink, _seedClusterRelCollName, evt);
  addCollToEvent(tightConversionColl, _conversionCollName+"Tight", evt);
  addCollToEvent(looseConversionColl, _conversionCollName+"Loose", evt);
  addCollToEvent(photonTightPFOColl, _photonPFOCollName+"TightSel", evt);
  addCollToEvent(photonLoosePFOColl, _photonPFOCollName+"LooseSel", evt);
  addCollToEvent(photonVeryLoosePFOColl, _photonPFOCollName+"VeryLooseSel", evt);
  addCollToEvent(electronTightPFOColl, _electronPFOCollName+"TightSel", evt);
  addCollToEvent(electronLoosePFOColl, _electronPFOCollName+"LooseSel", evt);
  addCollToEvent(electronVeryLoosePFOColl, _electronPFOCollName+"VeryLooseSel", evt);

  cleanup(preClusVec);
  cleanup(trackVec);

  streamlog_out ( DEBUG1 ) << "Finished Event " << evt->getEventNumber() << endl;

  return;
}


void GarlicProcessor::addCollToEvent(LCCollection* col, string colname, LCEvent* evt) {
  if (col) {
    if (evt && colname.find("IGNORE")==string::npos) {
      evt->addCollection(col, colname);
    } else {
      delete col;
    }
  }
  return;
}

void GarlicProcessor::cleanup(vector<GarlicExtendedCluster* > &preClusVec) {
  if(preClusVec.size()>0) {
    for (size_t i = 0; i < preClusVec.size(); i++) {
      GarlicExtendedCluster *a_ext_preCluster = preClusVec[i];
      if (a_ext_preCluster) {
        vector<GarlicExtendedHit* > ext_preHit_vec = *(a_ext_preCluster->getHits());
        for (size_t pHit_i=0; pHit_i<ext_preHit_vec.size(); pHit_i++) {
          GarlicExtendedHit *a_ext_pHit = dynamic_cast<GarlicExtendedHit*>(ext_preHit_vec[pHit_i]);
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

void GarlicProcessor::cleanup(vector<GarlicExtendedTrack* > &trackVec) {
  for (size_t i=0; i<trackVec.size(); i++)
    if (trackVec[i]) {delete trackVec[i]; trackVec[i]=NULL;}
  trackVec.clear();
  streamlog_out ( DEBUG2 ) << "Cleanup finished" << endl;
}



void GarlicProcessor::check(LCEvent * evt)
{
}



void GarlicProcessor::end()
{
  streamlog_out ( MESSAGE ) << endl
                            << "GarlicProcessor Report: " << endl
                            << _nEvents << " events processed" << endl
                            << "Found " << _nClusters << " Cluster" << endl;
  if(_nEvents>0)
    streamlog_out ( MESSAGE ) << "= " << ((float)_nClusters)/_nEvents  << " /event" << endl;

  if (_fhistos) {
    _fhistos->Write(0);
    _fhistos->Close();
  }

  delete _convFinder;

}



void GarlicProcessor::PrepareMCTracks(vector<GarlicExtendedTrack* > &trackVec)
{
  streamlog_out ( DEBUG2 ) << "Getting MC tracks..." << endl
                           << "BField: " << GarlicGeometryParameters::Instance().Get_bField() << " T" << endl;

  if (!_mcParticleColl) {
    streamlog_out ( DEBUG1 ) << "no MC collection found..." << endl;
    return;
  }

  for(int jmcp=0; jmcp<_mcParticleColl->getNumberOfElements(); jmcp++) {
    MCParticle* mcp = dynamic_cast <MCParticle*> (_mcParticleColl->getElementAt(jmcp));
    if (mcp && mcp->getCharge()!=0 && !mcp->isDecayedInTracker()) {
      GarlicExtendedTrack* trk = new GarlicExtendedTrack(mcp);
      trackVec.push_back(trk);
    }
  }
  streamlog_out ( DEBUG2 ) << "Have prepared " << trackVec.size() << " MC tracks" << endl;
}


void GarlicProcessor::PrepareTracks(const LCEvent *evt, vector<GarlicExtendedTrack* > &trackVec)
{
  streamlog_out ( DEBUG2 ) << "Getting tracks..." << endl
                           << "BField: " << GarlicGeometryParameters::Instance().Get_bField() << " T" << endl;

  TrackVec allSubTracks;

  const float hitDistCutBarrel=200; // mm - how close track hits should be to ECAL front face in order to use it to to veto ecal hits
  const float hitDistCutEndcap=300; // mm - larger in endcap to deal with large gap for TPC readout and cables etc

  try {
    LCCollection* trackColl = evt->getCollection(_TrackCollectionName);
    streamlog_out ( DEBUG2 )  << trackColl->getNumberOfElements() << " tracks in " << _TrackCollectionName << " collection" << endl;
    for(int track_i=0; track_i<trackColl->getNumberOfElements(); track_i++) {
      streamlog_out ( DEBUG2 ) << "Track " << track_i << endl;
      Track *a_track = dynamic_cast<Track*>(trackColl->getElementAt(track_i));

      // check that it has hits not too far from ECAL....
      bool hitsNearEcal=false;
      const TrackerHitVec hits = a_track->getTrackerHits();
      float mostForwardHit[3]={0,0,-999};
      for (size_t i=hits.size()-1; i>0; i--) {

        if ( fabs(hits[i]->getPosition()[2]) > fabs(mostForwardHit[2]) ) {
          for (int kk=0; kk<3; kk++) {
            mostForwardHit[kk] = hits[i]->getPosition()[kk];
          }
        }

        float zdist = GarlicGeometryParameters::Instance().Get_zOfBarrel()-fabs(hits[i]->getPosition()[2]);
        if ( zdist < hitDistCutEndcap ) {
          hitsNearEcal=true;
          break;
        }
        float rdist(0);
        for (int j=0; j<2; j++)
          rdist+=pow(hits[i]->getPosition()[j],2);
        rdist=sqrt(rdist);
        rdist = GarlicGeometryParameters::Instance().Get_rOfBarrel()-rdist;
        if ( rdist < hitDistCutBarrel ) {
          hitsNearEcal=true;
          break;
        }
      }

      // if the track is very forward, allow a larger distance: we don't expect any tpc hits
      // if the pt is low, track curls, may get broken
      if ( !hitsNearEcal ) {
        float mod(0);
        for (int kk=0; kk<3; kk++) {
          mod+=pow(mostForwardHit[kk],2);
        }
        mod=sqrt(mod);
        float angle = acos(fabs(mostForwardHit[2])/mod);
        if ( angle < GarlicAlgorithmParameters::Instance().GetForwardTrackAngle() ) {
          hitsNearEcal=true;
          streamlog_out ( DEBUG2 ) << "..accepting forward track " << angle << endl;
        }

        float radiusOfCurv = 1./fabs(a_track->getOmega());
        if ( radiusOfCurv < GarlicGeometryParameters::Instance().Get_rMaxOfBarrel() ) {
          hitsNearEcal=true;
          streamlog_out ( DEBUG2 ) << "...accpting lowPt track: radius = " << radiusOfCurv << endl;
        }
      }

      // this often cuts broken curlers...
      if ( !hitsNearEcal ) {
        streamlog_out ( DEBUG2 ) << "WARNING, rejecting track with no hit near ECAL!" << endl;
        for (size_t i=hits.size()-1; i>0; i--) {
          streamlog_out ( DEBUG2 ) << "hit " << i << " : ";
          for (int j=0; j<3; j++)
            streamlog_out ( DEBUG2 ) << hits[i]->getPosition()[j] << " ";
          streamlog_out ( DEBUG2 ) << endl;
        }
        continue;
      }

      TrackVec subtracks = a_track->getTracks();
      for (size_t f=0; f<subtracks.size(); f++) {
        allSubTracks.push_back(subtracks[f]);
      }


      GarlicExtendedTrack *a_ext_track = new GarlicExtendedTrack(a_track);

      streamlog_out ( DEBUG2 ) << "making new ex trk from " << _TrackCollectionName << " " << a_ext_track << " " << a_track->getOmega() << " " <<
        a_ext_track->getMomentum()[0] << " " << a_ext_track->getMomentum()[1] << " " << a_ext_track->getMomentum()[2] << endl;

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

      streamlog_out ( DEBUG2 ) << "making new ex trk from " << _TPCTrackCollectionName << " " << a_track->getOmega() << endl;

      GarlicExtendedTrack *a_ext_track = new GarlicExtendedTrack(a_track);

      streamlog_out ( DEBUG2 ) << "made! " << a_ext_track << endl;

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



void GarlicProcessor::PreparePreClusters(LCEvent *evt, vector<GarlicExtendedCluster* > &preClusVec) {

  LCCollection *preClusColl = 0;
  try {
    preClusColl = evt->getCollection(GarlicAlgorithmParameters::Instance().GetEcalPreClusterCollectionName());
  } catch(DataNotAvailableException &exc) {
    if (GarlicAlgorithmParameters::Instance().GetDebug()>1)
      std::cerr << "In GarlicProcessor::preparePreClusters " << exc.what() << std::endl;
    return;
  }

  int NPreCluster = preClusColl->getNumberOfElements();
  if (!GarlicGeometryParameters::Instance().Get_defaultDecoder() && NPreCluster>0) {
    CellIDDecoder<CalorimeterHit>* dec = new CellIDDecoder<CalorimeterHit> (preClusColl);
    GarlicGeometryParameters::Instance().Set_defaultDecoder(dec);
  }
  streamlog_out ( DEBUG2 ) << "Number of PreClusters in original collection: " << preClusColl->getNumberOfElements() << endl;

  for (int preC_i = 0; preC_i < NPreCluster; preC_i++) {
    ClusterImpl *a_cluster = dynamic_cast<ClusterImpl* >( preClusColl->getElementAt(preC_i) );

    GarlicExtendedCluster *a_ext_cluster = new GarlicExtendedCluster( );
    const CalorimeterHitVec &preHitVec=a_cluster->getCalorimeterHits();
    streamlog_out ( DEBUG2 ) << "Reading PreCluster " << preC_i << " " << preHitVec.size() << endl;
    vector<GarlicExtendedHit *> extHitVec;

    for (EVENT::CalorimeterHitVec::const_iterator hit_it=preHitVec.begin(); hit_it!=preHitVec.end(); hit_it++) {
      CalorimeterHitImpl *a_hit = dynamic_cast<CalorimeterHitImpl*>( *hit_it );
      GarlicExtendedHit *a_ext_hit = new GarlicExtendedHit(a_hit);
      a_ext_hit->setCluster(NULL);
      a_ext_hit->setPreshower(0);
      extHitVec.push_back(a_ext_hit);
    }

    // sort hits by pseudo layer
    sort(extHitVec.begin(),extHitVec.end(),GarlicExtendedHit::lowerPseudoLayer);
    a_ext_cluster->addHits(extHitVec);
    preClusVec.push_back(a_ext_cluster);
  }

  return;
}

void GarlicProcessor::RemoveElectronHits(vector < pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > electrons, GarlicExtendedCluster* preClus ) {
  for (size_t iel=0; iel<electrons.size(); iel++) {
    std::vector <GarlicExtendedHit*> * electronHits = electrons[iel].second->getHits();
    //    cout << "removing electron hits " << electronHits->size() << endl;
    for (size_t ih=0; ih<electronHits->size(); ih++) {
      preClus->removeHit( electronHits->at(ih) );
      //cout << "removing electron hit " << ih << " " << electronHits->at(ih)->getCaloHit()->getEnergy() << " " << electronHits->at(ih)->getCaloHit()->getCellID0() << endl;
    }
  }
  return;
}


//void
std::map < GarlicExtendedTrack*, GarlicExtendedCluster*> GarlicProcessor::RemoveHitsNearExtrapolatedTracks(vector<GarlicExtendedTrack* > &trackVec,
                                                                                           vector<GarlicExtendedCluster*> &preClusVec) {

  // find hits near track extrapolations.
  // put them into a cluster (later used to check if it's mip-like or showering)
  // remove hits from preclusters

  float min_dist_cut = GarlicAlgorithmParameters::Instance().GetTrackVetoWindow();
  float dist_cut = max( GarlicGeometryParameters::Instance().Get_padSizeEcal()[1], min_dist_cut);

  std::map < GarlicExtendedTrack*, GarlicExtendedCluster*> track_clusters;

  for (size_t track_i=0; track_i<trackVec.size(); track_i++) {

    streamlog_out ( DEBUG2 )
      << " Checking track " << track_i << " : " ;

    GarlicExtendedTrack *a_ext_track = dynamic_cast<GarlicExtendedTrack* >(trackVec[track_i]);

    if ( a_ext_track->getElectronSel() > 0 ) {
      streamlog_out ( DEBUG2 ) << "it's an electron track, ignore..." << endl;
      continue;
    }

    GarlicExtendedCluster* track_cluster = new GarlicExtendedCluster();
    track_cluster->setAssociatedTrack(a_ext_track);

    // check that extrapolation to ecal looks reasonably well controlled
    if (a_ext_track->getMinDist_HitEcalEntry() > 400) {
      streamlog_out ( WARNING ) <<
        "weird looking track, don't use to veto hits..." << endl;
      continue;
    }

    for(size_t pclus_i=0; pclus_i<preClusVec.size(); pclus_i++) {
      streamlog_out ( DEBUG2 ) << "Removing hits from PreCluster " << pclus_i << endl;
      vector<GarlicExtendedHit* > *hit_vec = preClusVec[pclus_i]->getHits();
      int nhits = hit_vec->size();
      for (int ih=nhits-1; ih>=0; ih--) {
        GarlicExtendedHit *a_ext_hit = dynamic_cast<GarlicExtendedHit* > (hit_vec->at(ih));

        float dist = a_ext_track->getDistanceToPoint( a_ext_hit->getCaloHit()->getPosition() );
        if(dist<dist_cut) {
          hit_vec->erase(hit_vec->begin()+ih);
          track_cluster->addHit(a_ext_hit);
        } else {
          //streamlog_out ( DEBUG2 ) << "Hit not removed: " << a_ext_hit->getCaloHit()->getCellID0() <<
          //  ", with dist = " << dist << " while asked: "<< dist_cut <<  endl;
        }
      }
      streamlog_out ( DEBUG2 ) << "Done with PreCluster " << pclus_i << endl;
    }

    if ( track_cluster->getHits()->size()>0 )
      track_clusters[a_ext_track]=track_cluster;
    else
      delete track_cluster;

  }

  return track_clusters;
}


vector < GarlicExtendedTrack* > GarlicProcessor::selectNearbyTracks(GarlicExtendedCluster* preCluster,  vector <GarlicExtendedTrack* >* trackVec) {
  float distCut = 20*preCluster->getWidth(); // be really generous to make sure we get all tracks
  vector < GarlicExtendedTrack* > nearbytracks;
  for (size_t it=0; it<trackVec->size(); it++) {
    GarlicExtendedTrack* trk = trackVec->at(it);
    float dist(0);
    for (int j=0; j<3; j++) {
      dist += pow ( trk->getEcalEntryPos()[j] - preCluster->getProjectedPosition()[j] , 2);
    }
    dist = sqrt(dist);
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



void GarlicProcessor::printMrGarlic() {

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


void GarlicProcessor::setUpGeometry() {

  setup();

  // size of dead zones, ncells, si thickness, (not in GEAR file, specified in steering file...)
  GarlicGeometryParameters::Instance().Set_defaultDecoder  (NULL);

  // b field
  GarlicGeometryParameters::Instance().Set_bField( Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z() );

  // TPC
  // GarlicGeometryParameters::Instance().Set_tpcRmin    ( float x);
  // GarlicGeometryParameters::Instance().Set_tpcRmax    ( float x);
  // GarlicGeometryParameters::Instance().Set_tpcZmax    ( float x);
  // GarlicGeometryParameters::Instance().Set_tpcPadRows ( int   x);



  // Calorimeter geometry from GEAR
  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
  const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();

  GarlicGeometryParameters::Instance().Set_symmetry ( pEcalBarrel.getSymmetryOrder() );

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
  GarlicGeometryParameters::Instance().Set_barrelStaveDir(barrelStaveDir);

  // radius of barrel front face
  GarlicGeometryParameters::Instance().Set_rOfBarrel ( pEcalBarrel.getExtent()[0] );

  // maximum radius of ECAL
  // getExtent()[1] is the radius of last sensitive layer
  // hits in the corner can be beyond last sensitive layer "radius"
  // need to add the thickness of outer plate, etc, assume ~30mm
  GarlicGeometryParameters::Instance().Set_rMaxOfBarrel ( pEcalBarrel.getExtent()[1] + 30.);

  // z of barrel end
  GarlicGeometryParameters::Instance().Set_zOfBarrel ( pEcalBarrel.getExtent()[3] );

  // inner r of endcap
  GarlicGeometryParameters::Instance().Set_rInnerEcalEndcap ( pEcalEndcap.getExtent()[0] );

  // outer r of endcap
  GarlicGeometryParameters::Instance().Set_rOfEndcap ( pEcalEndcap.getExtent()[1] );

  // z of endcap front face
  GarlicGeometryParameters::Instance().Set_zOfEndcap ( pEcalEndcap.getExtent()[2] );

  // number of layers (include the preshower layer)
  //GarlicGeometryParameters::Instance().Set_nBarrelEcalLayers( ecalBarrelLayout.getNLayers()+1 );
  //GarlicGeometryParameters::Instance().Set_nEndcapEcalLayers( ecalEndcapLayout.getNLayers()+1 );
  //GarlicGeometryParameters::Instance().Set_nPseudoLayers( max(ecalBarrelLayout.getNLayers(), ecalEndcapLayout.getNLayers()) + 1 );
  // for mokka fix, aug2014
  GarlicGeometryParameters::Instance().Set_nBarrelEcalLayers( ecalBarrelLayout.getNLayers() );
  GarlicGeometryParameters::Instance().Set_nEndcapEcalLayers( ecalEndcapLayout.getNLayers() );
  GarlicGeometryParameters::Instance().Set_nPseudoLayers( max(ecalBarrelLayout.getNLayers(), ecalEndcapLayout.getNLayers()) + 1 );

  // layer0 is "preshower" layer
  // layer1 is "first" ECAL layer
  // absorber thickness arrays in barrel and endcap

  float absThick[MAX_NUMBER_OF_LAYERS]={0};
  float cellSize[MAX_NUMBER_OF_LAYERS]={0};
  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    absThick[i] = ecalBarrelLayout.getAbsorberThickness(i);
    cellSize[i] = ecalBarrelLayout.getCellSize0(i);
  }
  GarlicGeometryParameters::Instance().Set_absThicknessBarrelLayer(absThick);
  GarlicGeometryParameters::Instance().Set_padSizeEcal(cellSize);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) absThick[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    absThick[i] = ecalEndcapLayout.getAbsorberThickness(i);
  }
  GarlicGeometryParameters::Instance().Set_absThicknessEndcapLayer(absThick);

  // layer positions: this should be approx position of centre of silicon layer
  float positions[MAX_NUMBER_OF_LAYERS]={0};

  // if the parameter Ecal_barrel_gear_per_sensitiveLayer is set in the gear file, 
  //     what we get from ecalBarrelLayout.getDistance() is ~ the centre of the silicon layer
  // otherwise, it is ~centre of absorber layer
  int getLayPosBySensitive=0;
  try {
    getLayPosBySensitive=pEcalBarrel.getIntVal("Ecal_barrel_gear_per_sensitiveLayer");
  }
  catch(gear::UnknownParameterException) {}

  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    if ( getLayPosBySensitive==1 ) {
      positions[i] = ecalBarrelLayout.getDistance(i); // new defn of layer position
    } else {
      positions[i] = ecalBarrelLayout.getDistance(i) - ecalBarrelLayout.getThickness(i)/2; // the older one (e.g. v05 models)
    }
  }
  if ( getLayPosBySensitive!=1 ) { // we need some messy trickery for the last sensitive layer...
    positions[ecalBarrelLayout.getNLayers()] = positions[ecalBarrelLayout.getNLayers()-1] + ecalBarrelLayout.getThickness(ecalBarrelLayout.getNLayers()-1);
  }
  GarlicGeometryParameters::Instance().Set_positionBarrelLayer(positions);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) positions[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    if ( getLayPosBySensitive==1 ) {
      positions[i] = ecalEndcapLayout.getDistance(i); // + ecalEndcapLayout.getThickness(i)/2; for mokka fix aug2014
    } else {
      positions[i] = ecalEndcapLayout.getDistance(i) + ecalEndcapLayout.getThickness(i)/2;
    }
  }
  if ( getLayPosBySensitive!=1 ) { // we need some messy trickery for the last sensitive layer...
    positions[ecalEndcapLayout.getNLayers()] = positions[ecalEndcapLayout.getNLayers()-1] + ecalEndcapLayout.getThickness(ecalEndcapLayout.getNLayers()-1);
  }
  GarlicGeometryParameters::Instance().Set_positionEndcapLayer(positions);

  float rad_length=3.5; // mm, for tungsten. if you want to test a different radiator, you'll need to expand this part
  if ( find( pEcalBarrel.getStringKeys().begin(), pEcalBarrel.getStringKeys().end(), "Ecal_radiator_material" )!=pEcalBarrel.getStringKeys().end() ) {
    if ( pEcalBarrel.getStringVal ( "Ecal_radiator_material" )!="tungsten" )
      streamlog_out( WARNING ) << "do not know radiation length of " << pEcalBarrel.getStringVal ( "Ecal_radiator_material" ) << " , assuming " << rad_length << " mm " << endl;
  }
  GarlicGeometryParameters::Instance().Set_absorberRadiationLength(rad_length);

  return;
}

