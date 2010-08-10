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

#include "ECALGarlicClusterHelpers.hh"

#include "ECALGarlicCluster.hh"

#include "ECALGarlicEnergyEstimator.hh"

#include "ECALGarlicConstants.hh"
using namespace ECALGarlicConstants;


#include <Rtypes.h>

// From MarlinUtil
//For Tracks and Track projections
#include <HelixClass.h>
#include <ECALGarlic.hh>

#include "ECALGarlicGeometryParameters.hh"

using namespace marlin;
using namespace lcio;
using namespace std;

ECALGarlic anECALGarlic;

ECALGarlic::ECALGarlic() 
  : Processor("ECALGarlic") {
  // processor description
  _description = "Clustering and photon recognition";

  string mcParticleCollName = "MCParticle";
  registerProcessorParameter("MCParticleCollection",
			     "MC Particle collection name",
			     _mcParticleCollectionName,
			     mcParticleCollName);

  string ecalBarrelPreShowerHitCollName = "EcalBarrelPreShowerCollection";
  registerProcessorParameter("EcalBarrelPreShowerHitCollection",
			     "PreShower Hit collection name of ECAL Barrel",
			     _ecalBarrelPreShowerHitCollectionName,
			     ecalBarrelPreShowerHitCollName);

  string ecalEndcapPreShowerHitCollName = "EcalEndcapPreShowerCollection";
  registerProcessorParameter("EcalEndcapPreShowerHitCollection",
			     "PreShower Hit collection name of ECAL Endcap",
			     _ecalEndcapPreShowerHitCollectionName,
			     ecalEndcapPreShowerHitCollName);


  // input pre-clusters
  string ecalPreClusterCollName = "ECAL_PreClusters";
  registerProcessorParameter("EcalPreClusterCollection",
			     "Name of PreCluster collection",
			     _x_ecalPreClusterCollectionName,
			     ecalPreClusterCollName);

  // input track collection
  string LDCTrackCollName = "LDCTracks";
  registerProcessorParameter( "LDCTrackCollection",
			   "LDC track collection name",
			   _LDCTrackCollectionName,
			   LDCTrackCollName);

  // output collections
  string ecalPhotonClusterCollName = "ECAL_PhotonClusters";
  registerProcessorParameter("EcalPhotonClusterCollection",
			     "Name of Photon Cluster collection",
			     _ecalPhotonClusterCollectionName,
			     ecalPhotonClusterCollName);

  string RemovedHitsCollName = "RemovedHits";
  registerProcessorParameter( "RemovedHitsCollection",
			   "collection name of removed hits (near tracks)",
			   _removedHitsCollectionName,
			   RemovedHitsCollName);

  string ParticleCollName = "GARLIC_Photons";
  registerProcessorParameter( "ParticleCollName",
			   "collection name for the reconstructed particles",
			   _particleCollectionName,
			      ParticleCollName);


  // parameters

  registerProcessorParameter("NLayersForSeeding",
			     "Number of ECAL pseudo layers used for projecting to obtain a seed, typically equivalent to 5 X0.",
			     _x_nLayersForSeeding,
			     int(8));

  registerProcessorParameter("MinHitsForSeeding",
			     "Minimum number of hits to accept a seed.",
			     _x_minHitsForSeeding,
			     int(4));

  registerProcessorParameter("MinHits10X0",
			     "Minimum number of hits after 10X0.",
			     _x_minHits10X0,
			     int(5));

 registerProcessorParameter("NIterations",
			     "Number of Iterations to apply the neighbpuring criterion.",
			     _x_nIterations,
			     int(3));

 registerProcessorParameter("DistanceFactor",
			     "Factor applied to CellSize used for the maximal colecting distance.",
			     _x_distanceFactor,
			     double(1.41422));

  registerProcessorParameter("MinClusterHits",
			     "Minimum Number of Hits to form a cluster",
			     _x_nHitsMin,
			     int(5));

  registerProcessorParameter("MinClusterEnergy",
			     "Minimum energy form a cluster",
			     _x_minEnergy,
			     float(0.15));

  registerProcessorParameter("ApplyGapCorrection",
			     "Try to correct for wafer and alveola gaps or alveola gaps only?",
			     _x_applyGapCorrection,
			     true);

  registerProcessorParameter("MergeSatellites",
			     "Merge satellite clusters?",
			     _x_mergeSatellites,
			     true);

  registerProcessorParameter("IncludePreShowerHits",
			     "Try to reintegrate PreShower Hits in PreClusters",
			     _x_includePreShower,
			     true);

  registerProcessorParameter("RemoveHitsNearToTracks",
			     "Should remove Hits near extrapolated tracks to reject pions?",
			     _x_removeHitsNearTracks,
			     true);

  registerProcessorParameter("RejectClustersByMLPCut",
			     "Should remove clusters that do not pass MLP cuts?",
			     _x_rejectMLP,
			     true);

  registerProcessorParameter("CorrectLeakage",
			     "Include correction for longitudinal leakage?",
			     _x_correctLeakage,
			     false);

  registerProcessorParameter("CorrectPhi",
			     "Include correction for phi?",
			     _x_correctPhi,
			     true);

  registerProcessorParameter("CorrectTheta",
			     "Include correction for theta?",
			     _x_correctTheta,
			     true);

  registerProcessorParameter("CheatTracks",
			     "Initialise HelixFit with MC particle momentum (1), or read from text file (2), or do nothing(0) ",
			     _x_cheatTracks,
			     int(0));

  registerProcessorParameter("XYGapTransitFactor",
			     "In order to correct for the influence of the magnetic field the shower direction in x and y will be corrected by this factor when following a cluster over the cable gap between barrel and endcap",
			     _x_xy_gap_transit_factor,
			     float(1));


  // parameters for photon energy estimation
  //  including corrections in theta and phi

  _x_toGeVParameters.clear();
  _x_toGeVParameters.push_back(4.47401e-02);
  _x_toGeVParameters.push_back(1.18857);
  _x_toGeVParameters.push_back(-1.81962e-03);
  _x_toGeVParameters.push_back(1.00757e-04);
  registerProcessorParameter("ToGeVParameters",
			     "The parameters for the function to go from energy deposit to GeV value",
			     _x_toGeVParameters,
			     _x_toGeVParameters,
			     _x_toGeVParameters.size() 
			     );

  _x_toGeVParameters_EC.clear();
  _x_toGeVParameters_EC.push_back(4.44222e-02);
  _x_toGeVParameters_EC.push_back(1.20706);
  _x_toGeVParameters_EC.push_back(-1.77230e-03);
  _x_toGeVParameters_EC.push_back(9.47984e-05);

  registerProcessorParameter("ToGeVParameters_EC",
			     "The parameters_EC for the function to go from energy deposit to GeV value",
			     _x_toGeVParameters_EC,
			     _x_toGeVParameters_EC,
			     _x_toGeVParameters_EC.size() 
			     );

  _x_corrThParameters0.clear();
  _x_corrThParameters0.push_back(-2.415974e-02);
  _x_corrThParameters0.push_back(-3.106519e-04);
  registerProcessorParameter("ThetaCorrPar0",
			     "The parameters for the function to estimate par0 of the function for theta correction",
			     _x_corrThParameters0,
			     _x_corrThParameters0,
			     _x_corrThParameters0.size() 
			     );

  _x_corrThParameters1.clear();
  _x_corrThParameters1.push_back(-2.512019e-02);
  _x_corrThParameters1.push_back(1.253120e-04);
  registerProcessorParameter("ThetaCorrPar1",
			     "The parameters for the function to estimate par1 of the function for theta correction",
			     _x_corrThParameters1,
			     _x_corrThParameters1,
			     _x_corrThParameters1.size() 
			     );



  _x_corrPhiParameters0.clear();
  _x_corrPhiParameters0.push_back(2.22114);
  _x_corrPhiParameters0.push_back(-5.96095e-01);
  _x_corrPhiParameters0.push_back(-4.85018e-07);
  _x_corrPhiParameters0.push_back(2.44816e-04);
  _x_corrPhiParameters0.push_back(-4.98479);
  _x_corrPhiParameters0.push_back(2.57302e-03);
  registerProcessorParameter("PhiCorrPar0",
			     "Phie parameters for Phie function to estimate par0 of Phie function for Phieta correction",
			     _x_corrPhiParameters0,
			     _x_corrPhiParameters0,
			     _x_corrPhiParameters0.size() 
			     );

  _x_corrPhiParameters1.clear();
  _x_corrPhiParameters1.push_back(3.93864e+02);
  _x_corrPhiParameters1.push_back(9.50289e+02);
  _x_corrPhiParameters1.push_back(3.12544e-05);
  _x_corrPhiParameters1.push_back(8.61194e-03);
  registerProcessorParameter("PhiCorrPar1",
			     "Phie parameters for Phie function to estimate par0 of Phie function for Phieta correction",
			     _x_corrPhiParameters1,
			     _x_corrPhiParameters1,
			     _x_corrPhiParameters1.size() 
			     );

  _x_corrPhiParameters2.clear();
  _x_corrPhiParameters2.push_back(1.53569e+03);
  _x_corrPhiParameters2.push_back(1.63601e+05);
  _x_corrPhiParameters2.push_back(-1.54682e-05);
  _x_corrPhiParameters2.push_back(1.63612e-01);
  registerProcessorParameter("PhiCorrPar2",
			     "Phie parameters for Phie function to estimate par1 of Phie function for Phieta correction",
			     _x_corrPhiParameters2,
			     _x_corrPhiParameters2,
			     _x_corrPhiParameters2.size() 
			     );

  _x_mlpCuts_B.clear();

  _x_mlpCuts_B.push_back(-.23);
  _x_mlpCuts_B.push_back(-.3);
  _x_mlpCuts_B.push_back(-.4); 
  _x_mlpCuts_B.push_back(-.2);
  _x_mlpCuts_B.push_back(-.1);
  _x_mlpCuts_B.push_back(-.2);
  _x_mlpCuts_B.push_back(-.25);
  _x_mlpCuts_B.push_back(-.2);
  _x_mlpCuts_B.push_back(-.4);
  _x_mlpCuts_B.push_back(0);
  _x_mlpCuts_B.push_back(-.3);
  _x_mlpCuts_B.push_back(.1);
  _x_mlpCuts_B.push_back(.5);
  registerProcessorParameter("MLPCutsBarrel",
			     "Cuts applied on the NN decision",
			     _x_mlpCuts_B,
			     _x_mlpCuts_B,
			     _x_mlpCuts_B.size() 
			     );

  _x_mlpCuts_EC.clear();
  _x_mlpCuts_EC.push_back(-.6);
  _x_mlpCuts_EC.push_back(-.5);
  _x_mlpCuts_EC.push_back(-.5);
  _x_mlpCuts_EC.push_back(-.3);
  _x_mlpCuts_EC.push_back(-.2);
  _x_mlpCuts_EC.push_back(-.35);
  _x_mlpCuts_EC.push_back(-.4);
  _x_mlpCuts_EC.push_back(0);
  _x_mlpCuts_EC.push_back(0);
  _x_mlpCuts_EC.push_back(0);
  _x_mlpCuts_EC.push_back(.1);
  _x_mlpCuts_EC.push_back(-.3);
  _x_mlpCuts_EC.push_back(.4);
  registerProcessorParameter("MLPCutsEndcap",
			     "Cuts applied on the NN decision",
			     _x_mlpCuts_EC,
			     _x_mlpCuts_EC,
			     _x_mlpCuts_EC.size() 
			     );


  _x_alp_params_B.push_back(6.10498);
  _x_alp_params_B.push_back(8.99890e-01);
  _x_alp_params_B.push_back(-1.80243);
  _x_alp_params_B.push_back(2.10402e+01);
  _x_alp_params_B.push_back(1.22358);
  _x_alp_params_B.push_back(2.26219);
  _x_alp_params_B.push_back(-2.52431);
  _x_alp_params_B.push_back(7.70250);
  registerProcessorParameter("AlphaBarrel",
			     "alpha parameter for barrel",
			     _x_alp_params_B,
			     _x_alp_params_B,
			     _x_alp_params_B.size()
			     );

  _x_bet_params_B.push_back(2.25252);
  _x_bet_params_B.push_back(7.44998e-01);
  _x_bet_params_B.push_back(5.87399e-01);
  registerProcessorParameter("BetaBarrel",
			     "beta parameter for barrel",
			     _x_bet_params_B,
			     _x_bet_params_B,
			     _x_bet_params_B.size()
			     );

  _x_g_params_B.push_back(2.52130e-02);
  _x_g_params_B.push_back(2.60046e-01);
  _x_g_params_B.push_back(3.72620e-01);
  _x_g_params_B.push_back(6.38886e-05);
  registerProcessorParameter("GammaBarrel",
			     "gamma parameter for barrel",
			     _x_g_params_B,
			     _x_g_params_B,
			     _x_g_params_B.size()
			     );

  _x_d_params_B.push_back(3.97137e-02);
  _x_d_params_B.push_back(1.69812e-01);
  _x_d_params_B.push_back(3.25349e-02);
  _x_d_params_B.push_back(2.89825e-01);
  _x_d_params_B.push_back(7.19802e-04);
  registerProcessorParameter("DeltaBarrel",
			     "delta parameter for barrel",
			     _x_d_params_B,
			     _x_d_params_B,
			     _x_d_params_B.size()
			     );

  _x_lam_params_B.push_back(2.25788);
  _x_lam_params_B.push_back(6.20427e-01);
  _x_lam_params_B.push_back(8.48215e-07);
  _x_lam_params_B.push_back(3.82175e-01);
  registerProcessorParameter("LambdaBarrel",
			     "lambda parameter for barrel",
			     _x_lam_params_B,
			     _x_lam_params_B,
			     _x_lam_params_B.size()
			     );

  _x_alp_params_EC.push_back(6.29132);
  _x_alp_params_EC.push_back(8.94917e-01);
  _x_alp_params_EC.push_back(-4.84368);
  _x_alp_params_EC.push_back(2.66329e+01);
  _x_alp_params_EC.push_back(1.03308);
  _x_alp_params_EC.push_back(2.30324);
  _x_alp_params_EC.push_back(-1.61760e-01);
  _x_alp_params_EC.push_back(8.90058);
  registerProcessorParameter("AlphaEndcap",
			     "alpha parameter for endcap",
			     _x_alp_params_EC,
			     _x_alp_params_EC,
			     _x_alp_params_EC.size()
			     );

  _x_bet_params_EC.push_back(2.36623);
  _x_bet_params_EC.push_back(8.08524e-01);
  _x_bet_params_EC.push_back(5.95689e-01);
  registerProcessorParameter("BetaEndcap",
			     "beta parameter for endcap",
			     _x_bet_params_EC,
			     _x_bet_params_EC,
			     _x_bet_params_EC.size()
			     );
  
  _x_g_params_EC.push_back(2.60573e-02);
  _x_g_params_EC.push_back(2.71003e-01);
  _x_g_params_EC.push_back(2.94263e-01);
  _x_g_params_EC.push_back(2.90116e-05);
  registerProcessorParameter("GammaEndcap",
			     "gamma parameter for endcap",
			     _x_g_params_EC,
			     _x_g_params_EC,
			     _x_g_params_EC.size()
			     );

  
  _x_d_params_EC.push_back(4.15537e-02);
  _x_d_params_EC.push_back(1.78909e-01);
  _x_d_params_EC.push_back(3.35461e-02);
  _x_d_params_EC.push_back(3.20522e-01);
  _x_d_params_EC.push_back(7.68040e-04);
  registerProcessorParameter("DeltaEndcap",
			     "delta parameter for endcap",
			     _x_d_params_EC,
			     _x_d_params_EC,
			     _x_d_params_EC.size()
			     );

  _x_lam_params_EC.push_back(1.50388);
  _x_lam_params_EC.push_back(8.17268e-01);
  _x_lam_params_EC.push_back(5.78497e-09);
  _x_lam_params_EC.push_back(5.84671e-01);
  registerProcessorParameter("LambdaEndcap",
			     "lambda parameter for endcap",
			     _x_lam_params_EC,
			     _x_lam_params_EC,
			     _x_lam_params_EC.size()
			     );



  _x_f1_params_b.push_back(0.4897884);
  _x_f1_params_b.push_back(-2.887818e-04);
  registerProcessorParameter("f1Barrel",
			     "f1 parameter for barrel",
			     _x_f1_params_b,
			     _x_f1_params_b,
			     _x_f1_params_b.size()
			     );


  _x_f2_params_b.push_back(0.4901761);
  _x_f2_params_b.push_back(-7.351988e-05);
  registerProcessorParameter("f2Barrel",
			     "f2 parameter for barrel",
			     _x_f2_params_b,
			     _x_f2_params_b,
			     _x_f2_params_b.size()
			     );

  _x_f1_params_e.push_back(0.4957159);
  _x_f1_params_e.push_back(-2.734185e-04);
  registerProcessorParameter("f1Endcap",
			     "f1 parameter for endcap",
			     _x_f1_params_e,
			     _x_f1_params_e,
			     _x_f1_params_e.size()
			     );

  _x_f2_params_e.push_back(0.4911473);
  _x_f2_params_e.push_back(-3.308976e-04);
  registerProcessorParameter("f2Endcap",
			     "f2 parameter for endcap",
			     _x_f2_params_e,
			     _x_f2_params_e,
			     _x_f2_params_e.size()
			     );


  _x_par0_f_params.push_back(34.6695);
  _x_par0_f_params.push_back(0.354309);
  _x_par0_f_params.push_back(-0.00095874);
  registerProcessorParameter("par0_f",
			     "par0f",
			     _x_par0_f_params, 
			     _x_par0_f_params, 
			     _x_par0_f_params.size()
			     );

  
  _x_par1_f_params.push_back(7.19559);
  _x_par1_f_params.push_back(0.0453275);
  _x_par1_f_params.push_back(-0.00010602);
  registerProcessorParameter("par1_f",
			     "par1f",
			     _x_par1_f_params, 
			     _x_par1_f_params, 
			     _x_par1_f_params.size()
			     );

  registerProcessorParameter("MinSatelliteEnergy",
			     "Lower Energy cut to define a satellite cluster",
			     _x_minSatelliteEn,
			     float(0.15));

  registerProcessorParameter("MaxSatelliteEnergy",
			     "Upper Energy cut to define a satellite cluster",
			     _x_maxSatelliteEn,
			     float(3.));

  registerProcessorParameter("DebugMode",
			     "Talk a lot? (0-3)",
			     _x_debug,
			     int(0));

  // detector parameters
  // should this not be taken from the gear file???

  registerProcessorParameter("GuardringSize",
			     "Size of the guard ring in mm",
			     _x_guardringSize,
			     float(0.5));

  registerProcessorParameter("FiberSize",
			     "Size of the fiber in mm",
			     _x_fiberSize,
			     float(0.95));

  registerProcessorParameter("FiberSizeModule",
			     "Size of the guard ring in mm",
			     _x_fiberSizeModule,
			     float(2.));

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

  registerProcessorParameter("NCellsPerWafer",
			     "Number of cells along one edge of a square wafer",
			     _x_nCellsPerWafer,
			     int(18));

}


void ECALGarlic::init()
{
  printMrGarlic();
  printParameters();


  _garPars = new ECALGarlicAlgorithmParameters();
  
  _garPars->SetEcalPreClusterCollectionName (_x_ecalPreClusterCollectionName);
  _garPars->SetNLayersForSeeding (_x_nLayersForSeeding);
  _garPars->SetCorrectPhi (_x_correctPhi);
  _garPars->SetCorrectTheta (_x_correctTheta);
  _garPars->SetDebug (_x_debug);
  _garPars->SetXY_gap_transit_factor(_x_xy_gap_transit_factor);
  _garPars->SetNHitsMin(_x_nHitsMin);
  _garPars->SetMinEnergy(_x_minEnergy);
  _garPars->SetDistanceFactor(_x_distanceFactor);
  _garPars->SetMinHits10X0(_x_minHits10X0);
  _garPars->SetNIterations(_x_nIterations);
  _garPars->SetMaxSatEn(_x_maxSatelliteEn);
  _garPars->SetMinSatEn(_x_minSatelliteEn);
  _garPars->SetApplyGapCorr(_x_applyGapCorrection);
  _garPars->SetRemoveHitsNearTracks(_x_removeHitsNearTracks);
  _garPars->SetRejectMLP(_x_rejectMLP);
  _garPars->SetMerSat(_x_mergeSatellites);
  _garPars->SetCheatTracks(_x_cheatTracks);
  _garPars->SetIncludePreshower(_x_includePreShower);
  _garPars->SetMinHitsForSeed(_x_minHitsForSeeding);
  _garPars->SetCorrectLeakage(_x_correctLeakage);


  _garPars->SetToGeVParsBarrel (_x_toGeVParameters);
  _garPars->SetToGeVParsEndcap (_x_toGeVParameters_EC);
  _garPars->SetMLPCutsBarrel   (_x_mlpCuts_B);
  _garPars->SetMLPCutsEndcap   (_x_mlpCuts_EC);
  _garPars->SetCorrThetaPars0  (_x_corrThParameters0);
  _garPars->SetCorrThetaPars1  (_x_corrThParameters1);
  _garPars->SetCorrPhiPars0    (_x_corrPhiParameters0);
  _garPars->SetCorrPhiPars1    (_x_corrPhiParameters1);
  _garPars->SetCorrPhiPars2    (_x_corrPhiParameters2);
  _garPars->SetAlphaParsBarrel (_x_alp_params_B);
  _garPars->SetAlphaParsEndcap (_x_alp_params_EC);
  _garPars->SetBetaParsBarrel  (_x_bet_params_B);
  _garPars->SetBetaParsEndcap  (_x_bet_params_EC);
  _garPars->SetGammaParsBarrel (_x_g_params_B);
  _garPars->SetGammaParsEndcap (_x_g_params_EC);
  _garPars->SetDeltaParsBarrel (_x_d_params_B);
  _garPars->SetDeltaParsEndcap (_x_d_params_EC);
  _garPars->SetLambdaParsBarrel(_x_lam_params_B);
  _garPars->SetLambdaParsEndcap(_x_lam_params_EC);
  _garPars->SetF1ParsBarrel    (_x_f1_params_b);
  _garPars->SetF2ParsBarrel    (_x_f2_params_b);
  _garPars->SetF1ParsEndcap    (_x_f1_params_e);
  _garPars->SetF2ParsEndcap    (_x_f2_params_e);
  _garPars->SetPar0FPars       (_x_par0_f_params);
  _garPars->SetPar1FPars       (_x_par1_f_params);



  _geomPars = new ECALGarlicGeometryParameters();



  _geomHelpers = new ECALGarlicGeometryHelpers(_garPars, _geomPars);

  _clusterHelpers = new ECALGarlicClusterHelpers(_garPars, _geomPars);

  _clusterer = new ECALGarlicCluster(_garPars, _geomPars);

  _energyEstimator = new ECALGarlicEnergyEstimator(_garPars, _geomPars);

  _nEvents=0;
  _nPhotonClusters=0;
  // set up NN variables

  return;
}



void ECALGarlic::processRunHeader(LCRunHeader * run)
{

  setUpGeometry();

  _geomPars->Set_rOfBarrel (_x_rOfBarrel);
  _geomPars->Set_zOfBarrel (_x_zOfBarrel);
  _geomPars->Set_rOfEndcap (_x_rOfEndcap);
  _geomPars->Set_zOfEndcap (_x_zOfEndcap);

  _geomPars->Set_activeThickness        (_x_activeThickness);
  _geomPars->Set_passiveThickness       (_x_passiveThickness);
  _geomPars->Set_firstEndcapLayerOffset (_x_firstEndcapLayerOffset);
  _geomPars->Set_firstBarrelLayerOffset (_x_firstBarrelLayerOffset);
  _geomPars->Set_bField                 (_x_bField);
  _geomPars->Set_rInnerEcalEndcap       (_x_rInnerEcalEndcap);
  _geomPars->Set_cosOfBarrel            (_x_cosOfBarrel);
  _geomPars->Set_guardringSize          (_x_guardringSize);
  _geomPars->Set_fiberSize              (_x_fiberSize);
  _geomPars->Set_fiberSizeModule        (_x_fiberSizeModule);
  _geomPars->Set_symmetry               (_x_symmetry);
  _geomPars->Set_nPseudoLayers          (_x_nPseudoLayers);
  _geomPars->Set_nCellsPerWafer         (_x_nCellsPerWafer);

  _geomPars->Set_positionBarrelLayer (_x_positionBarrelLayer);
  _geomPars->Set_absThicknessBarrelLayer (_x_absThicknessBarrelLayer);
  _geomPars->Set_absThicknessEndcapLayer (_x_absThicknessEndcapLayer);
  _geomPars->Set_padSizeEcal (_x_padSizeEcal);
  _geomPars->Set_positionEndcapLayer (_x_positionEndcapLayer);
  _geomPars->Set_barrelStaveDir (_x_barrelStaveDir);

  _geomPars->Set_defaultDecoder (NULL);

  return;
}


void ECALGarlic::processEvent(LCEvent * evt)   // main !
{
  if(evt->getEventNumber()%1000==0)
    cout << endl << "Event: " << evt->getEventNumber() << endl;

  _nEvents++;
  
  _mcParticleColl = 0;
  try {
    _mcParticleColl = evt->getCollection(_mcParticleCollectionName);
  }
  catch(DataNotAvailableException err) {};

  // get Track collection
  vector <ExtendedTrack* > trackVec;
  PrepareTracks(evt,trackVec);
      
  vector<ExtendedCluster* > preClusVec;
  map<ExtendedTrack*,map<int,vector<ExtendedHit*> > > allRemovedHits;

  // 1.) transform PreClusters to ExtendedClusters (pseudoIDs etc)
  if(_garPars->GetDebug()>1) cout << "Preparing PreClusters... " << endl;
  PreparePreClusters(evt, preClusVec);

  map<int,vector<ExtendedCluster*> *> cluster_map;
  IMPL::LCCollectionVec *seed_col = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  seed_col->setFlag(seed_col->getFlag()|( 1 << LCIO::RCHBIT_LONG));

  // 1b.) clear hits near to extrapolated tracks
  if(_garPars->GetRemoveHitsNearTracks()) {
    if(_garPars->GetDebug()>1) cout << "Removing hits near tracks..." << endl;
    if(trackVec.size()>0) 
      RemoveHitsNearExtrapolatedTracks(evt,trackVec,preClusVec,&allRemovedHits);
  }


  // the next steps are done per ROI (preCluster)
  int n_roi=preClusVec.size();
  if(n_roi>0) {
    if(_garPars->GetDebug()>1) cout << n_roi << " PreClusters" << endl;
 
    for(int roi_i=0;roi_i<n_roi;roi_i++) {
      ExtendedCluster *preCluster=preClusVec[roi_i];
      if(_garPars->GetDebug()>1) cout << "Building histogram for seeding for PreCluster " << roi_i << endl;

      // 2.) build histogram for seedID
      vector <vec3> possSeeds;
      ProjectForSeeding(*preCluster, possSeeds, seed_col);
    
      // 3.) cluster and check cluster criteria
      if(_garPars->GetDebug()>1) cout << "Clustering PreCluster " << roi_i << endl;
      vector<ExtendedCluster* > *clusters = new vector <ExtendedCluster*>;

      _clusterer->Cluster(*preCluster, possSeeds, *clusters, trackVec);
      int NClusters = clusters->size(); 

      if(_garPars->GetDebug()>1) cout << "Clustering PreCluster... got " << NClusters << "clusters in ROI " << roi_i << endl;

      if(_garPars->GetApplyGapCorr()) {
	for(int clus_i=0;clus_i<NClusters;clus_i++) {
	  if(_garPars->GetDebug()>1) cout << "Now applying initial gap correction to cluster " << clus_i << endl;
	  ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);
	  GhostCluster *a_ghost_clus = new GhostCluster();
	  a_ghost_clus->ofCluster=a_clus;
	  a_clus->Ghosts=a_ghost_clus;
	  _energyEstimator->ApplyFullGapCorrection(evt,a_clus,a_ghost_clus); 
	  if(_garPars->GetDebug()>1) cout << "Found " << (a_ghost_clus->ghostHitVec).size() << " ghost hits " << endl;
	}
      }

      for(int clus_i=0;clus_i<NClusters;clus_i++) {
	ClusterParameters *cPar = new ClusterParameters();
	ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);
	a_clus->Parameters = cPar;
	_clusterHelpers->FillClusterParameters(evt,a_clus,clusters,clus_i,cPar,trackVec); 
      }
      sort(clusters->begin(),clusters->end(),ExtendedCluster::higherEnergy);


      // 5.) calculate photon probabilities	
      for(int clus_i=0;clus_i<NClusters;clus_i++) {
	if(_garPars->GetDebug()>1) cout << "Calculating photon probability for cluster " << clus_i << endl;
	ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);
	ClusterParameters *clPar = a_clus->Parameters;
	_clusterHelpers->CalculatePhotonProbability(evt,a_clus,clusters,clus_i,clPar,trackVec,&allRemovedHits); 
      }

      if(_garPars->GetRejectMLP()) {
	vector<ExtendedCluster*>::iterator clus_it;
	int clus_cnt=-1;
	for(clus_it=clusters->begin();clus_it!=clusters->end();) {
	  clus_cnt++;
	  bool is_photon = RejectByMLPCut(*clus_it);
	  if(!is_photon) {
	    if(_garPars->GetDebug()>1) cout << "Cluster " <<  clus_cnt << " failed MLP cut, removing!" << endl;
	    _clusterHelpers->FreeHits(*clus_it); 
	    clusters->erase(clus_it);
	  }
	  else
	    clus_it++;
	}
      }
      NClusters = clusters->size(); 
      if(_garPars->GetDebug()>1) cout << NClusters << " Clusters in ROI " << roi_i << endl;

      if(_garPars->GetMerSat()) {
	//	  if(NClusters>1) { // why only >1??? 
	if(NClusters>0) { // changed by dj

	  if(_garPars->GetDebug()>0) cout << "Merging satellites" << endl;
	  _clusterer->MergeSatellites(clusters);
	  sort(clusters->begin(),clusters->end(),ExtendedCluster::higherEnergy);
	  NClusters = clusters->size(); // size of clusters changed

	  if(_garPars->GetDebug()>0) cout << "Having " << NClusters << " clusters after merging" << endl;

	  if(_garPars->GetApplyGapCorr()) { // why do we do this again, before and after merging???
	    for(int clus_i=0;clus_i<NClusters;clus_i++) {
	      if(_garPars->GetDebug()>1) cout << "Now applying gap correction after merging to cluster " << clus_i << endl;
	      ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);
	      GhostCluster *a_ghost_clus = a_clus->Ghosts;
	      delete a_ghost_clus;
	      GhostCluster *another_ghost_clus = new GhostCluster();
	      another_ghost_clus->ofCluster=a_clus;
	      a_clus->Ghosts=another_ghost_clus;
	      _energyEstimator->ApplyFullGapCorrection(evt,a_clus,another_ghost_clus); 
	      if(_garPars->GetDebug()>1) cout << "Found " << (another_ghost_clus->ghostHitVec).size() << " ghost hits " << endl;
	    }
	  }

	  for(int clus_i=0;clus_i<NClusters;clus_i++) {
	    ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);
	    ClusterParameters *clPar = a_clus->Parameters;
	    _clusterHelpers->CalculatePhotonProbability(evt,a_clus,clusters,clus_i,clPar,trackVec,&allRemovedHits); 
	  }

	}
      }

      NClusters = clusters->size(); 
      // recalculate Photon Probabilities
      // why we need to recalculate yet again??? dj
      for(int clus_i=0;clus_i<NClusters;clus_i++) {
	if(_garPars->GetDebug()>1) cout << "ReCalculating photon probability for cluster " << clus_i << endl;
	ExtendedCluster *a_clus=dynamic_cast<ExtendedCluster*>((*clusters)[clus_i]);	  
	ClusterParameters *clPar = a_clus->Parameters;
	_clusterHelpers->CalculatePhotonProbability(evt,a_clus,clusters,clus_i,clPar,trackVec,&allRemovedHits);
	if(_garPars->GetCorrectLeakage() && clPar->E_GeV>5.) {
	  _energyEstimator->CorrectLeakage(clPar);
	}
      }

      _nClusters+=clusters->size();
      if(clusters->size()==0) {
	clusters->clear();
      }
      cluster_map.insert(make_pair(roi_i,clusters));
      
    }
      
    // redefine cluster IDs
    map<int,vector<ExtendedCluster* > *>::iterator cluster_map_it;
    vector<ExtendedCluster* >::iterator cl_it;
    int ID = 0;
    for(cluster_map_it=cluster_map.begin();cluster_map_it!=cluster_map.end();cluster_map_it++) {
      vector<ExtendedCluster* > *clusVec = (cluster_map_it->second);
      for(cl_it=clusVec->begin();cl_it!=clusVec->end();cl_it++) {
	(*cl_it)->Parameters->ID = ID;
	ID++;
      }
    }

    // 6.) write Clusters as LC collection
    if(_garPars->GetDebug()>1) {
      cout << "Writing the clusters... " << endl;
      for(unsigned int roi_i=0;roi_i<cluster_map.size();roi_i++) {
	cout << "ROI " << roi_i << ": " <<  (cluster_map[roi_i])->size() << " clusters" << endl;
      }
    }
    WriteClusters(evt, &cluster_map);
      
  } // nRoI>0

  if(seed_col->getNumberOfElements()>0)
    evt->addCollection(seed_col, "ClusterSeeds");

  cleanup(&cluster_map,preClusVec,trackVec);

  if(_garPars->GetDebug()>0) cout << "Finished Event " << evt->getEventNumber() << endl;
  
}


void ECALGarlic::cleanup(map<int,vector<ExtendedCluster* > *> *clusMap, vector<ExtendedCluster* > &preClusVec, vector<ExtendedTrack* > &trackVec)
{
  map<int,vector<ExtendedCluster* > *>::iterator clusMap_it;
  int NRois = clusMap->size();
  if(_garPars->GetDebug()>2) cout << "Deleting clusters in " << NRois << " ROIs" << endl;
  if(NRois>0) {
    for(clusMap_it=clusMap->begin();clusMap_it!=clusMap->end();clusMap_it++) {
      vector<ExtendedCluster* > *clusVec = (clusMap_it->second);
      int NClusters = clusVec->size();
      if(NClusters>0) {
	for (int i = 0; i < NClusters; i++) {
	  ExtendedCluster *a_ext_cluster = (*clusVec)[i];
	  Ellipsoid *sh = a_ext_cluster->Shape;
	  delete sh;
	  ClusterParameters *cp = a_ext_cluster->Parameters;
	  delete cp;
	  delete a_ext_cluster;
	  //vector<ExtendedHit* > ext_hit_vec = a_ext_cluster->hitVec;
	  /*
	    int NHitsInCluster = ext_hit_vec.size();
	    for(int hit_i=0;hit_i<NHitsInCluster;hit_i++) {
	    ExtendedHit *a_ext_hit = dynamic_cast<ExtendedHit*>(ext_hit_vec[hit_i]);
	    CalorimeterHit *a_hit = a_ext_hit->hit();
	    //delete a_hit;
	    delete a_ext_hit;
	    }
	  */
	}
      }
    }
  }
  int NPreClusters = preClusVec.size();
  if(_garPars->GetDebug()>2) cout << "Deleting " << NPreClusters << " PreClusters" << endl;
  if(NPreClusters>0) {
    for (int i = 0; i < NPreClusters; i++) {
      ExtendedCluster *a_ext_preCluster = preClusVec[i];
      vector<ExtendedHit* > ext_preHit_vec = a_ext_preCluster->hitVec;
      int NHitsInPreCluster = ext_preHit_vec.size();
      for(int pHit_i=0;pHit_i<NHitsInPreCluster;pHit_i++) {
	ExtendedHit *a_ext_pHit = dynamic_cast<ExtendedHit*>(ext_preHit_vec[pHit_i]);
	delete a_ext_pHit;
      }
      delete a_ext_preCluster;
    }
  }
  int NTracks = trackVec.size();
  if(_garPars->GetDebug()>2) cout << "Deleting " << NTracks << " Tracks" << endl;
  if(NTracks>0) {
    trackVec.erase(trackVec.begin(),trackVec.end());
  }
  if(_garPars->GetDebug()>2) cout << "Cleanup finished" << endl;
}



void ECALGarlic::check(LCEvent * evt)
{
}



void ECALGarlic::end()
{
  cout << endl;
  cout << "ECALGarlic Report: " << endl;
  cout << _nEvents << " events processed" << endl;
  cout << "Found " << _nClusters << " Cluster" << endl;
  if(_nEvents>0)
    cout << "= " << ((float)_nClusters)/_nEvents  << " /event" << endl;
}

bool ECALGarlic::RejectByMLPCut(ExtendedCluster *myCluster)
{
  ClusterParameters *par = myCluster->Parameters;
  double MLP = par->MLP;
  double cut=0;

  if(par->E_GeV > 20. && _geomPars->Get_nCellsPerWafer()==18) { // distance cut
    double dist = (par->smallestDistToTrack);
    if(dist<12.)
      return 0;
  }

  if(MLP==1)
    return 1;

  if(par->zone==1 || par->zone==3) {
    if(par->E_GeV <= 0.25)
      cut = _garPars->GetMLPCutsBarrel()[0];
    if(par->E_GeV > 0.25 && par->E_GeV <=0.35)
      cut = _garPars->GetMLPCutsBarrel()[1];
    if(par->E_GeV > 0.35 && par->E_GeV <=0.5)
      cut = _garPars->GetMLPCutsBarrel()[2];
    if(par->E_GeV > 0.5 && par->E_GeV <=0.75)
      cut = _garPars->GetMLPCutsBarrel()[3];
    if(par->E_GeV > 0.75 && par->E_GeV <=1.)
      cut = _garPars->GetMLPCutsBarrel()[4];
    if(par->E_GeV > 1. && par->E_GeV <=1.25)
      cut = _garPars->GetMLPCutsBarrel()[5];
    if(par->E_GeV > 1.25 && par->E_GeV <=1.5)
      cut = _garPars->GetMLPCutsBarrel()[6];
    if(par->E_GeV > 1.5 && par->E_GeV <=1.75)
      cut = _garPars->GetMLPCutsBarrel()[7];
    if(par->E_GeV > 1.75 && par->E_GeV <=2.25)
      cut = _garPars->GetMLPCutsBarrel()[8];
    if(par->E_GeV > 2.25 && par->E_GeV <=3.)
      cut = _garPars->GetMLPCutsBarrel()[9];
    if(_geomPars->Get_nCellsPerWafer()==18) { 
      if(par->E_GeV > 3. && par->E_GeV <=5.)
	cut = _garPars->GetMLPCutsBarrel()[10];
      if(par->E_GeV > 5. && par->E_GeV <=10.)
	cut = _garPars->GetMLPCutsBarrel()[11];
      if(par->E_GeV > 10. && par->E_GeV <=20.)
	cut = _garPars->GetMLPCutsBarrel()[12];
    }
    else {
      if(par->E_GeV > 3. && par->E_GeV <=10.)
	cut = _garPars->GetMLPCutsBarrel()[10];
    }
    //if(par->E_GeV > 20. && par->E_GeV <=50.)
    //cut = _garPars->GetMLPCutsBarrel()[12];
    if(par->E_GeV > 20.)
      return 1;
  }
  if(par->zone==2) {
    if(par->E_GeV <= 0.25)
      cut = _garPars->GetMLPCutsEndcap()[0];
    if(par->E_GeV > 0.25 && par->E_GeV <=0.35)
      cut = _garPars->GetMLPCutsEndcap()[1];
    if(par->E_GeV > 0.35 && par->E_GeV <=0.5)
      cut = _garPars->GetMLPCutsEndcap()[2];
    if(par->E_GeV > 0.5 && par->E_GeV <=0.75)
      cut = _garPars->GetMLPCutsEndcap()[3];
    if(par->E_GeV > 0.75 && par->E_GeV <=1.)
      cut = _garPars->GetMLPCutsEndcap()[4];
    if(par->E_GeV > 1. && par->E_GeV <=1.25)
      cut = _garPars->GetMLPCutsEndcap()[5];
    if(par->E_GeV > 1.25 && par->E_GeV <=1.5)
      cut = _garPars->GetMLPCutsEndcap()[6];
    if(par->E_GeV > 1.5 && par->E_GeV <=1.75)
      cut = _garPars->GetMLPCutsEndcap()[7];
    if(par->E_GeV > 1.75 && par->E_GeV <=2.25)
      cut = _garPars->GetMLPCutsEndcap()[8];
    if(par->E_GeV > 2.25 && par->E_GeV <=3.)
      cut = _garPars->GetMLPCutsEndcap()[9];
    if(_geomPars->Get_nCellsPerWafer()==18) { 
      if(par->E_GeV > 3. && par->E_GeV <=5.)
	cut = _garPars->GetMLPCutsEndcap()[10];
      if(par->E_GeV > 5. && par->E_GeV <=10.)
	cut = _garPars->GetMLPCutsEndcap()[11];
      if(par->E_GeV > 10. && par->E_GeV <=20.)
	cut = _garPars->GetMLPCutsEndcap()[12];
    }
    else {
      if(par->E_GeV > 3. && par->E_GeV <=10.)
	cut = _garPars->GetMLPCutsEndcap()[10];
    }
    //if(par->E_GeV > 20. && par->E_GeV <=50.)
    //cut = _garPars->GetMLPCutsEndcap()[12];
    if(par->E_GeV > 20.)
      return 1;
  }

  if(MLP>=cut)
    return 1;
  else
    return 0;
}




void ECALGarlic::PrepareTracks(const LCEvent *evt, vector<ExtendedTrack* > &trackVec)
{
  if(_garPars->GetDebug()>2) {
    cout << "Getting tracks..." << endl
         << "BField: " << _geomPars->Get_bField() << " T" << endl;
  }

  if(_garPars->GetCheatTracks()==1 ) { // take MC momentum

    if (_mcParticleColl->getNumberOfElements()>0) {
      MCParticle *a_primary = dynamic_cast<MCParticle*>(_mcParticleColl->getElementAt(0));
      float primary_dir[3];
      primary_dir[0] = a_primary->getMomentum()[0];
      primary_dir[1] = a_primary->getMomentum()[1];
      primary_dir[2] = a_primary->getMomentum()[2];
      float origin[3] = {0,0,0};
      float charge = a_primary->getCharge();

      if(_garPars->GetDebug()>2) 
	cout << "Primary particle " << endl
	     << "Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << endl
	     << "Momentum: " << primary_dir[0] << ", " << primary_dir[1] << ", " << primary_dir[2] << endl
	     << "Charge: "<< charge << endl;

      ExtendedTrack *a_ext_track = new ExtendedTrack();
      HelixClass *a_helix = new HelixClass();
      a_helix->Initialize_VP(origin,primary_dir,charge,_geomPars->Get_bField());
      a_ext_track->helix=a_helix;
      float R = a_helix->getRadius();
      float entryPoint[3]={0,0,0};
      float refPoint[3];
      float z0=a_helix->getZ0();
      refPoint[0] = (a_helix->getReferencePoint())[0];
      refPoint[1] = (a_helix->getReferencePoint())[1];
      refPoint[2] = (a_helix->getReferencePoint())[2];
      if(_garPars->GetDebug()>2) {
	cout << "RefPoint: " << refPoint[0] << " "  << refPoint[2]<< " " << refPoint[2] <<  endl;
	cout << "Momentum: " << a_helix->getMomentum()[0] << " "<< a_helix->getMomentum()[1] << " "<< a_helix->getMomentum()[2] << endl;
      }

      float time=0;
      float shortest_time=1.0e+10;
      if(R < _geomPars->Get_rOfBarrel()) { // Endcap
	if(_garPars->GetDebug()>2)
	  cout << " enters in Endcap " << endl;
	if(z0>=0)
	  time = a_helix->getPointInZ(_geomPars->Get_zOfEndcap(), refPoint, entryPoint);
	else
	  time = a_helix->getPointInZ(-_geomPars->Get_zOfEndcap(), refPoint, entryPoint);
      } else { // Barrel : for each stave calculate a possible entry point, then take the smallest one in time
	if(_garPars->GetDebug()>2) cout << " enters in Barrel " << endl;
	for(int stave_i=0;stave_i<_geomPars->Get_symmetry();stave_i++) {
	  float point[2];
	  point[0]=_geomPars->Get_rOfBarrel()*_geomPars->Get_barrelStaveDir()[stave_i].x;
	  point[1]=_geomPars->Get_rOfBarrel()*_geomPars->Get_barrelStaveDir()[stave_i].y;
	  float posEntry[3];
	  time = a_helix->getPointInXY(point[0],point[1],_geomPars->Get_barrelStaveDir()[stave_i].x,_geomPars->Get_barrelStaveDir()[stave_i].y,refPoint,posEntry);
	  if(time>0 && time<1.0e+10 && time<shortest_time) {
	    shortest_time=time;
	    entryPoint[0] = posEntry[0];
	    entryPoint[1] = posEntry[1];
	    entryPoint[2] = posEntry[2];
	  }
	}
      }
      vec3 ecalEntry;
      ecalEntry.x=entryPoint[0];
      ecalEntry.y=entryPoint[1];
      ecalEntry.z=entryPoint[2];
      if(_garPars->GetDebug()>2) cout << "Track enters ECAL at " << ecalEntry.x << ", "  << ecalEntry.y << ", " << ecalEntry.z << endl;
      a_ext_track->ecalEntryPoint=ecalEntry;
      trackVec.push_back(a_ext_track);
    }

  } else if (_garPars->GetCheatTracks()==0) {      // take tracker momentum

    LCCollection *trackColl = 0;
    try {
      trackColl = evt->getCollection(_LDCTrackCollectionName);
    }
    catch(DataNotAvailableException err) {};

    if(trackColl) {
      int NTracks=trackColl->getNumberOfElements();

      if(_garPars->GetDebug()>0) cout << NTracks << " tracks in track collection" << endl;

      for(int track_i=0;track_i<NTracks;track_i++) {

	if(_garPars->GetDebug()>1)cout << "Track " << track_i << endl;

	Track *a_track = dynamic_cast<Track*>(trackColl->getElementAt(track_i));
	ExtendedTrack *a_ext_track = new ExtendedTrack();

	float phi0=a_track->getPhi();
	float d0=a_track->getD0();
	float z0=a_track->getZ0();
	float omega=a_track->getOmega();
	float tanlambda=a_track->getTanLambda();

	if(_garPars->GetDebug()>2) cout << "Parameters: " << phi0 << " "  << d0 << " " << z0 << " " << omega << " " << tanlambda << endl;

	// create helix
	HelixClass *a_helix = new HelixClass();
	a_helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _geomPars->Get_bField());     
	a_ext_track->helix=a_helix;

	// decide if particle will enter ECAL in Endcap or Barrel FIXME: crude approach
	float R = a_helix->getRadius();
	float entryPoint[3]={0,0,0};
	float refPoint[3];
	refPoint[0] = (a_helix->getReferencePoint())[0];
	refPoint[1] = (a_helix->getReferencePoint())[1];
	refPoint[2] = (a_helix->getReferencePoint())[2];
	if(_garPars->GetDebug()>2) cout << "RefPoint: " << refPoint[0] << " "  << refPoint[2]<< " " << refPoint[2] <<  endl;

	if(_garPars->GetDebug()>2) cout << "Momentum: " << a_helix->getMomentum()[0] << " "<< a_helix->getMomentum()[1] << " "<< a_helix->getMomentum()[2] << endl;

	if(  abs(a_helix->getMomentum()[0]) < 0.0000000001 && abs(a_helix->getMomentum()[1]) < 0.0000000001 ) {
	  cout << "Strange Track!" << endl;
	}

	float time=0;
	float shortest_time=1.0e+10;

	if(R < _geomPars->Get_rOfBarrel()) { // Endcap
	  if(_garPars->GetDebug()>2) cout << " enters in Endcap " << endl;
	  if(z0>=0) time = a_helix->getPointInZ(_geomPars->Get_zOfEndcap(), refPoint, entryPoint);
	  else      time = a_helix->getPointInZ(-_geomPars->Get_zOfEndcap(), refPoint, entryPoint);
	} else { // Barrel : for each stave calculate a possible entry point, then take the smallest one in time
	  if(_garPars->GetDebug()>2) cout << " enters in Barrel " << endl;
	  for(int stave_i=0;stave_i<_geomPars->Get_symmetry();stave_i++) {
	    float point[2];
	    point[0]=_geomPars->Get_rOfBarrel()*_geomPars->Get_barrelStaveDir()[stave_i].x;
	    point[1]=_geomPars->Get_rOfBarrel()*_geomPars->Get_barrelStaveDir()[stave_i].y;
	    float posEntry[3];
	    time = a_helix->getPointInXY(point[0],point[1],_geomPars->Get_barrelStaveDir()[stave_i].x,_geomPars->Get_barrelStaveDir()[stave_i].y,refPoint,posEntry);
	    if(time>0 && time<1.0e+10 && time<shortest_time) {
	      shortest_time=time;
	      entryPoint[0] = posEntry[0];
	      entryPoint[1] = posEntry[1];
	      entryPoint[2] = posEntry[2];
	    }
	  }
	}
	vec3 ecalEntry;
	ecalEntry.x=entryPoint[0];
	ecalEntry.y=entryPoint[1];
	ecalEntry.z=entryPoint[2];
	if(_garPars->GetDebug()>2) cout << "Track " << track_i << " enters ECAL at " << ecalEntry.x << ", "  << ecalEntry.y << ", " << ecalEntry.z << endl;
	  
	a_ext_track->ecalEntryPoint=ecalEntry;
	trackVec.push_back(a_ext_track);
      }
    }
  }

  if (_garPars->GetDebug()>1) cout << "Having " << trackVec.size() << " tracks" << endl;
}



void ECALGarlic::PreparePreClusters(LCEvent *evt, vector<ExtendedCluster* > &preClusVec) 
{
  vector<ExtendedHit *> extPreShowerHits;
  if(_garPars->GetIncludePreshower()) {

    IMPL::LCCollectionVec *ECAL_ps_coll = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
    ECAL_ps_coll->setFlag(ECAL_ps_coll->getFlag()|( 1 << LCIO::RCHBIT_LONG));

    // Get and convert PreShower hits
    LCCollection *barrelPSColl = 0;
    try {
      barrelPSColl = evt->getCollection(_ecalBarrelPreShowerHitCollectionName);
    }
    
    catch(DataNotAvailableException err) {};
    if(barrelPSColl) {
      for(int element_i=0; element_i < (int) barrelPSColl->getNumberOfElements(); element_i++)  {
	SimCalorimeterHit *a_sim_hit=dynamic_cast<SimCalorimeterHit *>(barrelPSColl->getElementAt(element_i));
	float sim_en=a_sim_hit->getEnergy();
	if(sim_en>0.0000724) {
	  CalorimeterHitImpl *a_hit=new CalorimeterHitImpl;
	  a_hit->setCellID0(a_sim_hit->getCellID0());
	  a_hit->setEnergy(a_sim_hit->getEnergy());
	  a_hit->setPosition(a_sim_hit->getPosition());
	  ExtendedHit *a_ext_hit = new ExtendedHit();
	  a_ext_hit->hit = a_hit;
	  _geomHelpers->AssignPseudoLayer(a_ext_hit);
	  a_ext_hit->clusterHitVec=0;
	  a_ext_hit->preShower=1;
	  extPreShowerHits.push_back(a_ext_hit);
	  ECAL_ps_coll->addElement(a_hit);
	}
      }
    }
    LCCollection *endcapPSColl = 0;
    try {
      endcapPSColl = evt->getCollection(_ecalEndcapPreShowerHitCollectionName);
    }
    catch(DataNotAvailableException err) {};
    if(endcapPSColl) {
      for(int element_i=0; element_i < (int) endcapPSColl->getNumberOfElements(); element_i++)  {
	SimCalorimeterHit *a_sim_hit=dynamic_cast<SimCalorimeterHit *>(endcapPSColl->getElementAt(element_i));
	float sim_en=a_sim_hit->getEnergy();
	if(sim_en>0.0000724) {
	  CalorimeterHitImpl *a_hit=new CalorimeterHitImpl;
	  a_hit->setCellID0(a_sim_hit->getCellID0());
	  a_hit->setEnergy(a_sim_hit->getEnergy());
	  a_hit->setPosition(a_sim_hit->getPosition());
	  ExtendedHit *a_ext_hit = new ExtendedHit();
	  a_ext_hit->hit = a_hit;
	  _geomHelpers->AssignPseudoLayer(a_ext_hit);
	  a_ext_hit->clusterHitVec=0;
	  a_ext_hit->preShower=1;
	  extPreShowerHits.push_back(a_ext_hit);
	  ECAL_ps_coll->addElement(a_hit);
	}
      }
    }
    evt->addCollection(ECAL_ps_coll,"ECAL_PS");
  }
 

  LCCollection *preClusColl = 0;
  try {
    preClusColl = evt->getCollection(_garPars->GetEcalPreClusterCollectionName());
  }
  catch(DataNotAvailableException &exc) {
    if (_garPars->GetDebug()>1)
      std::cerr << "In ECALGarlic::preparePreClusters "	<< exc.what() << std::endl;
    return;
  }
  
  int NPreCluster = preClusColl->getNumberOfElements();
  
  if(!_geomPars->Get_defaultDecoder() && NPreCluster>0) {
    CellIDDecoder<CalorimeterHit>* dec = new CellIDDecoder<CalorimeterHit> (preClusColl);
    _geomPars->Set_defaultDecoder(dec);
    //    ECALGarlicGeometryParameters::_defaultDecoder = new CellIDDecoder<CalorimeterHit> (preClusColl);
  }

  if(_garPars->GetDebug()>2)
    cout << "Number of PreClusters in original collection: " << preClusColl->getNumberOfElements() << endl;
  for (int preC_i = 0; preC_i < NPreCluster; preC_i++) {
    ClusterImpl *a_cluster = dynamic_cast<ClusterImpl* >( preClusColl->getElementAt(preC_i) );
    bool is_in_barrel=0;
    bool is_in_endcap=0;
    ExtendedCluster *a_ext_cluster = new ExtendedCluster();
    const CalorimeterHitVec &preHitVec=a_cluster->getCalorimeterHits();
    //CellIDDecoder<CalorimeterHit> decoder(preClusColl);
    if (_garPars->GetDebug()>1) cout << "Reading PreCluster " << preC_i << " " << preHitVec.size() << endl;
    vector<ExtendedHit *> extHitVec;
    for (EVENT::CalorimeterHitVec::const_iterator hit_it=preHitVec.begin();
	 hit_it!=preHitVec.end();
	 hit_it++) {
      CalorimeterHitImpl *a_hit = dynamic_cast<CalorimeterHitImpl*>( *hit_it );
      //for histogramming
      ExtendedHit *a_ext_hit = new ExtendedHit();

      int module = (*(_geomPars->Get_defaultDecoder()))(a_hit)["M"];
      //int stave = decoder(a_hit)["S-1"];

      // daniel removes this...its not really a good idea i think
      // also requires digi to be done at same time as garlic...
      //      a_hit->setEnergy((a_hit->getEnergy()/44.02)*40.91); // go to old calibration

      if(0 < module && module < 6)
	is_in_barrel=1;
      else {
	is_in_endcap=1;
	//double old_energy = a_hit->getEnergy();

	// daniel removes this...its not really a good idea i think
	// also requires digi to be done at same time as garlic...
	// a_hit->setEnergy((a_hit->getEnergy()/102.5)*100);

	//cout << "Energy changed for Endcap hit: " << old_energy << " -> " << a_hit->getEnergy() << " or " << (*hit_it)->getEnergy() << endl;
      }
      a_ext_hit->hit = a_hit;
      _geomHelpers->AssignPseudoLayer(a_ext_hit);
      //a_ext_hit->clusterHitVec=&extHitVec;
      a_ext_hit->clusterHitVec=0;
      a_ext_hit->preShower=0;
      extHitVec.push_back(a_ext_hit);
      //cout << "Layer: " << a_ext_hit->pseudoLayer << ", Hit " << a_ext_hit->hit->getCellID0() << endl ;

      if(extPreShowerHits.size()>0) {   // include PreShower hits
	vec3 hitPos;
	hitPos.x=a_hit->getPosition()[0];
	hitPos.y=a_hit->getPosition()[1];
	hitPos.z=a_hit->getPosition()[2];
	for(vector<ExtendedHit*>::iterator psh_i=extPreShowerHits.begin();psh_i!=extPreShowerHits.end();) {
	  ExtendedHit *a_ext_hit = dynamic_cast<ExtendedHit*>(*psh_i);
	  vec3 psHitPos;
	  psHitPos.x=a_ext_hit->hit->getPosition()[0];
	  psHitPos.y=a_ext_hit->hit->getPosition()[1];
	  psHitPos.z=a_ext_hit->hit->getPosition()[2];
	  double dist = _geomHelpers->Get3dDistance(&hitPos,&psHitPos);
	  if(dist<300) { // include in PreCluster
	    extHitVec.push_back(a_ext_hit);
	    extPreShowerHits.erase(psh_i);
	  }
	  else
	    psh_i++;
	}
      }
      
    }
    // sort hits by pseudo layer 
    sort(extHitVec.begin(),extHitVec.end(),ExtendedHit::lowerPseudoLayer);
    
    if(_garPars->GetDebug()>2) { //excessive debugging of hit in psuedo layers
      int nHits=extHitVec.size();
      if(_garPars->GetDebug()>2)
	cout << "PreCluster " << preC_i << " : " << nHits << " hits" << endl;
      ExtendedHit *a_ext_hit;
      for(int i=0;i<nHits;i++) {
	a_ext_hit=extHitVec[i];
	vec3 hitPos;
	hitPos.x=a_ext_hit->hit->getPosition()[0];
	hitPos.y=a_ext_hit->hit->getPosition()[1];
	hitPos.z=a_ext_hit->hit->getPosition()[2];
      }
    }
    a_ext_cluster->hitVec=extHitVec;
    if(is_in_barrel==1 && is_in_endcap==1)
      a_ext_cluster->location=CLUS_LOCATION_OVERLAP;
    else 
      if(is_in_barrel==1 && is_in_endcap==0)
	a_ext_cluster->location=CLUS_LOCATION_BARREL;
      else
	a_ext_cluster->location=CLUS_LOCATION_ENDCAP;
    preClusVec.push_back(a_ext_cluster);
  }
}


void ECALGarlic::RemoveHitsNearExtrapolatedTracks(LCEvent *evt,
						  vector<ExtendedTrack* > &trackVec, vector<ExtendedCluster* > &preClusVec, 
						  map<ExtendedTrack*,map<int,vector<ExtendedHit*> > > *allRemovedHits)
{
  IMPL::LCCollectionVec *removedHits = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
  removedHits->setFlag(removedHits->getFlag()|( 1 << LCIO::RCHBIT_LONG));
  int NTracks = trackVec.size();
  int NPreClusters = preClusVec.size();
  for(int track_i=0;track_i<NTracks;track_i++) {
    if(_garPars->GetDebug()>2) 
      cout << " Checking track " << track_i << " :" << endl;
    ExtendedTrack *a_ext_track = dynamic_cast<ExtendedTrack* >(trackVec[track_i]);
    /*
    Track *a_track = dynamic_cast<Track*>(a_ext_track->track);
    float phi0=a_track->getPhi();
    float d0=a_track->getD0();
    float z0=a_track->getZ0();
    float omega=a_track->getOmega();
    float tanlambda=a_track->getTanLambda();
    // create helix
    HelixClass *a_helix = new HelixClass();
    a_helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _geomPars->Get_bField());
    */
    map<int,vector<ExtendedHit*> > removedHere;
    vector<ExtendedHit*> hitsPerLayer[29];  // dismiss preShower hits
    HelixClass *a_helix = a_ext_track->helix;
    for(int pclus_i=0;pclus_i<NPreClusters;pclus_i++) {
      if(_garPars->GetDebug()>2) 
	cout << "Removing hits from PreCluster " << pclus_i << endl;
      ExtendedCluster *a_pre_cluster = preClusVec[pclus_i];
      vector<ExtendedHit* > *hit_vec = &(a_pre_cluster->hitVec);
      vector<ExtendedHit* >::iterator hit_it;
      for(hit_it=hit_vec->begin();hit_it!=hit_vec->end();) {
	ExtendedHit *a_ext_hit = dynamic_cast<ExtendedHit* >(*hit_it);
	float hitPos[3];
	hitPos[0]=a_ext_hit->hit->getPosition()[0];
	hitPos[1]=a_ext_hit->hit->getPosition()[1];
	hitPos[2]=a_ext_hit->hit->getPosition()[2];
	int hit_layer=(*(_geomPars->Get_defaultDecoder()))(a_ext_hit->hit)["K-1"];
	//int psLayer = a_ext_hit->pseudoLayer;
	float dist[3];
	float dummy;
	dummy = a_helix->getDistanceToPoint(hitPos,dist);
	double dist_cut=1;
	if(_geomPars->Get_padSizeEcal()[1]>6)
	  dist_cut = _geomPars->Get_padSizeEcal()[1];
	else
	  dist_cut = 2.*_geomPars->Get_padSizeEcal()[1];
	if(dist[2]<dist_cut) {       
	  hit_it = hit_vec->erase(hit_it);
	  CalorimeterHitImpl *a_hit = new CalorimeterHitImpl;
	  CalorimeterHitImpl *b_hit  = dynamic_cast<CalorimeterHitImpl* >(a_ext_hit->hit);
	  (*a_hit) = (*b_hit);	  
	  removedHits->addElement(a_hit);
	  if(a_ext_hit->preShower==0 && hit_layer<29) // dismiss preShower hits, and exclude last layer in ECAL plug for convenience
	    (hitsPerLayer[hit_layer]).push_back(a_ext_hit);
	  if(_garPars->GetDebug()>2) 
	    cout << "Removed hit "<< a_ext_hit->hit->getCellID0() << " at "  << hitPos[0] << ", " << hitPos[1] << ", " << hitPos[2] << ", with dist = " << dist[2] << " in Layer " << hit_layer << endl;
	}
	else {
	  if(_garPars->GetDebug()>2)
	    cout << "Hit not removed: " << a_ext_hit->hit->getCellID0() << ", with dist = " << dist[2] << " while asked: "<< dist_cut <<  endl;
	  hit_it++;
	}
      }
      if(_garPars->GetDebug()>2)
	cout << "Done with PreCluster " << pclus_i << endl;
    }
    for(int i=0;i<29;i++) {
      if(_garPars->GetDebug()>2)
	cout << "Layer " << i << ": " << hitsPerLayer[i].size() << " hits removed" << endl;
      //if(hitsPerLayer[i].size()>0)
      removedHere.insert(make_pair(i,hitsPerLayer[i]));
    }
    allRemovedHits->insert(make_pair(a_ext_track,removedHere));
  }
  evt->addCollection(removedHits,_removedHitsCollectionName);
}


void ECALGarlic::ProjectForSeeding(ExtendedCluster &preClus, vector<vec3> &possibleSeeds, LCCollectionVec *seed_col)
{
  
  //TFile *rootfile=TFile::Open("ProjectionHistos.root","RECREATE");
  //rootfile->cd();
  double hitSizeInPhi=fabs(atan2((double)_geomPars->Get_padSizeEcal()[1], (double)_geomPars->Get_rOfBarrel()));
  double hitSizeInTheta=fabs(atan2((double)_geomPars->Get_padSizeEcal()[1], (double)_geomPars->Get_zOfBarrel()));
  if(_garPars->GetDebug()>2) {
    cout << "Hitsize in phi is " << hitSizeInPhi << endl;
    cout << "Hitsize in theta is " << hitSizeInTheta << endl;
  }
  vector<ExtendedHit* > consideredHits;
  vector<ExtendedHit* > *allHits=&(preClus.hitVec);
  int NHits=allHits->size();
  
  multimap<int,ExtendedHit* > binHitMap;
  multimap<int,ExtendedHit* >::iterator bHM_it;
  
  switch(preClus.location)
    {
    case CLUS_LOCATION_BARREL:   // hits are projected on a cylindrical surface in z-phi with r=innerRadius of ECAL
      {
	if(_garPars->GetDebug()>2)
	  cout << "Projecting Barrel Hits" << endl;
	for(int hit_i=0;hit_i<NHits;hit_i++) {
	  ExtendedHit *a_ext_hit=(*allHits)[hit_i];
	  if(a_ext_hit->pseudoLayer<_garPars->GetNLayersForSeeding()+2){
	    consideredHits.push_back(a_ext_hit);
	  }
	}
	bool phiSeam=false;
	int NConsideredHits=consideredHits.size();
	if(NConsideredHits<_garPars->GetMinHitsForSeed()) {
	  if(_garPars->GetDebug()>2)
	    cout << "not enough hits" << endl;
	  break;
	}
	//get dimension of projection area
	float phiMin=2*pi;
	float phiMax=-2*pi;

	float zMin=_geomPars->Get_zOfBarrel();
	float zMax=-_geomPars->Get_zOfBarrel();
	
	//float max_en=0;
	//ExtendedHit *max_en_hit;
	for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	  
	  ExtendedHit *a_ext_hit=consideredHits[hit_i];
	  //project  
	  float x = (a_ext_hit->hit)->getPosition()[0];
	  float y = (a_ext_hit->hit)->getPosition()[1];
	  float z = (a_ext_hit->hit)->getPosition()[2];
	  //float r = sqrt(x*x+y*y+z*z);
	  double z_proj = z*_geomPars->Get_rOfBarrel()/sqrt(x*x+y*y);
	  double phi = atan2((double)x, (double)y);
	  //  float hit_en = (a_ext_hit->hit)->getEnergy();
	  //if(hit_en>max_en)
	  //  ma_en_hit = a_ext_hit;
	  if(phi<0)
	    phi=twopi+phi;

	  if (z_proj>zMax)
	    zMax=z_proj;
	  if (z_proj<zMin)
	    zMin=z_proj;
	  if (phi>phiMax)
	    phiMax=phi;
	  if (phi<phiMin)
	    phiMin=phi;
	}

	zMin=zMin-(_geomPars->Get_padSizeEcal()[0]/2);
	zMax=zMax+(_geomPars->Get_padSizeEcal()[0]/2);
	double zMean = (zMin+zMax)/2;
	double sinTh = abs(_geomPars->Get_rOfBarrel()/(sqrt(_geomPars->Get_rOfBarrel()*_geomPars->Get_rOfBarrel()+zMean*zMean)));
	double newBinSizeZ = _geomPars->Get_padSizeEcal()[0]/sinTh;
	phiMin=phiMin-(hitSizeInPhi/2);
	phiMax=phiMax+(hitSizeInPhi/2);
	int zBins=(int)(fabs((zMax-zMin))/newBinSizeZ);
	int phiBins=(int)(fabs((phiMax-phiMin))/hitSizeInPhi);
	if(phiMax>6.2 && phiMin<0.1) { // spread over seam in phi, need to build new histogram with seamless phi values
	  phiSeam=true;
	  if(_garPars->GetDebug()>2)
	    cout << "Building histogram without seam in PHI..." << endl;
	  phiMin=2*pi;
	  phiMax=-2*pi;
	  for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	    ExtendedHit *a_ext_hit=consideredHits[hit_i];
	    //project  
	    float x = (a_ext_hit->hit)->getPosition()[0];
	    float y = (a_ext_hit->hit)->getPosition()[1];
	    //float z = (a_ext_hit->hit)->getPosition()[2];
	    //float r = sqrt(x*x+y*y+z*z);
	    //float z_proj = z*_geomPars->Get_rOfBarrel()/sqrt(x*x+y*y);
	    double phi = atan2((double)x, (double)y);
	    
	    if(phi<0)
	      phi=twopi+phi;

	    if(phi>=pi)
	      phi=phi-pi;
	    else
	      phi=phi+pi;
	    //if (z_proj>zMax)
	    //  zMax=z_proj;
	    //if (z_proj<zMin)
	    //  zMin=z_proj;
	    if (phi>phiMax)
	      phiMax=phi;
	    if (phi<phiMin)
	      phiMin=phi;
	  }
	  phiMin=phiMin-(hitSizeInPhi/2);
	  phiMax=phiMax+(hitSizeInPhi/2);
	  phiBins=(int)(fabs((phiMax-phiMin))/hitSizeInPhi);
	}
	if(_garPars->GetDebug()>2) {
	  cout << "Histogram parameters: " <<  phiBins << ", " << phiMin<< ", " <<phiMax<< ", " <<zBins<< ", " <<zMin<< ", " <<zMax << endl;
	}
	if(_garPars->GetDebug()>2 && phiSeam==true)
	  cout << "..having a seam in PHI" << endl;
	TH2F *en_proj_histo = new TH2F("en_proj_histo","en_proj_histo",phiBins,phiMin,phiMax,zBins,zMin,zMax);
	TH2I *hit_proj_histo = new TH2I("hit_proj_histo","hit_proj_histo",phiBins,phiMin,phiMax,zBins,zMin,zMax);

	for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	  ExtendedHit *a_ext_hit=consideredHits[hit_i];
	  //project  
	  float x = (a_ext_hit->hit)->getPosition()[0];
	  float y = (a_ext_hit->hit)->getPosition()[1];
	  float z = (a_ext_hit->hit)->getPosition()[2];
	  double z_proj = z*_geomPars->Get_rOfBarrel()/sqrt((x*x+y*y));
	  double phi = atan2((double)x, (double)y);
	  if(phi<0)
	    phi=twopi+phi;
	  if(phiSeam==true) {
	    if(phi>=pi)
	      phi=phi-pi;
	    else
	      phi=phi+pi;
	  }
	  float en = (a_ext_hit->hit)->getEnergy();
	  en_proj_histo->Fill(phi,z_proj,en);
	  hit_proj_histo->Fill(phi,z_proj);
	  int bin = en_proj_histo->FindBin(phi,z_proj);
	  binHitMap.insert(make_pair(bin,a_ext_hit));
	  if(_garPars->GetDebug()>2) 
	    cout << "Histogrammed hit ; " <<  phi << ", " << z_proj << ", with E=" << en << endl;
	} 
	
	// order bins per energy deposit
	multimap<float,int>  proj_bins;
	proj_bins.clear();
	multimap<float,int>::reverse_iterator proj_bins_it;
	int nX  = en_proj_histo->GetXaxis()->GetNbins()+2;
	int NBins = en_proj_histo->GetBin(en_proj_histo->GetNbinsX(),en_proj_histo->GetNbinsY());
		  
	for(int bin_i=1;bin_i <= NBins;bin_i++) {
	  if (en_proj_histo->GetBinContent(bin_i) != 0) 
	    proj_bins.insert(make_pair(en_proj_histo->GetBinContent(bin_i),bin_i ));
	}
	
	vector<vector<int> > seed_clusters;
	vector<vector<int> >::iterator seed_clusters_it;
	vector<int>::iterator cl_it;
	
	for(proj_bins_it=proj_bins.rbegin();proj_bins_it!=proj_bins.rend();proj_bins_it++) {
	  int bin_i = (proj_bins_it->second);
	  bool in_seed=false;
	  for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	    for(cl_it=seed_clusters_it->begin();cl_it!=seed_clusters_it->end();cl_it++) {
	      if(bin_i == *cl_it) {
		in_seed=true;
		break;
	      }
	    }
	    if(in_seed==true)
	      break;
	    int old_bin_i = (seed_clusters_it->begin())[0];
	    int old_bin_x = old_bin_i%nX;
	    int old_bin_y = old_bin_i/nX;
	    int bin_x = bin_i%nX;
	    int bin_y = bin_i/nX;
	    if( (fabs(double(bin_x-old_bin_x))<3) && (fabs(double(bin_y-old_bin_y))<3) ) {
	      in_seed=true;
	      break;
	    }
	  }
	  if(in_seed==true)
	    continue;
	  else {
	    vector<int> a_seed;
	    float maxHitEn = 0;
	    double seedEn=0;
	    a_seed.push_back(bin_i);
	    //RecAddNeighbourBins(en_proj_histo, bin_i, &a_seed);
	    _clusterHelpers->IterAddNeighbourBins(en_proj_histo, bin_i, &a_seed,0);
	    // check if seed originated from more than 4 hit
	    int nOriginHits=0;
	    for(int a_bin=0; a_bin<((int)a_seed.size());a_bin++) {
	      int bin_j = a_seed[a_bin];
	      nOriginHits+=(int)(hit_proj_histo->GetBinContent(bin_j));
	      for(bHM_it=binHitMap.begin();bHM_it!=binHitMap.end();bHM_it++) {
		if(bHM_it->first!=bin_j)
		  continue;
		else {
		  float hit_en = (bHM_it->second)->hit->getEnergy();
		  if(hit_en>maxHitEn)
		    maxHitEn=hit_en;
		  seedEn+=hit_en;
		}
	      }
	    }
	    if(nOriginHits>3) {
	      seed_clusters.push_back(a_seed);
	      //if(seed_clusters.size()==1)
	      //	_seedEn_histo->Fill(seedEn);
	      //_seedHitEn_histo->Fill(maxHitEn);
	      //_seedEn_histo->Fill(seedEn);
	      if(_garPars->GetDebug()>2) 
		cout << "Possible seed found from  " << a_seed.size() << " bins" << endl;
	    }
	    else
	      if(_garPars->GetDebug()>2)
		cout << "Rejected seed from " << nOriginHits << " hits" << endl;
	  }
	}

	for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	  int bin_i = (seed_clusters_it->begin())[0];
	  double binPos[2];  // binPos[0]=phi, binPos[1]=z
	  binPos[0]=en_proj_histo->GetXaxis()->GetBinCenter(bin_i%nX);
	  binPos[1]=en_proj_histo->GetYaxis()->GetBinCenter(bin_i/nX);
	  // to get the seed position on the ECAL front plane we need to project back
	  double phi = binPos[0];
	  if(phiSeam==true) { // reverse the effect of the new Phi definition
	    if(phi>=0)
	      phi=phi+pi;
	    else
	      phi=phi-pi;
	  }
	  if(phi<0) {
	    cout << "Shouldnt be here: phi from hitso is negative!" << endl;
	    phi=twopi+phi;
	  }
	  int gamma=(int)(((phi-(twopi/(2*_geomPars->Get_symmetry())))/(twopi/_geomPars->Get_symmetry()))+1);
	  double GAMMA=gamma*(twopi/_geomPars->Get_symmetry());
	  double PHI = phi-gamma*(twopi/_geomPars->Get_symmetry());
	  double X=_geomPars->Get_rOfBarrel()*tan(PHI);
	  double Y=_geomPars->Get_rOfBarrel();
	  vec3 seedPos;	    
	  seedPos.x = X*cos(GAMMA)+Y*sin(GAMMA);   
	  seedPos.y = Y*cos(GAMMA)-X*sin(GAMMA);  
	  double z_proj = binPos[1];
	  seedPos.z = z_proj*sqrt(seedPos.x*seedPos.x+seedPos.y*seedPos.y)/_geomPars->Get_rOfBarrel();
	  
	  // automatically ordered by energy deposit!
	  possibleSeeds.push_back(seedPos);
	  CalorimeterHitImpl *a_seed=new CalorimeterHitImpl;
	  //a_hit->setCellID0(a_sim_hit->getCellID0());
	  float s_pos[3] = {seedPos.x,seedPos.y,seedPos.z};
	  a_seed->setEnergy(0);
	  a_seed->setPosition(s_pos);
	  seed_col->addElement(a_seed);
	}
	
	//en_proj_histo->Write();
	delete en_proj_histo;
	delete hit_proj_histo;
	break;
      }
      case CLUS_LOCATION_ENDCAP:   // hits are projected on the XY front plane of the Endcap
      {
	if(_garPars->GetDebug()>2)
	  cout << "Projecting Endcap Hits" << endl;
	for(int hit_i=0;hit_i<NHits;hit_i++) {
	  ExtendedHit *a_ext_hit=(*allHits)[hit_i];
	  if(a_ext_hit->pseudoLayer<_garPars->GetNLayersForSeeding()+2){
	    consideredHits.push_back(a_ext_hit);
	  }
	}
	int NConsideredHits=consideredHits.size();
	if(_garPars->GetDebug()>2)
	  cout << "Considering " << NConsideredHits << " hits" << endl;
	if(NConsideredHits<_garPars->GetMinHitsForSeed()) {
	if(_garPars->GetDebug()>2)
	  cout << "not enough hits" << endl;
	  break;
	}
	//to destinguish between the two endcaps
	bool zIsPos=true;
	//get dimension of projection area
	double xMin=2*_geomPars->Get_rOfEndcap();
	double xMax=-2*_geomPars->Get_rOfEndcap();
	double yMin=2*_geomPars->Get_rOfEndcap();
	double yMax=-2*_geomPars->Get_rOfEndcap();
	
	for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	  
	  ExtendedHit *a_ext_hit=consideredHits[hit_i];
	  //project on Endcap frontplane
	  float x_h = (a_ext_hit->hit)->getPosition()[0];
	  float y_h = (a_ext_hit->hit)->getPosition()[1];
	  float z_h = (a_ext_hit->hit)->getPosition()[2];
	  double scale = fabs(_geomPars->Get_zOfEndcap()/z_h);
	  
	  double x = x_h*scale;
	  double y = y_h*scale;

	  if (x>xMax)
	    xMax=x;
	  if (x<xMin)
	    xMin=x;
	  if (y>yMax)
	    yMax=y;
	  if (y<yMin)
	    yMin=y;
	}
	xMin=xMin-(_geomPars->Get_padSizeEcal()[0]/2);
	xMax=xMax+(_geomPars->Get_padSizeEcal()[0]/2);
	yMin=yMin-(_geomPars->Get_padSizeEcal()[0]/2);
	yMax=yMax+(_geomPars->Get_padSizeEcal()[0]/2);
	double xMean=(xMin+xMax)/2;
	double yMean=(yMin+yMax)/2;
	double cosThX = abs(_geomPars->Get_zOfEndcap()/(sqrt(_geomPars->Get_zOfEndcap()*_geomPars->Get_zOfEndcap()+xMean*xMean)));
	double cosThY = abs(_geomPars->Get_zOfEndcap()/(sqrt(_geomPars->Get_zOfEndcap()*_geomPars->Get_zOfEndcap()+yMean*yMean)));
	double newBinSizeX = _geomPars->Get_padSizeEcal()[0]/cosThX;
	double newBinSizeY = _geomPars->Get_padSizeEcal()[0]/cosThY;
	int xBins=(int)((xMax-xMin)/newBinSizeX);
	int yBins=(int)((yMax-yMin)/newBinSizeY);

	if(_garPars->GetDebug()>2) 
	  cout << "Histogram parameters: " <<  xBins << ", " << xMin<< ", " << xMax<< ", " <<yBins<< ", " <<yMin<< ", " <<yMax << endl;
	TH2F *en_proj_histo = new TH2F("en_proj_histo","en_proj_histo",xBins,xMin,xMax,yBins,yMin,yMax);
	TH2I *hit_proj_histo = new TH2I("hit_proj_histo","hit_proj_histo",xBins,xMin,xMax,yBins,yMin,yMax);

	for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	  
	  ExtendedHit *a_ext_hit=consideredHits[hit_i];
	  //project  
	  float x_h = (a_ext_hit->hit)->getPosition()[0];
	  float y_h = (a_ext_hit->hit)->getPosition()[1];
	  float z_h = (a_ext_hit->hit)->getPosition()[2];
	  double scale = fabs(_geomPars->Get_zOfEndcap()/z_h);

	  double x = x_h*scale;
	  double y = y_h*scale;

	  if((a_ext_hit->hit)->getPosition()[2] > 0)
	    zIsPos=true;
	  else
	    zIsPos=false;
	  float en = (a_ext_hit->hit)->getEnergy();
	  en_proj_histo->Fill(x,y,en);
	  hit_proj_histo->Fill(x,y);
	  int bin = en_proj_histo->FindBin(x,y);
	  binHitMap.insert(make_pair(bin,a_ext_hit));
	  if(_garPars->GetDebug()>2) 
	    cout << "Histogrammed hit ; " <<  x_h << ", " << y_h << ", with E=" << en << endl;
	}

	// order bins per energy deposit
	multimap<float,int>  proj_bins;
	proj_bins.clear();
	multimap<float,int>::reverse_iterator proj_bins_it;
	int nX  = en_proj_histo->GetXaxis()->GetNbins()+2;
	int NBins = en_proj_histo->GetBin(en_proj_histo->GetNbinsX(),en_proj_histo->GetNbinsY());

	for(int bin_i=1;bin_i <= NBins;bin_i++) {
	  if (en_proj_histo->GetBinContent(bin_i) != 0) 
	    proj_bins.insert(make_pair(en_proj_histo->GetBinContent(bin_i),bin_i ));
	}
	

	vector<vector<int> > seed_clusters;
	vector<vector<int> >::iterator seed_clusters_it;
	vector<int>::iterator cl_it;

	for(proj_bins_it=proj_bins.rbegin();proj_bins_it!=proj_bins.rend();proj_bins_it++) {
	  int bin_i = (proj_bins_it->second);
	  bool in_seed=false;
	  for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	    for(cl_it=seed_clusters_it->begin();cl_it!=seed_clusters_it->end();cl_it++) {
	      if(bin_i == *cl_it) {
		in_seed=true;
		break;
	      }
	    }
	    if(in_seed==true)
	      break;
	    int old_bin_i = (seed_clusters_it->begin())[0];
	    int old_bin_x = old_bin_i%nX;
	    int old_bin_y = old_bin_i/nX;
	    int bin_x = bin_i%nX;
	    int bin_y = bin_i/nX;
	    if( (fabs(double(bin_x-old_bin_x))<3) && (fabs(double(bin_y-old_bin_y))<3) ) {
	      in_seed=true;
	      break;
	    }
	  }

	  if(in_seed==true)
	    continue;
	  else {
	    vector<int> a_seed;
	    double maxHitEn = 0;
	    double seedEn = 0;
	    a_seed.push_back(bin_i);
	    //RecAddNeighbourBins(en_proj_histo, bin_i, &a_seed);
	    _clusterHelpers->IterAddNeighbourBins(en_proj_histo, bin_i, &a_seed,0);
	    // check if seed originated from more than 4 hits
	    int nOriginHits=0;
	    for(int a_bin=0; a_bin<((int)a_seed.size());a_bin++) {
	      int bin_j = a_seed[a_bin];
	      nOriginHits+=(int)(hit_proj_histo->GetBinContent(bin_j));
	      for(bHM_it=binHitMap.begin();bHM_it!=binHitMap.end();bHM_it++) {
		if(bHM_it->first!=bin_j)
		  continue;
		else {
		  float hit_en = (bHM_it->second)->hit->getEnergy();
		  if(hit_en>maxHitEn)
		    maxHitEn=hit_en;
		  seedEn+=hit_en;
		}
	      }
	    }

	    if(nOriginHits>3) {
	      seed_clusters.push_back(a_seed);
	      //_seedHitEn_histo->Fill(maxHitEn);
	      //if(seed_clusters.size()==1)
	      //		_seedEn_histo->Fill(seedEn);
	      if(_garPars->GetDebug()>2) 
		cout << "Possible seed found from  " << a_seed.size() << " bins" << endl;
	    }
	    else
	      if(_garPars->GetDebug()>2) 
		cout << "Rejected seed from " << nOriginHits << " hits" << endl;
	  }
	}

	for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	  int bin_i = (seed_clusters_it->begin())[0];
	  double binPos[2];  // binPos[0]=phi, binPos[1]=z
	  binPos[0]=en_proj_histo->GetXaxis()->GetBinCenter((bin_i)%nX);
	  binPos[1]=en_proj_histo->GetYaxis()->GetBinCenter((bin_i)/nX);
	  // to get the seed posittion on the ECAL front plane we need to project back
	  vec3 seedPos;
	  seedPos.x = binPos[0];
	  seedPos.y = binPos[1];
	  if (zIsPos==true)
	    seedPos.z = _geomPars->Get_zOfEndcap();
	  else
	    seedPos.z =-(_geomPars->Get_zOfEndcap());
	    
	  // automatically ordered by energy deposit!
	  possibleSeeds.push_back(seedPos);
	  CalorimeterHitImpl *a_seed=new CalorimeterHitImpl;
	  //a_hit->setCellID0(a_sim_hit->getCellID0());
	  float s_pos[3] = {seedPos.x,seedPos.y,seedPos.z};
	  a_seed->setEnergy(0);
	  a_seed->setPosition(s_pos);
	  seed_col->addElement(a_seed);
	}

	delete hit_proj_histo;
	delete en_proj_histo;
	break;
      }
    case CLUS_LOCATION_OVERLAP:  { // determine X0 in front of each hit!
      if(_garPars->GetDebug()>2) {
	cout << "Projecting Barrel-Endcap Overlap Region" << endl;
	cout << "Cos of Barrel is " << _geomPars->Get_cosOfBarrel() << endl;
      }
      
      //padsizes....
      //double projectedSizeTheta = _geomPars->Get_cosOfBarrel()*_geomPars->Get_padSizeEcal()[1];
      //double projectedSizePhi = cos(0.392699082)*_geomPars->Get_padSizeEcal()[1]; // cos(22.5deg)
      double effectiveRadius = sqrt(_geomPars->Get_zOfBarrel()*_geomPars->Get_zOfBarrel()+_geomPars->Get_rOfBarrel()*_geomPars->Get_rOfBarrel());
      //double newSizeTheta = fabs(atan2(projectedSizeTheta, effectiveRadius));
      double newSizeTheta = fabs(atan2((double)_geomPars->Get_padSizeEcal()[1], effectiveRadius));
      //double newSizePhi = fabs(atan2(projectedSizePhi, (double)_geomPars->Get_rOfBarrel()));
      double newSizePhi = fabs(atan2((double)_geomPars->Get_padSizeEcal()[1], (double)_geomPars->Get_rOfBarrel()));

      // estimate padsize in theta
      hitSizeInTheta = newSizeTheta; 
      // ...and phi
      hitSizeInPhi = newSizePhi;     // is cos 22.5 degrees at worst
      
      if(_garPars->GetDebug()>2) {
	cout << "New hitsize in phi is " << hitSizeInPhi << endl;
	cout << "New hitsize in theta is " << hitSizeInTheta << endl;
      }

      bool phiSeam=false;
      bool zIsPos=true;
      for(int hit_i=0;hit_i<NHits;hit_i++) {
	ExtendedHit *a_ext_hit=(*allHits)[hit_i];
	vec3 hitPos;
	hitPos.x = (a_ext_hit->hit)->getPosition()[0];
	hitPos.y = (a_ext_hit->hit)->getPosition()[1];
	hitPos.z = (a_ext_hit->hit)->getPosition()[2];
	if(hitPos.z<_geomPars->Get_zOfBarrel() && hitPos.z>-_geomPars->Get_zOfBarrel()) { //Barrel hit
	  if(a_ext_hit->pseudoLayer<_garPars->GetNLayersForSeeding()+2){
	    consideredHits.push_back(a_ext_hit);
	  }
	}
	else { //Endcap hit
	  double phi_pos = atan2((double)hitPos.x,(double)hitPos.y);
	  int gamma;
	  if(phi_pos<0) 
	    phi_pos=twopi+phi_pos;
	  gamma=(int)(((phi_pos-(twopi/(2*_geomPars->Get_symmetry())))/(twopi/_geomPars->Get_symmetry()))+1);
	  double GAMMA=gamma*(twopi/_geomPars->Get_symmetry());
	  double xPos=hitPos.x*cos(GAMMA)-hitPos.y*sin(GAMMA);
	  double yPos=hitPos.y*cos(GAMMA)+hitPos.x*sin(GAMMA);
	  double zPos=hitPos.z;
	  double cosHit=zPos/sqrt(yPos*yPos+zPos*zPos);
	  bool barrel_proj = 1;
	  //if(!(cosHit<_geomPars->Get_cosOfBarrel() && cosHit>-_geomPars->Get_cosOfBarrel())) { // Endcap hit without Barrel projection

	  if((abs(xPos)*_geomPars->Get_zOfBarrel()/zPos) < (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry())))) {
	    if((yPos*_geomPars->Get_zOfBarrel()/zPos) < _geomPars->Get_rOfBarrel())
	      barrel_proj=0;
	  }
	  else {
	    double y_offset = ((abs(xPos) - (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry())))) * tan(twopi/_geomPars->Get_symmetry()) );
	    if((yPos*_geomPars->Get_zOfBarrel()/zPos) < (_geomPars->Get_rOfBarrel()-y_offset))
	      barrel_proj=0;
	  }

	  if(barrel_proj==0) { // Endcap hit without Barrel projection
	    if(a_ext_hit->pseudoLayer<_garPars->GetNLayersForSeeding()+2){
	      if(_garPars->GetDebug()>2)
		cout << "Considering hit " << a_ext_hit->hit->getCellID0() << " without barrel projection" << endl;
	      consideredHits.push_back(a_ext_hit);
	    }
	    else
	      if(_garPars->GetDebug()>2)
		cout << "Rejected hit " << a_ext_hit->hit->getCellID0() << " without barrel projection" << endl;
	  }
	  else { // Endcap hit with Barrel projection: determine "passed material"!
	    if(a_ext_hit->pseudoLayer>_garPars->GetNLayersForSeeding()+1){
	      if(_garPars->GetDebug()>2)
		cout << "Rejected hit " << a_ext_hit->hit->getCellID0() << " without barrel projection" << endl;
	      continue;
	    }
	    vec3 entryBarrel;
	    entryBarrel.y=_geomPars->Get_rOfBarrel();
	    entryBarrel.x=xPos*entryBarrel.y/yPos;
	    entryBarrel.z=zPos*entryBarrel.y/yPos;
	    vec3 exitBarrel;
	    if(zPos>0)
	      exitBarrel.z=_geomPars->Get_zOfBarrel();
	    else
	      exitBarrel.z=-_geomPars->Get_zOfBarrel();
	    exitBarrel.x=xPos*exitBarrel.z/zPos;
	    exitBarrel.y=yPos*exitBarrel.z/zPos;
	    double materialPassed=_geomHelpers->Get3dDistance(&entryBarrel,&exitBarrel);
	    double X0Passed=materialPassed/3.5;
	    double X0ToPass=_garPars->GetNLayersForSeeding()*_geomPars->Get_absThicknessBarrelLayer()[_garPars->GetNLayersForSeeding()]/3.5;
	    int leavesLayers=(int)((X0ToPass-X0Passed)/(_geomPars->Get_absThicknessBarrelLayer()[_garPars->GetNLayersForSeeding()]/3.5));
	    if(_garPars->GetDebug()>2)
	      cout << "Material Passed, X0: " << materialPassed << " , " << X0Passed << " , X0ToPass, layers): " << X0ToPass << "  " << leavesLayers << endl;
	    if((a_ext_hit->pseudoLayer) < (leavesLayers+1) && leavesLayers+1>0){
	      if(_garPars->GetDebug()>2)
		cout << "Considering hit " << a_ext_hit->hit->getCellID0() <<  " in layer " << a_ext_hit->pseudoLayer << " with barrel projection, cos = " << cosHit  << endl;
	      consideredHits.push_back(a_ext_hit);
	    }
	    else
	      if(_garPars->GetDebug()>2)
		cout << "Rejected hit " << a_ext_hit->hit->getCellID0() << " in layer " << a_ext_hit->pseudoLayer << " with barrel projection, cos = " << cosHit << endl;

	  }
	}
      }
      
      double phiMin=twopi;
      double phiMax=-twopi;
      double thetaMin=pi;
      double thetaMax=-pi;
	
      int NConsideredHits=consideredHits.size();
      if(NConsideredHits<_garPars->GetMinHitsForSeed()) {
	if(_garPars->GetDebug()>2)
	  cout << "not enough hits" << endl;
	break;
      }
      for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	ExtendedHit *a_ext_hit=consideredHits[hit_i];
	if((a_ext_hit->hit)->getPosition()[2] > 0)
	  zIsPos=true;
	else
	  zIsPos=false;
	//project  
	float x = (a_ext_hit->hit)->getPosition()[0];
	float y = (a_ext_hit->hit)->getPosition()[1];
	float z = (a_ext_hit->hit)->getPosition()[2];
	double r_phi=sqrt(x*x+y*y);
	double theta = atan2((double)r_phi,(double)z);
	double phi = atan2((double)x, (double)y);
	
	//
	if(phi<0)
	  phi=twopi+phi;
	//

	if (theta>thetaMax)
	  thetaMax=theta;
	if (theta<thetaMin)
	  thetaMin=theta;
	if (phi>phiMax)
	  phiMax=phi;
	if (phi<phiMin)
	  phiMin=phi;
      }
      thetaMin=thetaMin-(hitSizeInTheta/2);
      thetaMax=thetaMax+(hitSizeInTheta/2);
      phiMin=phiMin-(hitSizeInPhi/2);
      phiMax=phiMax+(hitSizeInPhi/2);
      int phiBins=(int)(fabs((phiMax-phiMin))/hitSizeInPhi);
      int thetaBins=(int)(fabs((thetaMax-thetaMin))/hitSizeInTheta);
      if(phiMax>6.2 && phiMin<0.1) { // spread over seam in phi, need to build new histogram with seamless phi values
	phiSeam=true;
	if(_garPars->GetDebug()>2)
	  cout << "Building histogram without seam in PHI..." << endl;
	phiMin=2*pi;
	phiMax=-2*pi;
	for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	  ExtendedHit *a_ext_hit=consideredHits[hit_i];
	  //project  
	  float x = (a_ext_hit->hit)->getPosition()[0];
	  float y = (a_ext_hit->hit)->getPosition()[1];
	  //float z = (a_ext_hit->hit)->getPosition()[2];
	  //float r = sqrt(x*x+y*y+z*z);
	  //float z_proj = z*_geomPars->Get_rOfBarrel()/sqrt(x*x+y*y);
	  double phi = atan2((double)x, (double)y);
	  //
	  if(phi<0)
	    phi=twopi+phi;
	  //
	  
	  if(phi>=pi)
	    phi=phi-pi;
	  else
	    phi=phi+pi;

	  
	
	  //if (z_proj>zMax)
	  //  zMax=z_proj;
	  //if (z_proj<zMin)
	  //  zMin=z_proj;
	  if (phi>phiMax)
	    phiMax=phi;
	  if (phi<phiMin)
	    phiMin=phi;
	}
	phiMin=phiMin-(hitSizeInPhi/2);
	phiMax=phiMax+(hitSizeInPhi/2);
	phiBins=(int)(fabs((phiMax-phiMin))/hitSizeInPhi);
      }
      if(_garPars->GetDebug()>2 && phiSeam==true)
	  cout << "..having a seam in PHI" << endl;
      if(_garPars->GetDebug()>2) 
	cout << "Histogram parameters: " <<  phiBins << ", " << phiMin<< ", " << phiMax<< ", " <<thetaBins<< ", " <<thetaMin<< ", " <<thetaMax << endl;
      TH2F *en_proj_histo = new TH2F("en_proj_histo","en_proj_histo",phiBins,phiMin,phiMax,thetaBins,thetaMin,thetaMax);
      TH2I *hit_proj_histo = new TH2I("hit_proj_histo","hit_proj_histo",phiBins,phiMin,phiMax,thetaBins,thetaMin,thetaMax);
      
      for(int hit_i=0;hit_i<NConsideredHits;hit_i++) {
	ExtendedHit *a_ext_hit=consideredHits[hit_i];
	//project  
	float x = (a_ext_hit->hit)->getPosition()[0];
	float y = (a_ext_hit->hit)->getPosition()[1];
	float z = (a_ext_hit->hit)->getPosition()[2];
	double r_phi=sqrt(x*x+y*y);
	double theta = atan2((double)r_phi,(double)z);
	double phi = atan2((double)x, (double)y);
	//
	if(phi<0)
	  phi=twopi+phi;
	//


	if(phiSeam==true) {
	  if(phi>=pi)
	    phi=phi-pi;
	  else
	    phi=phi+pi;
	}

	float en = (a_ext_hit->hit)->getEnergy();
	en_proj_histo->Fill(phi,theta,en);
	hit_proj_histo->Fill(phi,theta);
	int bin = en_proj_histo->FindBin(phi,theta);
	  binHitMap.insert(make_pair(bin,a_ext_hit));
	if(_garPars->GetDebug()>2) 
	  cout << "Histogrammed hit ; " <<  phi << ", " << theta << ", with E=" << en << endl;
      }

      // order bins per energy deposit
      multimap<float,int>  proj_bins;
      proj_bins.clear();
      multimap<float,int>::reverse_iterator proj_bins_it;
      int nX  = en_proj_histo->GetXaxis()->GetNbins()+2;
      int NBins = en_proj_histo->GetBin(en_proj_histo->GetNbinsX(),en_proj_histo->GetNbinsY());
      for(int bin_i=1;bin_i <= NBins;bin_i++) {
	if (en_proj_histo->GetBinContent(bin_i) != 0) 
	  proj_bins.insert(make_pair(en_proj_histo->GetBinContent(bin_i),bin_i ));
      }	

      vector<vector<int> > seed_clusters;
      vector<vector<int> >::iterator seed_clusters_it;
      vector<int>::iterator cl_it;
      
      for(proj_bins_it=proj_bins.rbegin();proj_bins_it!=proj_bins.rend();proj_bins_it++) {
	int bin_i = (proj_bins_it->second);
	bool in_seed=false;
	for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	  for(cl_it=seed_clusters_it->begin();cl_it!=seed_clusters_it->end();cl_it++) {
	    if(bin_i == *cl_it) {
	      in_seed=true;
	      break;
	    }
	  }
	  if(in_seed==true)
	    break;
	  int old_bin_i = (seed_clusters_it->begin())[0];
	  int old_bin_x = old_bin_i%nX;
	  int old_bin_y = old_bin_i/nX;
	  int bin_x = bin_i%nX;
	  int bin_y = bin_i/nX;
	  if( (fabs(double(bin_x-old_bin_x))<3) && (fabs(double(bin_y-old_bin_y))<3) ) {
	    in_seed=true;
	    break;
	  }
	}
	if(in_seed==true)
	  continue;
	else {
	  vector<int> a_seed;
	  double maxHitEn = 0;
	  double seedEn = 0;
	  a_seed.push_back(bin_i);
	  //RecAddNeighbourBins(en_proj_histo, bin_i, &a_seed);
	  _clusterHelpers->IterAddNeighbourBins(en_proj_histo, bin_i, &a_seed,0);
	  // check if seed originated from more than 4 hits
	    int nOriginHits=0;
	    for(int a_bin=0; a_bin<((int)a_seed.size());a_bin++) {
	      int bin_j = a_seed[a_bin];
	      nOriginHits+=(int)(hit_proj_histo->GetBinContent(bin_j));
	      for(bHM_it=binHitMap.begin();bHM_it!=binHitMap.end();bHM_it++) {
		if(bHM_it->first!=bin_j)
		  continue;
		else {
		  float hit_en = (bHM_it->second)->hit->getEnergy();
		  if(hit_en>maxHitEn)
		    maxHitEn=hit_en;
		  seedEn+=hit_en;
		}
	      }
	    }
	    if(nOriginHits>3) {
	      seed_clusters.push_back(a_seed);
	      //	      if(seed_clusters.size()==1)
	      //	_seedEn_histo->Fill(seedEn);
	      //_seedHitEn_histo->Fill(maxHitEn);
	      //_seedEn_histo->Fill(seedEn);
	      if(_garPars->GetDebug()>2) 
		cout << "Possible seed found from  " << a_seed.size() << " bins" << endl;
	    }
	    else
	      if(_garPars->GetDebug()>2) 
		cout << "Rejected seed from " << nOriginHits << " hits" << endl;
	}
      }

      for(seed_clusters_it=seed_clusters.begin();seed_clusters_it!=seed_clusters.end();seed_clusters_it++) {
	int bin_i = (seed_clusters_it->begin())[0];
	double binPos[2];  // binPos[0]=phi, binPos[1]=theta
	binPos[0]=en_proj_histo->GetXaxis()->GetBinCenter((bin_i)%nX);
	binPos[1]=en_proj_histo->GetYaxis()->GetBinCenter((bin_i)/nX);
	// to get the seed position on the ECAL front plane we need to project back
	vec3 seedPos;	    
	//if(!(cos(binPos[1])>_geomPars->Get_cosOfBarrel() && cos(binPos[1])>-_geomPars->Get_cosOfBarrel())) { 
	//if(cos(binPos[1])<_geomPars->Get_cosOfBarrel() && cos(binPos[1])>(cos(pi-acos(_geomPars->Get_cosOfBarrel())))) {  // FIXME: this has to be more precise!

	vec3 point;
	point.z = _geomPars->Get_zOfBarrel();
	double p_phi = binPos[0];
	if(phiSeam==true) { // reverse the effect of the new Phi definition
	  if(p_phi>=pi)
	    p_phi=p_phi+pi;
	  else
	    p_phi=p_phi-pi;
	}

	point.y = abs(((point.z)*tan(binPos[1]))/sqrt((1+(tan(p_phi)*tan(p_phi)))));
	if(p_phi > (pi/2) && p_phi < (3*pi/2))
	  point.y=-1*(point.y);
	point.x = abs(point.y * tan(p_phi));
	if(p_phi > (pi) )
	  point.x=-1*(point.x);

	if(_garPars->GetDebug()>2)
	  cout << "Reference point at: " << point.x << ", " << point.y << ", " << point.z << endl;
	
	int p_gamma=(int)(((p_phi-(twopi/(2*_geomPars->Get_symmetry())))/(twopi/_geomPars->Get_symmetry()))+1);
	double p_GAMMA=p_gamma*(twopi/_geomPars->Get_symmetry());
	double p_xPos=point.x*cos(p_GAMMA)-point.y*sin(p_GAMMA);
	double p_yPos=point.y*cos(p_GAMMA)+point.x*sin(p_GAMMA);
	double p_zPos=point.z;
	
	bool barrel_proj = 1;
	
	if((abs(p_xPos)*_geomPars->Get_zOfBarrel()/p_zPos) < (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry())))) {
	  if(_garPars->GetDebug()>2) 
	    cout << "p_xPos: " << abs(p_xPos)*_geomPars->Get_zOfBarrel()/p_zPos << " , (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry()))): " << (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry()))) << endl;
	  if((p_yPos*_geomPars->Get_zOfBarrel()/p_zPos) < _geomPars->Get_rOfBarrel()) {
	    if(_garPars->GetDebug()>2) 
	      cout << "p_yPos: " << p_yPos*_geomPars->Get_zOfBarrel()/p_zPos << " , _geomPars->Get_rOfBarrel(): " << _geomPars->Get_rOfBarrel() << endl;
	    barrel_proj=0;
	  }
	}
	else {
	  double y_offset = ((abs(p_xPos) - (_geomPars->Get_rOfBarrel()*tan(twopi/(2*_geomPars->Get_symmetry())))) * tan(twopi/_geomPars->Get_symmetry()) );
	  if((p_yPos*_geomPars->Get_zOfBarrel()/p_zPos) < (_geomPars->Get_rOfBarrel()-y_offset)) {
	    if(_garPars->GetDebug()>2) 
	      cout << "p_yPos: " << abs(p_xPos)*_geomPars->Get_zOfBarrel()/p_zPos << " , _geomPars->Get_rOfBarrel() - y_offset: " << _geomPars->Get_rOfBarrel() << " - " << y_offset << endl;
	    barrel_proj=0;
	  }
	}
	
	if(barrel_proj==1) { 
	  
	    //	if(binPos[1]>acos(_geomPars->Get_cosOfBarrel())) { // Project to barrel
	  if(_garPars->GetDebug()>2) 
	    cout << "Projecting to barrel because theta=" << (binPos[1])  << ", wrt to (Barrel)=" <<acos(_geomPars->Get_cosOfBarrel()) << endl;
	    double phi = binPos[0];
	    if(phi<0) {
	      cout << "Shouldnt be here: phi from hist is negative (OVERLAP-BARREL)" << endl;
	      phi=twopi+phi;
	    }
	    if(phiSeam==true) { // reverse the effect of the new Phi definition
	      if(phi>=pi)
		phi=phi+pi;
	      else
		phi=phi-pi;
	    }
	  int gamma=(int)(((phi-(twopi/(2*_geomPars->Get_symmetry())))/(twopi/_geomPars->Get_symmetry()))+1);
	  double GAMMA=gamma*(twopi/_geomPars->Get_symmetry());
	  double PHI = phi-gamma*(twopi/_geomPars->Get_symmetry());
	  double X = _geomPars->Get_rOfBarrel()*tan(PHI);
	  double Y = _geomPars->Get_rOfBarrel();
	  double R_PHI = sqrt(X*X+Y*Y);
	  double Z = R_PHI/tan(binPos[1]);	  
	  seedPos.x = X*cos(GAMMA)+Y*sin(GAMMA);   
	  seedPos.y = Y*cos(GAMMA)-X*sin(GAMMA);  
	  seedPos.z = Z;
	}
	else { //Project to Endcap
	  if(_garPars->GetDebug()>2) 
	    cout << "Projecting to endcap because theta=" << binPos[1]  << ", wrt to (Barrel)=" <<acos(_geomPars->Get_cosOfBarrel()) << endl;
	  double theta = binPos[1];
	  double Z=_geomPars->Get_zOfEndcap();
	  if (zIsPos==true)
	    Z = _geomPars->Get_zOfEndcap();
	  else
	    Z =-(_geomPars->Get_zOfEndcap());
	  double Y=Z*tan(theta);
	  double phi = binPos[0];
	  if(phi<0){
	    cout << "Shouldnt be here: phi from hist is negative (OVERLAP-ENDCAP)" << endl;
	    phi=twopi+phi;
	  }
	  if(phiSeam==true) { // reverse the effect of the new Phi definition
	      if(phi>=pi)
		phi=phi+pi;
	      else
		phi=phi-pi;
	    }

	  //  cout << "Phi: " << phi << " , tan phi: " << tan(phi) << endl;
	  //CHANGE
	  Y = abs(((Z)*tan(theta))/sqrt((1+(tan(phi)*tan(phi)))));
	  // determine sign of y
	  if(phi > (pi/2) && phi < (3*pi/2))
	    Y=-Y;
	  double X = abs(Y * tan(phi));
	  if(phi > (pi) )
	    X=-X;
	  // CHANGE
	  /*
	  int gamma=(int)(((phi-(twopi/(2*_geomPars->Get_symmetry())))/(twopi/_geomPars->Get_symmetry()))+1);
	  float GAMMA=gamma*(twopi/_geomPars->Get_symmetry());
	  float PHI = phi-gamma*(twopi/_geomPars->Get_symmetry());
	  //float X = _geomPars->Get_rOfBarrel()*tan(PHI);
	  float X = Y*tan(PHI); // CHANGE
	  seedPos.x = X*cos(GAMMA)+Y*sin(GAMMA);   
	  seedPos.y = Y*cos(GAMMA)-X*sin(GAMMA);
	  */
	  // CHANGE
	  seedPos.x = X;
	  seedPos.y = Y;
	  // CHANGE
	  
	  if (zIsPos==true)
	    seedPos.z = _geomPars->Get_zOfEndcap();
	  else
	    seedPos.z =-(_geomPars->Get_zOfEndcap());
	  if(_garPars->GetDebug()>2)
	    cout << "z is pos: " << zIsPos << endl;;
	}
	// automatically ordered by energy deposit!
	possibleSeeds.push_back(seedPos);
	if(_garPars->GetDebug()>2) 
	  cout << "Added Seed position " << seedPos.x << ", " <<  seedPos.y << ", " << seedPos.z << ",  from bin, phi: " << binPos[0] << " , theta: " << binPos[1] << endl;
	CalorimeterHitImpl *a_seed=new CalorimeterHitImpl;
	//a_hit->setCellID0(a_sim_hit->getCellID0());
	float s_pos[3] = {seedPos.x,seedPos.y,seedPos.z};
	a_seed->setEnergy(0);
	a_seed->setPosition(s_pos);
	seed_col->addElement(a_seed);
      }

      delete en_proj_histo;
      delete hit_proj_histo;
      break;
    }
    case CLUS_LOCATION_UNKNOWN:  {
      cout << "ERROR: Cluster location is unknown!!!" << endl;
      break;
    }
    }
  if(_garPars->GetDebug()>1 && possibleSeeds.size()==0)
    cout << "No seeds found!" << endl;
}






void ECALGarlic::WriteClusters(LCEvent *evt, map<int,vector<ExtendedCluster* > *> *clusMap)
{

  map<int,vector<ExtendedCluster* > *>::iterator clusMap_it;
  LCCollectionVec *photonClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  LCCollectionVec *recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  photonClusterColl->setFlag( photonClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  if(_garPars->GetDebug()>2) 
    cout << "Created collection" << endl;
  for(clusMap_it=clusMap->begin();clusMap_it!=clusMap->end();clusMap_it++) {
    vector<ExtendedCluster* > *clusVec = (clusMap_it->second);
    _nPhotonClusters += clusVec->size();
    if(_garPars->GetDebug()>2) 
      cout << "...to hold " << clusVec->size() << " clusters from ROI " << clusMap_it->first << ", total = " << _nPhotonClusters << endl;
    for (unsigned int i = 0; i < clusVec->size(); i++) {
      // File << endl;
      ClusterImpl *a_cluster = new ClusterImpl();
      ExtendedCluster *a_ext_cluster = (*clusVec)[i];
      vector<ExtendedHit* > ext_hit_vec = a_ext_cluster->hitVec;
      int NHitsInCluster = ext_hit_vec.size();
      for(int hit_i=0;hit_i<NHitsInCluster;hit_i++) {
	ExtendedHit *a_ext_hit = dynamic_cast<ExtendedHit*>(ext_hit_vec[hit_i]);
	CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit*>(a_ext_hit->hit);
	if(a_hit)
	  a_cluster->addHit(a_hit,1);
	else 
	  cout << "The hit is not there" << endl;
      }
      a_cluster->setEnergy((a_ext_cluster->Parameters)->E_GeV);
      float position[3];
      position[0] = (a_ext_cluster->Parameters)->COGx;
      position[1] = (a_ext_cluster->Parameters)->COGy;
      position[2] = (a_ext_cluster->Parameters)->COGz;
      a_cluster->setPosition(position);
      photonClusterColl->addElement(a_cluster);
      if (_garPars->GetDebug()>1) cout << "Added cluster with " << NHitsInCluster << " hits, E= " << a_cluster->getEnergy() << endl;
      ReconstructedParticleImpl* recoPart = new ReconstructedParticleImpl();
      TVector3 mom(position);
      mom = mom.Unit();
      mom *= (a_ext_cluster->Parameters)->E_GeV;
      double momArr[3];
      momArr[0] = mom.X();
      momArr[1] = mom.Y();      
      momArr[2] = mom.Z();
      recoPart->setMomentum( momArr );
      recoPart->setEnergy((a_ext_cluster->Parameters)->E_GeV );
      recoPart->setMass( 0 );
      recoPart->setCharge( 0. );
      recoPart->setType( 22 );
      recoPart->addCluster(a_cluster);
      recparcol->addElement( recoPart );
    }

  }
  evt->addCollection(photonClusterColl, _ecalPhotonClusterCollectionName);
  evt->addCollection( recparcol , _particleCollectionName);

  if (_garPars->GetDebug()>0) {
    if(photonClusterColl->getNumberOfElements()>0)
      cout << "Appended PhotonClusterCollection" << endl;
    else
      cout << "No clusters in Event!" << endl;
  }
}



void ECALGarlic::printMrGarlic() {

  cout <<"                       `.+MY'WMa," << endl;
  cout <<"                       .Mt.4&JnHTN," << endl;
  cout <<"                       d@,o..Z,'`4WN," << endl;
  cout <<"                       .WMHHM#m,  JyMx" << endl;
  cout <<"     `     `    `       `    ,Mp4. j,M,  `        `   `    `   `    `" << endl;
  cout <<"                             `J#.t ,|`N.     `" << endl;
  cout <<"                    ` `  `` .d#!J`  S ZN. `" << endl;
  cout <<"                `     ...gMY'`,=    .h.7N,      `" << endl;
  cout <<"     `     `   ` ..gMB'7`  `         `?&,THNaJ,         `    `   `" << endl;
  cout <<"          ` `..H#'!..v^                  ?4..`'WNa, ` `" << endl;
  cout <<"       `` .JM9^  .J^ ..J^`         `  ?7i.. `G, ..TMm,  `             `" << endl;
  cout <<"     `  .H#^ ` .J' .Y` ..v''3J.  ..v7=i.  .S. Jl``n.`WN," << endl;
  cout <<"     `.H@`    .Z  .^  J'     ``nJ!      7,  ^  .L` T, .TN.`    `" << endl;
  cout <<"     .H       J`     .bJ.......dG.......JR      ,|  ?,  ?N.         `" << endl;
  cout <<"     d#      .%       4.dN#NJbCTYOd##bJSJt       4   4   M|7" << endl;
  cout <<"     db `    ,!       `J4#W'  `  `  ?#WY! ..     `.  J   M       `" << endl;
  cout <<"     JN. j   ,|         JJ%          ,.$  .Y``   ,:  Z  .#`" << endl;
  cout <<"      Wb .h`  4.       ` ?S,     `  .Y= `.M:   ` J  `! .M      `" << endl;
  cout <<"     `.Wh. 4.  S.`   `.J,   73J...v'`  .JMF    `.% ?``.Mt" << endl;
  cout <<"      ` 7MJ ?o  4,       ?=&.......Jg#N#M^    `.% `  .M^" << endl;
  cout <<"        `.TMa.7+.?o.`        ` ?7''YBB'!    ` .= `..M'  `    `      `" << endl;
  cout <<"            ?WNadG,?a.`   `     `       ``  .Y...HB^" << endl;
  cout <<"               .7'HMMMNJ,.` ``  `  `  ` ..gMHMB'!" << endl;
  cout <<"     `               `.dMMMMHW+....JWMHMZMN,                     `" << endl;
  cout <<"                  ` .d#TdNWWpWWNWWMkWHWHZ1ITN,`" << endl;
  cout <<"               `  .d#OydNWt?YC!?I!?1?7dHR+zzvWb       `    `" << endl;
  cout <<"             `   `J#Jd4HW#.zz+.+l:+z1:dWKvNOzvUb  `            `" << endl;
  cout <<"     `       `.Mm.MtdDdHW%?!v!?1v:+v?1HW$,Jy.+.M|                   `" << endl;
  cout <<"               dNdBJdxdpWr11z.++z+J1udHH1v.D1I:Jb" << endl;
  cout <<"           `  `.dBnOVHdHWb..+.`?!`vjHpWtjaJO++.d%  `    `    `" << endl;
  cout <<"     `       `.M^. ?7M##HHJvz:1z+x1HpWDzJRTmo:JN," << endl;
  cout <<"              dh,  .Gd@;;??7TYY9VYTYYWH'! 4JzGM'                 `" << endl;
  cout <<"               .M,.,d#;?;;;;;;;;;;;;+SJ.`  .B=`" << endl;
  cout <<"       `   `     ?!dM1;;?;;;;;;;;;;;?+b..,l,r         `    `          `" << endl;
  cout <<"              `   .M@;;?;??;??;++?;;?;+vYWM= `" << endl;
  cout <<"     `       ..MMMMMmx;;;;;;;?jM#;?;;;;;;d@                    `" << endl;
  cout <<"          `  dM9VOVTHdHJ?;?;jd#4Nx;?;;;;;?Hm.                       `" << endl;
  cout <<"            .MZlltltltWHp;;jMF` WN+;?;??;;?JMe`       `    `" << endl;
  cout <<"          ` ,MZllllltltd@c?d#    TNx;;;;?;jHM#    `              `" << endl;
  cout <<"     `     ` WNyttlltwAdHIjM'    `.TMNkmWHUuXML `" << endl;
  cout <<"              ?WMNUUZlltdMMF       JMSuuzuuuzM#                `" << endl;
  cout <<"                .MNOlttOtd#        .WMNkQkWHNHD`        `" << endl;
  cout <<"                  TMNgggg#^ `     `   ?''''!                 `" << endl;
  cout <<"     `               ??;                                            `" << endl;
  cout <<"                                                  `   `" << endl;
  cout <<"           `                                               `" << endl;
  cout <<"     `                                                           `" << endl;
  cout <<"                 `                `" << endl;
  cout <<"" << endl;

}


void ECALGarlic::setUpGeometry() {

  float _x_thicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _x_thicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];

  // determine geometry and define pseudo layers in ECAL only
  // this code is an exect copy of M.Thomson's PandoraPFA code
  // FIXME: hardcoded number of cells along wafer row = 18 - seems to be mostly fixed (dj)
  
  _x_bField = Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();

  // Calorimeter geometry from GEAR
  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
  const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();

  // determine pseudo-geometry of detector (determined by ECAL barrel)
  // symmetry 0 =cylinder, 1=prorotype, >1 = polygon
  _x_symmetry = pEcalBarrel.getSymmetryOrder();
  _x_rOfBarrel = (float)pEcalBarrel.getExtent()[0];
  _x_zOfBarrel = (float)pEcalBarrel.getExtent()[3];

  _x_zOfEndcap = (float)pEcalEndcap.getExtent()[2];
  _x_rInnerEcalEndcap = (float)pEcalEndcap.getExtent()[0];
  _x_rOfEndcap = (float)pEcalEndcap.getExtent()[1];

  // Determine ECAL polygon angles
  // Store radial vectors perpendicular to stave layers in _barrelStaveDir 
  // this is modified to be conform with the MOKKA Stave encoding

  if(_x_symmetry>1){
    float nFoldSymmetry = static_cast<float>(_x_symmetry);
    float phi0 = pEcalBarrel.getPhi0();

    for(int i=0;i<_x_symmetry;++i){
      float phi  = phi0 + i*twopi/nFoldSymmetry;
      vec3 staveDir = {-sin(phi),cos(phi),0.};
      if(_garPars->GetDebug()>2)
	cout << "Stave: " << i << " , Phi = " << phi << ", X = " << staveDir.x << ", Y = " << staveDir.y << endl;
      _x_barrelStaveDir.push_back(staveDir);
    }
  }  

  // Using GEAR information define PSEUDOLAYERS in BARREL
  // Internally PandoraPFA defines pseudolayers starting from 1
  //   with layer zero reserved for track projections
  // for us layer 0 is the PreShower layer so this will be the same approach
  int layer=1;
  for(int i=0;i<ecalBarrelLayout.getNLayers();++i){
    _x_positionBarrelLayer[layer] = static_cast<float>(ecalBarrelLayout.getDistance(i));
    _x_thicknessBarrelLayer[layer] = ecalBarrelLayout.getThickness(i);
    _x_absThicknessBarrelLayer[layer] = ecalBarrelLayout.getAbsorberThickness(i);
    _x_padSizeEcal[layer] = ecalBarrelLayout.getCellSize0(i);
    // now position the active layer
    double delta =(_x_absThicknessBarrelLayer[layer]+_x_thicknessBarrelLayer[layer])/2.0;
    _x_positionBarrelLayer[layer]+= delta;
    layer++;
  }

  // define PreShower layer
  _x_positionBarrelLayer[0] = _x_positionBarrelLayer[1]-_x_thicknessBarrelLayer[1];
  _x_thicknessBarrelLayer[0] = _x_thicknessBarrelLayer[1];
  _x_absThicknessBarrelLayer[0] = 0;//_absThicknessBarrelLayer[1];
  _x_padSizeEcal[0] = _x_padSizeEcal[1];

  if (_garPars->GetDebug()>2) {
    cout << "Positioning Barrel layers: " << endl;
    cout << "Barrel RMin: " << _x_rOfBarrel << endl;
    cout << "First Layer Distance: " << ecalBarrelLayout.getDistance(0) << endl;
  }
  _x_positionBarrelLayer[0] = ecalBarrelLayout.getDistance(0)+_x_firstBarrelLayerOffset; // this equals rMin
  _x_absThicknessBarrelLayer[0] = 0;
  _x_thicknessBarrelLayer[0] = ecalBarrelLayout.getAbsorberThickness(1)+_x_passiveThickness+_x_activeThickness;
  _x_positionBarrelLayer[0] += (ecalBarrelLayout.getAbsorberThickness(1)+_x_passiveThickness+(_x_activeThickness/2));
  layer=0;
  _x_padSizeEcal[layer] = ecalBarrelLayout.getCellSize0(1);

  if (_garPars->GetDebug()>2) {
    cout << "Layer " << layer << endl 
	 <<  "Thickness: " << _x_thicknessBarrelLayer[layer] << ", absThick: " << _x_absThicknessBarrelLayer[layer] << ", Position: " << _x_positionBarrelLayer[layer] << endl;
  }
  for(int i=0;i<ecalBarrelLayout.getNLayers();++i){
    layer = i+1;
    _x_positionBarrelLayer[layer] = _x_positionBarrelLayer[layer-1];
    _x_absThicknessBarrelLayer[layer] = ecalBarrelLayout.getAbsorberThickness(i);
    _x_thicknessBarrelLayer[layer] = _x_absThicknessBarrelLayer[layer]+_x_passiveThickness+_x_activeThickness;
    double delta = 0;
    if(i%2 == 0)
      delta = (_x_absThicknessBarrelLayer[layer]+_x_activeThickness);
    else
      delta = (_x_absThicknessBarrelLayer[layer]+2*_x_passiveThickness+_x_activeThickness);
    _x_positionBarrelLayer[layer]+= delta;
    _x_padSizeEcal[layer] = ecalBarrelLayout.getCellSize0(i);
  }

  int nBarrelLayers = layer+1;

  // Using GEAR information define PSEUDOLAYERS in ENDCAP
  // Internally PandoraPFA defines pseudolayers starting from 1
  //   with layer zero reserved for track projections
  // for us layer 0 is the PreShower layer so this will be the same approach
  layer=1;
  if (_garPars->GetDebug()>2) {
    cout << "Positioning Endcap layers: " << endl;
    cout << "Z of Endcap: " << _x_zOfEndcap << endl;
  }
  for(int i=0;i<ecalEndcapLayout.getNLayers();++i){
    _x_positionEndcapLayer[layer] = ecalEndcapLayout.getDistance(i);
    _x_thicknessEndcapLayer[layer] = ecalEndcapLayout.getThickness(i);
    _x_absThicknessEndcapLayer[layer] = ecalEndcapLayout.getAbsorberThickness(i);
    double delta =(_x_absThicknessEndcapLayer[layer]+_x_thicknessEndcapLayer[layer])/2.0;
    _x_positionEndcapLayer[layer]+= delta;
    if(_garPars->GetDebug()>2) 
      cout << "Layer " << layer << ", i=" << i << ":" << endl 
	   <<  "Distance from gear: " <<  ecalEndcapLayout.getDistance(i) << ", Thickness: " << ecalEndcapLayout.getThickness(i) << ", absThick: " << ecalEndcapLayout.getAbsorberThickness(i) << ", Delta: " << delta << ", Position: " << _x_positionEndcapLayer[layer] << endl;
    layer++;
  }

  // define PreShower layer
  _x_positionEndcapLayer[0] = _x_positionEndcapLayer[1]-_x_thicknessEndcapLayer[1];
  _x_thicknessEndcapLayer[0] = _x_thicknessEndcapLayer[1];
  _x_absThicknessEndcapLayer[0] = 0;//_absThicknessEndcapLayer[1];
  _x_padSizeEcal[0] = _x_padSizeEcal[1];

  // more precise calculation of ECAL layer positions
  
  //build PS layer  
  _x_positionEndcapLayer[0] = ecalEndcapLayout.getDistance(0)+_x_firstEndcapLayerOffset; // this equals _zOfEndcap
  _x_absThicknessEndcapLayer[0] = 0;
  _x_thicknessEndcapLayer[0] = ecalEndcapLayout.getAbsorberThickness(1)+_x_passiveThickness+_x_activeThickness;
  _x_positionEndcapLayer[0] += (ecalEndcapLayout.getAbsorberThickness(1)+_x_passiveThickness+(_x_activeThickness/2));
  layer=0;
  if (_garPars->GetDebug()>2) {
    cout << "Layer " << layer << endl 
	 <<  "Thickness: " << _x_thicknessEndcapLayer[layer] << ", absThick: " << _x_absThicknessEndcapLayer[layer] << ", Position: " << _x_positionEndcapLayer[layer] << endl;
  }
  for(int i=0;i<ecalEndcapLayout.getNLayers();++i){
    layer = i+1;
    _x_positionEndcapLayer[layer] = _x_positionEndcapLayer[layer-1];
    _x_absThicknessEndcapLayer[layer] = ecalEndcapLayout.getAbsorberThickness(i);
    _x_thicknessEndcapLayer[layer] = _x_absThicknessEndcapLayer[layer]+_x_passiveThickness+_x_activeThickness;
    double delta = 0;
    if(i%2 == 0)
      delta = (_x_absThicknessEndcapLayer[layer]+_x_activeThickness);
    else
      delta = (_x_absThicknessEndcapLayer[layer]+2*_x_passiveThickness+_x_activeThickness);
    _x_positionEndcapLayer[layer]+= delta;
    if (_garPars->GetDebug()>2) {
      cout << "Layer " << layer << endl 
	   <<  "Delta: " << delta <<  ", Thickness: " << _x_thicknessEndcapLayer[layer] << ", absThick: " << _x_absThicknessEndcapLayer[layer] << ", Position: " << _x_positionEndcapLayer[layer] << endl;
    }
  }

  int nEndcapLayers = layer+1;

  // Number of Pseudo Layers
  _x_nPseudoLayers = nBarrelLayers;
  if(nEndcapLayers>nBarrelLayers)_x_nPseudoLayers=nEndcapLayers;

  if (_garPars->GetDebug()>2) {
    cout << "N Pseudo Layers: " << _x_nPseudoLayers << endl;
  }

  if(_garPars->GetDebug()>2){
    for(int i=0;i<nBarrelLayers;++i){
      std::cout << "Barrel Layer " << i << " : " << _x_positionBarrelLayer[i] << ", pad Size: " << _x_padSizeEcal[i] <<std::endl; 
    }
    for(int i=0;i<nBarrelLayers;++i){
      std::cout << "EndCap Layer " << i << " : " << _x_positionEndcapLayer[i] << ", pad Size: " << _x_padSizeEcal[i] <<std::endl; 
    }
  }
    
  //redifinition of _zOfEndcap and _rOfBarrel to respective PS layers
  _x_zOfEndcap = _x_positionEndcapLayer[0];    
  _x_rOfBarrel = _x_positionBarrelLayer[0];

  if (_garPars->GetDebug()>2) {
    cout << "Corrected Z of Endcap: " << _x_zOfEndcap << endl;
    cout << "Corrected R of Barrel: " << _x_rOfBarrel << endl;
  }

  _x_cosOfBarrel = _x_zOfBarrel/(sqrt(_x_rOfBarrel*_x_rOfBarrel+_x_zOfBarrel*_x_zOfBarrel));
  if (_garPars->GetDebug()>2) {
    cout << "Cos Of Barrel: " << _x_cosOfBarrel << endl;
  }

  return;
}

//  LocalWords:  ECALGarlicAlgorithmParameters
