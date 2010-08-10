#include "ECALGarlicConstants.hh"
#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryParameters.hh"
#include "ECALGarlicClusterHelpers.hh"
#include "ECALGarlicMLP.hh"

#include "TH2F.h"
#include "TVector3.h"
#include "TMatrixD.h"

using namespace ECALGarlicConstants;

#include <ClusterShapes.h>

#include <cmath>

#include <math.h>

//using namespace ECALGarlicAlgorithmParameters;

using std::cout;
using std::endl;
using std::vector;
using std::map;

void ECALGarlicClusterHelpers::setupNN() {

  vector<std::string> varnames;
  varnames.push_back("ClusterParameters.start");
  varnames.push_back("ClusterParameters.end");
  varnames.push_back("ClusterParameters.depth");
  varnames.push_back("ClusterParameters.E1C/ClusterParameters.Etot_g");
  varnames.push_back("ClusterParameters.E4C/ClusterParameters.Etot_g");
  varnames.push_back("ClusterParameters.E9C/ClusterParameters.Etot_g");
  varnames.push_back("ClusterParameters.dirErr");
  varnames.push_back("ClusterParameters.seedDirErr");
  varnames.push_back("ClusterParameters.Width");
  varnames.push_back("ClusterParameters.Eccentricity");
  varnames.push_back("ClusterParameters.nHits");
  varnames.push_back("ClusterParameters.Volume");
  varnames.push_back("ClusterParameters.depth-ClusterParameters.start");

  vector<std::string> varnames2;
  varnames2.push_back("ClusterParameters.start");
  varnames2.push_back("ClusterParameters.end");
  varnames2.push_back("ClusterParameters.depth");
  varnames2.push_back("ClusterParameters.E1C/ClusterParameters.Etot_g");
  varnames2.push_back("ClusterParameters.E4C/ClusterParameters.Etot_g");
  varnames2.push_back("ClusterParameters.E9C/ClusterParameters.Etot_g");
  varnames2.push_back("ClusterParameters.dirErr");
  varnames2.push_back("ClusterParameters.seedDirErr");
  varnames2.push_back("ClusterParameters.Es3/ClusterParameters.Etot_g");
  varnames2.push_back("ClusterParameters.Width");
  varnames2.push_back("ClusterParameters.Eccentricity");
  varnames2.push_back("ClusterParameters.nHits");
  varnames2.push_back("ClusterParameters.Volume");
  varnames2.push_back("ClusterParameters.depth-ClusterParameters.start");

  int nCellsPerWafer = _geomParams->Get_nCellsPerWafer();

  if(nCellsPerWafer == 18) {
    MLPResponse0_25_B = new ReadMLP0_25_B(varnames);
    MLPResponse0_35_B = new ReadMLP0_35_B(varnames);
    MLPResponse0_5_B = new ReadMLP0_5_B(varnames);
    MLPResponse0_75_B = new ReadMLP0_75_B(varnames);
    MLPResponse1_B = new ReadMLP1_B(varnames);
    MLPResponse1_25_B = new ReadMLP1_25_B(varnames);
    MLPResponse1_5_B = new ReadMLP1_5_B(varnames);
    MLPResponse1_75_B = new ReadMLP1_75_B(varnames);
    MLPResponse2_25_B = new ReadMLP2_25_B(varnames);
    MLPResponse3_B = new ReadMLP3_B(varnames);
    MLPResponse5_B = new ReadMLP5_B(varnames2);
    MLPResponse10_B = new ReadMLP10_B(varnames2);
    MLPResponse20_B = new ReadMLP20_B(varnames2);
    //MLPResponse50_B = new ReadMLP50_B(varnames2);

    MLPResponse0_25_EC = new ReadMLP0_25_EC(varnames);
    MLPResponse0_35_EC = new ReadMLP0_35_EC(varnames);
    MLPResponse0_5_EC = new ReadMLP0_5_EC(varnames);
    MLPResponse0_75_EC = new ReadMLP0_75_EC(varnames);
    MLPResponse1_EC = new ReadMLP1_EC(varnames);
    MLPResponse1_25_EC = new ReadMLP1_25_EC(varnames);
    MLPResponse1_5_EC = new ReadMLP1_5_EC(varnames);
    MLPResponse1_75_EC = new ReadMLP1_75_EC(varnames);
    MLPResponse2_25_EC = new ReadMLP2_25_EC(varnames);
    MLPResponse3_EC = new ReadMLP3_EC(varnames);
    MLPResponse5_EC = new ReadMLP5_EC(varnames2);
    MLPResponse10_EC = new ReadMLP10_EC(varnames2);
    MLPResponse20_EC = new ReadMLP20_EC(varnames2);
    //MLPResponse50_EC = new ReadMLP50_EC(varnames2);
  } else if (nCellsPerWafer == 5) {
    MLPResponse0_25_B = new ReadMLP0_25_B_1x1(varnames);
   
    MLPResponse0_35_B = new ReadMLP0_35_B_1x1(varnames);
    MLPResponse0_5_B = new ReadMLP0_5_B_1x1(varnames);
    MLPResponse0_75_B = new ReadMLP0_75_B_1x1(varnames);
    MLPResponse1_B = new ReadMLP1_B_1x1(varnames);
    MLPResponse1_25_B = new ReadMLP1_25_B_1x1(varnames);
    MLPResponse1_5_B = new ReadMLP1_5_B_1x1(varnames);
    MLPResponse1_75_B = new ReadMLP1_75_B_1x1(varnames);
    MLPResponse2_25_B = new ReadMLP2_25_B_1x1(varnames);
    MLPResponse3_B = new ReadMLP3_B_1x1(varnames);
    MLPResponse10_B = new ReadMLP10_B_1x1(varnames2);

    MLPResponse0_25_EC = new ReadMLP0_25_EC_1x1(varnames);
    MLPResponse0_35_EC = new ReadMLP0_35_EC_1x1(varnames);
    MLPResponse0_5_EC = new ReadMLP0_5_EC_1x1(varnames);
    MLPResponse0_75_EC = new ReadMLP0_75_EC_1x1(varnames);
    MLPResponse1_EC = new ReadMLP1_EC_1x1(varnames);
    MLPResponse1_25_EC = new ReadMLP1_25_EC_1x1(varnames);
    MLPResponse1_5_EC = new ReadMLP1_5_EC_1x1(varnames);
    MLPResponse1_75_EC = new ReadMLP1_75_EC_1x1(varnames);
    MLPResponse2_25_EC = new ReadMLP2_25_EC_1x1(varnames);
    MLPResponse3_EC = new ReadMLP3_EC_1x1(varnames);
    MLPResponse10_EC = new ReadMLP10_EC_1x1(varnames2);
  } else if (nCellsPerWafer == 4) {
    MLPResponse0_25_B = new ReadMLP0_25_B_2x2(varnames);
    MLPResponse0_35_B = new ReadMLP0_35_B_2x2(varnames);
    MLPResponse0_5_B = new ReadMLP0_5_B_2x2(varnames);
    MLPResponse0_75_B = new ReadMLP0_75_B_2x2(varnames);
    MLPResponse1_B = new ReadMLP1_B_2x2(varnames);
    MLPResponse1_25_B = new ReadMLP1_25_B_2x2(varnames);
    MLPResponse1_5_B = new ReadMLP1_5_B_2x2(varnames);
    MLPResponse1_75_B = new ReadMLP1_75_B_2x2(varnames);
    MLPResponse2_25_B = new ReadMLP2_25_B_2x2(varnames);
    MLPResponse3_B = new ReadMLP3_B_2x2(varnames);
    MLPResponse10_B = new ReadMLP10_B_2x2(varnames2);
      
    MLPResponse0_25_EC = new ReadMLP0_25_EC_2x2(varnames);
    MLPResponse0_35_EC = new ReadMLP0_35_EC_2x2(varnames);
    MLPResponse0_5_EC = new ReadMLP0_5_EC_2x2(varnames);
    MLPResponse0_75_EC = new ReadMLP0_75_EC_2x2(varnames);
    MLPResponse1_EC = new ReadMLP1_EC_2x2(varnames);
    MLPResponse1_25_EC = new ReadMLP1_25_EC_2x2(varnames);
    MLPResponse1_5_EC = new ReadMLP1_5_EC_2x2(varnames);
    MLPResponse1_75_EC = new ReadMLP1_75_EC_2x2(varnames);
    MLPResponse2_25_EC = new ReadMLP2_25_EC_2x2(varnames);
    MLPResponse3_EC = new ReadMLP3_EC_2x2(varnames);
    MLPResponse10_EC = new ReadMLP10_EC_2x2(varnames2);
  } else if (nCellsPerWafer == 9) {
    MLPResponse0_25_B = new ReadMLP0_25_B_10x10(varnames);
    MLPResponse0_35_B = new ReadMLP0_35_B_10x10(varnames);
    MLPResponse0_5_B = new ReadMLP0_5_B_10x10(varnames);
    MLPResponse0_75_B = new ReadMLP0_75_B_10x10(varnames);
    MLPResponse1_B = new ReadMLP1_B_10x10(varnames);
    MLPResponse1_25_B = new ReadMLP1_25_B_10x10(varnames);
    MLPResponse1_5_B = new ReadMLP1_5_B_10x10(varnames);
    MLPResponse1_75_B = new ReadMLP1_75_B_10x10(varnames);
    MLPResponse2_25_B = new ReadMLP2_25_B_10x10(varnames);
    MLPResponse3_B = new ReadMLP3_B_10x10(varnames);
    MLPResponse10_B = new ReadMLP10_B_10x10(varnames2);

    MLPResponse0_25_EC = new ReadMLP0_25_EC_10x10(varnames);
    MLPResponse0_35_EC = new ReadMLP0_35_EC_10x10(varnames);
    MLPResponse0_5_EC = new ReadMLP0_5_EC_10x10(varnames);
    MLPResponse0_75_EC = new ReadMLP0_75_EC_10x10(varnames);
    MLPResponse1_EC = new ReadMLP1_EC_10x10(varnames);
    MLPResponse1_25_EC = new ReadMLP1_25_EC_10x10(varnames);
    MLPResponse1_5_EC = new ReadMLP1_5_EC_10x10(varnames);
    MLPResponse1_75_EC = new ReadMLP1_75_EC_10x10(varnames);
    MLPResponse2_25_EC = new ReadMLP2_25_EC_10x10(varnames);
    MLPResponse3_EC = new ReadMLP3_EC_10x10(varnames);
    MLPResponse10_EC = new ReadMLP10_EC_10x10(varnames2);
  }

  _nn_is_setup=true;

  return;
}



void ECALGarlicClusterHelpers::FillClusterParameters(LCEvent *evt, 
						     ExtendedCluster *myCluster, 
						     vector<ExtendedCluster*> *clusters, 
						     int clus_ID, 
						     ClusterParameters *clusPar,
						     vector<ExtendedTrack*> tracks) {

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  clusPar->ID=clus_ID;
  LCCollection *preClusColl = 0; 
  preClusColl = evt->getCollection(_algoParams->GetEcalPreClusterCollectionName());
  CellIDDecoder<CalorimeterHit> decoder(preClusColl);  
  
  // 1.) energy distribution in structures
  double en_s1=0;
  double en_s2=0;
  double en_s3=0;
  // 2.) cluster center of gravity
  vec3 cluster_cog;
  cluster_cog.x=0;
  cluster_cog.y=0;
  cluster_cog.z=0;
  double en_norm=0;
  // 3.) longitudinal shower profile
  //TH1F *long_sh_profile=new TH1F("LongitudinalShowerProfile","LongitudinalShowerProfile",30,0,30);
  // 4.) start/end shower
  double min_X0=50;
  double max_X0=0;
  // 5.) E4C/Etot
  //  vec3 *seed_dir=&(myCluster->seededFrom);
  // 6.) E1C/Etot
  double E1C_en=0;
  double E4C_en=0;
  double E9C_en=0;

  // 7.) nearest track
  double dist_track=9999;
  double shortest_dist_track=9999;
  // 8.) Etot/Multiplicity
  // 9.) Mean Shower Depth
  // 10.) Cluster Zone
  clusPar->zone=0;
  
  bool hasBarrelHits=0;
  bool hasEndcapHits=0;
  vector<ExtendedHit* > *clusterHits = &(myCluster->hitVec);
  int NClusteredHits = clusterHits->size();
  // calculate COG
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    en_norm+=hit_en;
    cluster_cog.x+=(a_hit->getPosition()[0])*hit_en;
    cluster_cog.y+=(a_hit->getPosition()[1])*hit_en;
    cluster_cog.z+=(a_hit->getPosition()[2])*hit_en;
   int module = decoder(a_hit)["M"];
    if(0 < module && module < 6)
      hasBarrelHits=1;
    else
      hasEndcapHits=1;
  }
  int NCountedGhostHits=0;
  if(myCluster->Ghosts!=0) {
    GhostCluster *ghostHits = myCluster->Ghosts;
    if(ghostHits!=0) {
      vector<GhostHit*> ghostHitVec = (ghostHits->ghostHitVec);
      int NGhostHits = ghostHitVec.size();
      if(NGhostHits>0) {
	//	cout << "Ghost cluster hits: " << NGhostHits << endl;
	for(int g_i=0;g_i<NGhostHits;g_i++) {
	  GhostHit *g_hit= dynamic_cast<GhostHit*>(ghostHitVec[g_i]);
	  float hit_en=g_hit->Energy;
	  cluster_cog.x+=(g_hit->Position).x*hit_en;
	  cluster_cog.y+=(g_hit->Position).y*hit_en;
	  cluster_cog.z+=(g_hit->Position).z*hit_en;
	  en_norm+=hit_en;
	  if(g_hit->Count==1)
	    NCountedGhostHits++;
	}
      }
    }
  }
  cluster_cog.x=cluster_cog.x/en_norm;
  cluster_cog.y=cluster_cog.y/en_norm;
  cluster_cog.z=cluster_cog.z/en_norm;
  clusPar->COGx=cluster_cog.x;
  clusPar->COGy=cluster_cog.y;
  clusPar->COGz=cluster_cog.z;
  if(_algoParams->GetDebug()>2)
    cout << "Cluster COG is at : " << cluster_cog.x << ", " << cluster_cog.y << ", " << cluster_cog.z << endl;
  if(hasBarrelHits==1 && hasEndcapHits==0)
    clusPar->zone=1;
  if(hasBarrelHits==0 && hasEndcapHits==1)
    clusPar->zone=2;
  if(hasBarrelHits==1 && hasEndcapHits==1)
    clusPar->zone=3;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    //if(myHit->preShower==1)
    //  continue;
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    vec3 hitPos;
    hitPos.x=a_hit->getPosition()[0];
    hitPos.y=a_hit->getPosition()[1];
    hitPos.z=a_hit->getPosition()[2];
    //int cell_i = decoder(a_hit)["I"];
    //int cell_j = decoder(a_hit)["J"];
    int layer = decoder(a_hit)["K-1"];
    int psLayer = myHit->pseudoLayer;
    double lambda=(cluster_cog.x*hitPos.x)+(cluster_cog.y*hitPos.y)+(cluster_cog.z*hitPos.z);// calculate distance to line from IP through seed: "direction"
    double lambda_prime=fabs(lambda/((cluster_cog.x*cluster_cog.x)+(cluster_cog.y*cluster_cog.y)+(cluster_cog.z*cluster_cog.z)));
    vec3 lotpunkt;
    lotpunkt.x=cluster_cog.x*lambda_prime;
    lotpunkt.y=cluster_cog.y*lambda_prime;
    lotpunkt.z=cluster_cog.z*lambda_prime;
    double a_dist=_geomHelper->Get3dDistance(&lotpunkt,&hitPos);
    if(a_dist<(_geomParams->Get_padSizeEcal() [1]/2))
      E1C_en+=hit_en;
    if(a_dist<(sqrt(2.)*_geomParams->Get_padSizeEcal() [1]))
      E4C_en+=hit_en;
    if(a_dist<(sqrt(8.)*_geomParams->Get_padSizeEcal() [1]))
      E9C_en+=hit_en;   
    double X0_passed=0; // calculate X0 in front of hit
    if(hitPos.z<_geomParams->Get_zOfBarrel() && hitPos.z>-_geomParams->Get_zOfBarrel()) { //Barrel Hit
      int stave = decoder(a_hit)["S-1"];
      //float stave_angle=acos(_geomParams->Get_barrelStaveDir()[stave].x);
      double stave_x=_geomParams->Get_barrelStaveDir()[stave].x;
      double stave_y=_geomParams->Get_barrelStaveDir()[stave].y;
      float r_frontplane=_geomParams->Get_rOfBarrel();
      double lambda=r_frontplane*r_frontplane;
      if(_algoParams->GetDebug()>2) {
	cout << "Hit at " << hitPos.x << ", " <<  hitPos.y << ", " << hitPos.z << endl;
	//cout << "Stave angle: " << stave_angle << ", r_frontplane: " << r_frontplane << endl;
	cout << "Stave: " << stave << " , r_frontplane* stave.x/.y: " << r_frontplane*stave_x << ", " << r_frontplane*stave_y << endl;
      }
      double lambda_prime=fabs(lambda/(hitPos.x*r_frontplane*stave_x+hitPos.y*r_frontplane*stave_y));
      if(_algoParams->GetDebug()>2)
	cout << "lambda_prime: " << lambda_prime << endl;
      vec3 frontplaneProj;
      frontplaneProj.x=hitPos.x*lambda_prime;
      frontplaneProj.y=hitPos.y*lambda_prime;
      frontplaneProj.z=hitPos.z*lambda_prime;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
      if(_algoParams->GetDebug()>2)
	cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
      double abs_th=0;
      int lay_i=psLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	lay_i--;
      }
      double X0_ratio = 0;
      if(psLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[psLayer])-_geomParams->Get_rOfBarrel());
      else
	X0_ratio =0;
      X0_passed=(X0_ratio*dist)/3.5;
      if(_algoParams->GetDebug()>2)
	cout << "X0 passed: " << X0_passed << endl;
    }
    else { //Endcap hit...
      double phi_pos = atan2((double)hitPos.x,(double)hitPos.y);
      int gamma;
      if(phi_pos<0) 
	phi_pos=twopi+phi_pos;
      gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
      double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
      double xPos=hitPos.x*cos(GAMMA)-hitPos.y*sin(GAMMA);
      double yPos=hitPos.y*cos(GAMMA)+hitPos.x*sin(GAMMA);
      double zPos=hitPos.z;
      //      float cosHit=zPos/sqrt(yPos*yPos+zPos*zPos);
      //if(!(cosHit<_cosOfBarrel && cosHit>-_cosOfBarrel && clusPar->zone==2)) { 
      bool barrel_proj = 1;
      if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry()))) ) {
	if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
	  barrel_proj=0;
      }
      else {
	double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
	if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
	  barrel_proj=0;
      }
      
      if(barrel_proj==0 && clusPar->zone==2) { // Endcap hit without Barrel projection

	// if(cosHit<_cosOfBarrel) { // ...without Barrel Projection and cluster in overlap zone
	double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	//	if(z_ratio>1)
	//	  cout << "Warning: z_ratio >1: "<< z_ratio << endl;
	vec3 frontplaneProj;
	frontplaneProj.x=hitPos.x*z_ratio;
	frontplaneProj.y=hitPos.y*z_ratio;
	frontplaneProj.z=hitPos.z*z_ratio;
	double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	double abs_th=0;
	int lay_i=psLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	  lay_i--;
	}
	double X0_ratio = 0;
	if(psLayer>0)
	  X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	else
	  X0_ratio = 0;
	X0_passed=(X0_ratio*dist)/3.5;
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	if(_algoParams->GetDebug()>2)
	  cout << "X0 passed: " << X0_passed << endl;
      }
      else { //...with Barrel Projection
	if(_algoParams->GetDebug()>2)
	  cout << "hit with barrel projection" << endl;
	vec3 entryBarrel;   //  first get material passed in barrel
	entryBarrel.y=_geomParams->Get_rOfBarrel();
	entryBarrel.x=xPos*entryBarrel.y/yPos;
	entryBarrel.z=zPos*entryBarrel.y/yPos;
	vec3 exitBarrel;
	if(zPos>0)
	  exitBarrel.z=_geomParams->Get_zOfBarrel();
	else
	  exitBarrel.z=-_geomParams->Get_zOfBarrel();
	exitBarrel.x=xPos*exitBarrel.z/zPos;
	exitBarrel.y=yPos*exitBarrel.z/zPos;
	double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
	int   bestBarrelLayer=0;
	double barrelLayerDist = 9999.;
	double prodRadiusMax=0.;
	for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
	  double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
	  if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
	}
	for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	  if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	    barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	    bestBarrelLayer = ilayer;
	  }
	}
	double abs_th=0;
	int lay_i=bestBarrelLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	  lay_i--;
	}
	double X0_ratio = 0;
	if(bestBarrelLayer>0)
	  X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
	else
	  X0_ratio = 0;
	double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	// then adding "normal"endcap material passed
	double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	vec3 frontplaneProj;
	frontplaneProj.x=hitPos.x*z_ratio;
	frontplaneProj.y=hitPos.y*z_ratio;
	frontplaneProj.z=hitPos.z*z_ratio;
	double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	abs_th=0;
	lay_i=psLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	  lay_i--;
	}
	if(psLayer>0)
	  X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	else
	  X0_ratio = 0;
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << layer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	X0_passed=((X0_ratio*dist)/3.5)+barrelX0Passed;
	if(_algoParams->GetDebug()>2)
	  cout << "X0 passed: " << X0_passed << endl;
      }
    }
    if(X0_passed>max_X0)
      max_X0=X0_passed;
    if(X0_passed<min_X0) {
      min_X0=fabs(X0_passed);
    }
    //if(X0_passed>30)
    //  cout << "Still sth wrong in X0 passed, " << X0_passed << endl;
    //long_sh_profile->Fill(X0_passed,hit_en);
    if(layer<10)
      en_s1+=hit_en;
    else {
      if(layer<25)
	en_s2+=hit_en;
      else
	en_s3+=hit_en;
    }
  }
  if(myCluster->Ghosts!=0) {
    GhostCluster *ghostHits = myCluster->Ghosts;
    if(ghostHits!=0) {
      vector<GhostHit*> ghostHitVec = (ghostHits->ghostHitVec);
      int NGhostHits = ghostHitVec.size();
      if(NGhostHits>0) {
	//	cout << "Ghost cluster hits: " << NGhostHits << endl;
	for(int g_i=0;g_i<NGhostHits;g_i++) {
	  GhostHit *g_hit= dynamic_cast<GhostHit*>(ghostHitVec[g_i]);
	  float hit_en=g_hit->Energy;
	  int layer = g_hit->Layer;
	  int psLayer = g_hit->PseudoLayer;
	  vec3 hitPos;
	  hitPos.x=(g_hit->Position).x;
	  hitPos.y=(g_hit->Position).y;
	  hitPos.z=(g_hit->Position).z;
	  double lambda=(cluster_cog.x*hitPos.x)+(cluster_cog.y*hitPos.y)+(cluster_cog.z*hitPos.z);// calculate distance to line from IP through seed: "direction"
	  double lambda_prime=fabs(lambda/((cluster_cog.x*cluster_cog.x)+(cluster_cog.y*cluster_cog.y)+(cluster_cog.z*cluster_cog.z)));
	  vec3 lotpunkt;
	  lotpunkt.x=cluster_cog.x*lambda_prime;
	  lotpunkt.y=cluster_cog.y*lambda_prime;
	  lotpunkt.z=cluster_cog.z*lambda_prime;
	  double a_dist=_geomHelper->Get3dDistance(&lotpunkt,&hitPos);
	  if(a_dist<(_geomParams->Get_padSizeEcal() [1]/2))
	    E1C_en+=hit_en;
	  if(a_dist<(sqrt(2.)*_geomParams->Get_padSizeEcal() [1]))
	    E4C_en+=hit_en;
	  if(a_dist<(sqrt(8.)*_geomParams->Get_padSizeEcal() [1]))
	    E9C_en+=hit_en;
	  double X0_passed=0; // calculate X0 in front of hit
	  if(hitPos.z<_geomParams->Get_zOfBarrel() && hitPos.z>-_geomParams->Get_zOfBarrel()) { //Barrel Hit
	    int stave = g_hit->Stave;
	    double stave_x=_geomParams->Get_barrelStaveDir()[stave].x;
	    double stave_y=_geomParams->Get_barrelStaveDir()[stave].y;
	    double r_frontplane=_geomParams->Get_rOfBarrel();
	    double lambda=r_frontplane*r_frontplane;
	    if(_algoParams->GetDebug()>2) {
	      cout << "Hit at " << hitPos.x << ", " <<  hitPos.y << ", " << hitPos.z << endl;
	      //cout << "Stave angle: " << stave_angle << ", r_frontplane: " << r_frontplane << endl;
	      cout << "Stave: " << stave << " , r_frontplane* stave.x/.y: " << r_frontplane*stave_x << ", " << r_frontplane*stave_y << endl;
	    }
	    double lambda_prime=fabs(lambda/(hitPos.x*r_frontplane*stave_x+hitPos.y*r_frontplane*stave_y));
	    if(_algoParams->GetDebug()>2)
	      cout << "lambda_prime: " << lambda_prime << endl;
	    vec3 frontplaneProj;
	    frontplaneProj.x=hitPos.x*lambda_prime;
	    frontplaneProj.y=hitPos.y*lambda_prime;
	    frontplaneProj.z=hitPos.z*lambda_prime;
	    double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	    if(_algoParams->GetDebug()>2)
	      cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
	    double abs_th=0;
	    int lay_i=psLayer;
	    while(lay_i>0) {
	      abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	      lay_i--;
	    }
	    double X0_ratio = 0;
	    if(psLayer>0)
	      X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[psLayer])-_geomParams->Get_rOfBarrel());
	    else 
	      X0_ratio = 0;
	    X0_passed=(X0_ratio*dist)/3.5;
	    if(_algoParams->GetDebug()>2)
	      cout << "X0 passed: " << X0_passed << endl;
	  }
	  else { //Endcap hit...
	    double phi_pos = atan2((double)hitPos.x,(double)hitPos.y);
	    int gamma;
	    if(phi_pos<0) 
	      phi_pos=twopi+phi_pos;
	    gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
	    double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
	    double xPos=hitPos.x*cos(GAMMA)-hitPos.y*sin(GAMMA);
	    double yPos=hitPos.y*cos(GAMMA)+hitPos.x*sin(GAMMA);
	    double zPos=hitPos.z;
	    // float cosHit=zPos/sqrt(yPos*yPos+zPos*zPos);
	    //if(!(cosHit<_cosOfBarrel && cosHit>-_cosOfBarrel && clusPar->zone==2)) { 
	    // if(cosHit<_cosOfBarrel) { // ...without Barrel Projection and cluster in overlap zone
	    bool barrel_proj = 1;
	    
	    if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) {
	      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
		barrel_proj=0;
	    }
	    else {
	      double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
	      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
		barrel_proj=0;
	    }
	    
	    if(barrel_proj==0 && clusPar->zone==2) { // Endcap hit without Barrel projection
	      
	      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	      //	if(z_ratio>1)
	      //	  cout << "Warning: z_ratio >1: " << z_ratio << endl;
	      vec3 frontplaneProj;
	      frontplaneProj.x=hitPos.x*z_ratio;
	      frontplaneProj.y=hitPos.y*z_ratio;
	      frontplaneProj.z=hitPos.z*z_ratio;
	      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	      double abs_th=0;
	      int lay_i=psLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
		lay_i--;
	      }
	      double X0_ratio = 0;
	      if(psLayer>0)
		X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	      else
		X0_ratio = 0;
	      X0_passed=(X0_ratio*dist)/3.5;
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      if(_algoParams->GetDebug()>2)
		cout << "X0 passed: " << X0_passed << endl;
	    }
	    else { //...with Barrel Projection
	      if(_algoParams->GetDebug()>2)
		cout << "hit with barrel projection" << endl;
	      vec3 entryBarrel;   //  first get material passed in barrel
	      entryBarrel.y=_geomParams->Get_rOfBarrel();
	      entryBarrel.x=xPos*entryBarrel.y/yPos;
	      entryBarrel.z=zPos*entryBarrel.y/yPos;
	      vec3 exitBarrel;
	      if(zPos>0)
		exitBarrel.z=_geomParams->Get_zOfBarrel();
	      else
		exitBarrel.z=-_geomParams->Get_zOfBarrel();
	      exitBarrel.x=xPos*exitBarrel.z/zPos;
	      exitBarrel.y=yPos*exitBarrel.z/zPos;
	      double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
	      int   bestBarrelLayer=0;
	      double barrelLayerDist = 9999.;
	      double prodRadiusMax=0.;
	      for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
		double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
		if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
	      }
	      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
		if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
		  barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
		  bestBarrelLayer = ilayer;
		}
	      }
	      if(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer]<0 && bestBarrelLayer>0) {
		bestBarrelLayer--;   
	      }
	      double leftBarrelDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
	      
	      double abs_th=0;
	      int lay_i=bestBarrelLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
		lay_i--;
	      }
	      abs_th+=leftBarrelDist;
	      double X0_ratio = 0;
	      if(bestBarrelLayer>0)
		X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
	      else
		X0_ratio = 0;
	      double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      // then adding "normal"endcap material passed
	      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	      vec3 frontplaneProj;
	      frontplaneProj.x=hitPos.x*z_ratio;
	      frontplaneProj.y=hitPos.y*z_ratio;
	      frontplaneProj.z=hitPos.z*z_ratio;
	      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	      abs_th=0;
	      lay_i=psLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
		lay_i--;
	      }
	      if(psLayer>0)
		X0_ratio=abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	      else 
		X0_ratio = 0;
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      X0_passed=((X0_ratio*dist)/3.5)+barrelX0Passed;
	      if(_algoParams->GetDebug()>2)
		cout << "X0 passed: " << X0_passed << endl;
	    }
	  }
	  if(layer<10)
	    en_s1+=hit_en;
	  else {
	    if(layer<25)
	      en_s2+=hit_en;
	    else
	      en_s3+=hit_en;
	  }
	}
      }
    }
  }

  clusPar->start=min_X0;
  clusPar->end=max_X0;
  double Etot=en_s1+en_s2+en_s3;
  clusPar->Etot=Etot;
  clusPar->Es1=en_s1;
  clusPar->Es2=en_s2;
  clusPar->Es3=en_s3;
  clusPar->E1C=E1C_en;
  clusPar->E4C=E4C_en;
  clusPar->E9C=E9C_en;

  clusPar->distToBiggest=9999;
  clusPar->smallestDistToBiggest=9999;

  double r_phi=sqrt(cluster_cog.x*cluster_cog.x+cluster_cog.y*cluster_cog.y);
  double theta_cog = atan2((double)r_phi,(double)cluster_cog.z);
  clusPar->theta=theta_cog;
  // get mean shower depth from COG: same procedure as for hits, FIXME: be more precise?
  double phi_cog = atan2((double)cluster_cog.x,(double)cluster_cog.y);
  int gamma_cog;
  if(phi_cog<0) 
    phi_cog=twopi+phi_cog;
  clusPar->phi=phi_cog;
  gamma_cog=(int)(((phi_cog-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
  double GAMMA_cog=gamma_cog*(twopi/_geomParams->Get_symmetry());
  double X0_cog=0;
  if(cluster_cog.z>-_geomParams->Get_zOfBarrel() && cluster_cog.z<_geomParams->Get_zOfBarrel()) { // in Barrel 
    // find best barrel layer
    int   bestBarrelLayer=0;
    double barrelLayerDist = 9999.;
    double prodRadiusMax=0.;
    for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
      double prodRadius = cluster_cog.x*_geomParams->Get_barrelStaveDir()[istave].x+cluster_cog.y*_geomParams->Get_barrelStaveDir()[istave].y;
      if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
    }
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	bestBarrelLayer = ilayer;
      }
    }
    if((prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])<0 && bestBarrelLayer>0) {
      bestBarrelLayer--;   
    }
    double leftDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
    if(_algoParams->GetDebug()>2)
      cout << "Accounted layer " << bestBarrelLayer << " for COG, leaving distance: " << leftDist << endl;
    int layer=bestBarrelLayer;
    double r_frontplane=_geomParams->Get_rOfBarrel();
    double lambda=r_frontplane*r_frontplane;
    if(_algoParams->GetDebug()>2) {
      cout << "GAMMA_COG: " << GAMMA_cog << ", r_frontplane: " << r_frontplane << endl;
      cout << "r_frontplane* (sin,cos): " << r_frontplane*sin(GAMMA_cog) << ", " << r_frontplane*cos(GAMMA_cog) << endl;
    }
    double lambda_prime=fabs(lambda/(cluster_cog.x*r_frontplane*sin(GAMMA_cog)+cluster_cog.y*r_frontplane*cos(GAMMA_cog)));
    if(_algoParams->GetDebug()>2)
      cout << "lambda_prime: " << lambda_prime << endl;
    vec3 frontplaneProj;
    frontplaneProj.x=cluster_cog.x*lambda_prime;
    frontplaneProj.y=cluster_cog.y*lambda_prime;
    frontplaneProj.z=cluster_cog.z*lambda_prime;
    double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
    if(_algoParams->GetDebug()>2)
      cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
    double abs_th=0;
    int lay_i=layer;
    while(lay_i>0) {
      abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
      lay_i--;
    }
    abs_th+=leftDist;
    double X0_ratio = 0;
    if(bestBarrelLayer>0)
      X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
    else
      X0_ratio = abs_th;
    X0_cog=(X0_ratio*dist)/3.5;
    if(_algoParams->GetDebug()>2)
      cout << "Passed " << layer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
  }
  else { //in Endcap ...
    // determine ENDCAP layer
    int   bestEndcapLayer=0;
    double endcapLayerDist = 9999.;
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    if((fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])<0 && bestEndcapLayer>0) {
      bestEndcapLayer--;   
    }
    double leftDist=0;
    if((fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])<0 && bestEndcapLayer==0)
      leftDist = 0;
    else
      leftDist = fabs(fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer]);
    if(_algoParams->GetDebug()>2)
      cout << "Accounted Endcap layer " << bestEndcapLayer << " for COG, leaving distance: " << leftDist << endl;
    //int layer=bestEndcapLayer;
    double xPos=cluster_cog.x*cos(GAMMA_cog)-cluster_cog.y*sin(GAMMA_cog);
    double yPos=cluster_cog.y*cos(GAMMA_cog)+cluster_cog.x*sin(GAMMA_cog);
    double zPos=cluster_cog.z;
    //    float cosCOG=zPos/sqrt(yPos*yPos+zPos*zPos);
    //    if(cosCOG<_cosOfBarrel) { // ...without Barrel Projection
    bool barrel_proj = 1;
    
    if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) {
      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
	barrel_proj=0;
    }
    else {
      double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
	barrel_proj=0;
    }
    
    if(barrel_proj==0) { // Endcap hit without Barrel projection
    
      //if(!(cosCOG<_cosOfBarrel && cosCOG>-_cosOfBarrel)) { 
      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/cluster_cog.z);
      vec3 frontplaneProj;
      frontplaneProj.x=cluster_cog.x*z_ratio;
      frontplaneProj.y=cluster_cog.y*z_ratio;
      frontplaneProj.z=cluster_cog.z*z_ratio;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
      double abs_th=0;
      int lay_i=bestEndcapLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftDist;
      double X0_ratio = 0;
      if(bestEndcapLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])-_geomParams->Get_zOfEndcap());
      else
	X0_ratio = abs_th;
      X0_cog=(X0_ratio*dist)/3.5;
      if(_algoParams->GetDebug()>2)
	cout << "Passed " << bestEndcapLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
    }
    else { //...with Barrel Projection
      vec3 entryBarrel;   //  first get material passed in barrel
      entryBarrel.y=_geomParams->Get_rOfBarrel();
      entryBarrel.x=xPos*entryBarrel.y/yPos;
      entryBarrel.z=zPos*entryBarrel.y/yPos;
      vec3 exitBarrel;
      if(zPos>0)
	exitBarrel.z=_geomParams->Get_zOfBarrel();
      else
	exitBarrel.z=-_geomParams->Get_zOfBarrel();
      exitBarrel.x=xPos*exitBarrel.z/zPos;
      exitBarrel.y=yPos*exitBarrel.z/zPos;
      double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
      int   bestBarrelLayer=0;
      double barrelLayerDist = 9999.;
      double prodRadiusMax=0.;
      for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
	double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
	if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
      }
      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	  barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	  bestBarrelLayer = ilayer;
	}
      }

      if((prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])<0 && bestBarrelLayer>0) {
	bestBarrelLayer--;   
      }
      double leftBarrelDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];

      double abs_th=0;
      int lay_i=bestBarrelLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftBarrelDist;
      double X0_ratio = 0;
      if(bestBarrelLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
      else
	X0_ratio = abs_th;
      double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
      if(_algoParams->GetDebug()>2)
	cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 ..." << endl;
      // then adding "normal"endcap material passed
      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/cluster_cog.z);
      vec3 frontplaneProj;
      frontplaneProj.x=cluster_cog.x*z_ratio;
      frontplaneProj.y=cluster_cog.y*z_ratio;
      frontplaneProj.z=cluster_cog.z*z_ratio;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
      abs_th=0;
      lay_i=bestEndcapLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftDist;
      X0_ratio=abs_th/(fabs(cluster_cog.z)-_geomParams->Get_zOfEndcap());
      if(_algoParams->GetDebug()>2)
	cout << "...and " << bestEndcapLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 ..." << endl;
      X0_cog=((X0_ratio*dist)/3.5)+barrelX0Passed;
    }
  }
  clusPar->depth=X0_cog;
  //cout << "Depth: " << clusPar->depth << endl;

  int nTracks = tracks.size();
  
  for(int t_i=0;t_i<nTracks;t_i++) {
    ExtendedTrack *a_track = dynamic_cast<ExtendedTrack*> (tracks[t_i]);
    HelixClass *a_helix = a_track->helix;
    float dist[3];
    float cl_cog[3];
    cl_cog[0]=clusPar->COGx;
    cl_cog[1]=clusPar->COGy;
    cl_cog[2]=clusPar->COGz;
    float dummy;
    dummy = a_helix->getDistanceToPoint(cl_cog,dist);
    if(dist[2]<dist_track)
      dist_track=dist[2];
    for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
      float hitPos[3];
      hitPos[0]=myHit->hit->getPosition()[0];
      hitPos[1]=myHit->hit->getPosition()[1];
      hitPos[2]=myHit->hit->getPosition()[2];
      dummy = a_helix->getDistanceToPoint(hitPos,dist);
      if(dist[2]<shortest_dist_track)
	shortest_dist_track=dist[2];
    }
  }
  clusPar->distToTrack=dist_track;
  clusPar->smallestDistToTrack=shortest_dist_track;

  int nPreshowerHits=0;
  int nPL0Hits=0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->preShower)
      nPreshowerHits++;
    if(myHit->pseudoLayer==0)
      nPL0Hits++;
  }
  clusPar->nHits=(NClusteredHits+NCountedGhostHits);//-nPreshowerHits);
  clusPar->psHits=nPreshowerHits;
  clusPar->pl0Hits=nPL0Hits;
  //postion correction via phi and theta

  double clus_phi = clusPar->phi;
  int a_gamma;
  if(clus_phi<0) 
    clus_phi=twopi+clus_phi;
  a_gamma=(int)(((clus_phi-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
  double GAMMA=a_gamma*(twopi/_geomParams->Get_symmetry());
  clus_phi = fabs(clus_phi-GAMMA);
      
  if(clusPar->zone==1) {  // barrel
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    clusPar->Etot_g_theta = clusPar->Etot_g*sin(clusPar->theta);
    clusPar->Etot_g_pos = clusPar->Etot_g*sin(clusPar->theta)*cos(clus_phi);
  }

  if(clusPar->zone==2) {  // endcap
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    clusPar->Etot_g_theta = clusPar->Etot_g*fabs(cos(clusPar->theta));
    clusPar->Etot_g_pos = clusPar->Etot_g*fabs(cos(clusPar->theta))*cos(clus_phi);
  }

  if(clusPar->zone==3) {  // overlap: decide from position of cog 
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    if(cluster_cog.z<_geomParams->Get_zOfBarrel() && cluster_cog.z>-_geomParams->Get_zOfBarrel()) { //treat as in barrel
      clusPar->Etot_g_theta = clusPar->Etot_g*sin(clusPar->theta);
      clusPar->Etot_g_pos = clusPar->Etot_g*sin(clusPar->theta)*cos(clus_phi);
    }
    else { // treat as in endcap
      clusPar->Etot_g_theta = clusPar->Etot_g*fabs(cos(clusPar->theta));
      clusPar->Etot_g_pos = clusPar->Etot_g*fabs(cos(clusPar->theta))*cos(clus_phi);
    }
  }

//  if(clusPar->zone==1) {
//    //    clusPar->E_GeV_pre = toGeVFctn->Eval(clusPar->Etot_g);
//    clusPar->E_GeV_pre = ECALGarlicEnergyEstimator::Instance().getGeV(clusPar->Etot_g);
//  } else {
//    //clusPar->E_GeV_pre = toGeVFctn_EC->Eval(clusPar->Etot_g);
//    clusPar->E_GeV_pre = ECALGarlicEnergyEstimator::Instance().getGeV_EC(clusPar->Etot_g);
//  }

  bool isbarrel = clusPar->zone==1;
  clusPar->E_GeV_pre = _energyEstimator->getGeV(clusPar->Etot_g, isbarrel);


  
  clusPar->E_GeV_noLC = clusPar->E_GeV; 
  clusPar->Etot_g_noLC = clusPar->Etot_g;
  
  //  if(_optimiseResolution) {
  //   double estimated_energy = EstimateEnergy(clusPar,clusPar->E_GeV_pre );
  //  clusPar->E_GeV = estimated_energy;
  // }
  //else
    clusPar->E_GeV = clusPar->E_GeV_pre ;
    clusPar->E_GeV_noTheta = clusPar->E_GeV;
    clusPar->E_GeV_noPhi = clusPar->E_GeV;
  if(_algoParams->GetCorrectTheta()) {
    if(clusPar->zone==1 || clusPar->zone==3) {

      //      double corrected_energy = CorrectTheta(clusPar);

      double corrected_energy = _energyEstimator->correctEnergyTheta(clusPar->E_GeV, clusPar->theta);

      clusPar->E_GeV_noTheta = clusPar->E_GeV;
      clusPar->E_GeV = corrected_energy ;
    }
  }
  
  if(_algoParams->GetCorrectPhi()) {
    if(clusPar->zone==1 || clusPar->zone==3) {
      //      double corrected_energy = CorrectPhi(clusPar);

      double corrected_energy = _energyEstimator->correctEnergyPhi(clusPar->E_GeV, clusPar->phi);


      clusPar->E_GeV_noPhi = clusPar->E_GeV;
      clusPar->E_GeV = corrected_energy ;
    }
  }
  /*
  CalculateSeedDirection(myCluster);
  CalculateClusterDirection(myCluster);
  FitEllipsoid(myCluster);
  
  clusPar->Density = ((myCluster->Shape)->Volume)/(clusPar->nHits);

  float r_phi_dir=sqrt((myCluster->dir).x*(myCluster->dir).x+(myCluster->dir).y*(myCluster->dir).y);
  float theta_dir = atan2(r_phi_dir,(myCluster->dir).z);
  if(theta_dir<0) 
    theta_dir=twopi+theta_dir;
  float phi_dir = atan2((myCluster->dir).x,(myCluster->dir).y);
  if(phi_dir<0) 
    phi_dir=twopi+phi_dir;
  float dir_err = sqrt(pow(phi_dir-clus_phi,2)+pow((theta_dir-(clusPar->theta)),2));
  clusPar->dirErr = dir_err;

  r_phi_dir=sqrt((myCluster->seed_dir).x*(myCluster->seed_dir).x+(myCluster->seed_dir).y*(myCluster->seed_dir).y);
  theta_dir = atan2(r_phi_dir,(myCluster->seed_dir).z);
  if(theta_dir<0) 
    theta_dir=twopi+theta_dir;
  phi_dir = atan2((myCluster->seed_dir).x,(myCluster->seed_dir).y);
  if(phi_dir<0) 
    phi_dir=twopi+phi_dir;
  float seed_dir_err = sqrt(pow(phi_dir-clus_phi,2)+pow((theta_dir-(clusPar->theta)),2));
  clusPar->seedDirErr = seed_dir_err;
  */

  // Fit longitudinal shower shape
  //float b=0.5;
  //float Cj=0.5;   // Cj=-0.5, Cj=+0.5 for photons
  //float Z=74;
  //float Ec=0.61/(Z+1.24); 
  //float a=b*(TMath::Log(clusPar->Etot/Ec)-Cj)+1;
  //float a=long_sh_profile->GetBinCenter((long_sh_profile->GetMaximumBin()))*b+1;
  //TF1 *long_sh_shape = new TF1("longShSh","[0]*[1]*(pow([1]*(x),([2]-1))*exp(-1*[1]*(x)))/TMath::Gamma([2])",0,30);
  //long_sh_shape->SetParameter(0,Etot);
  //long_sh_shape->FixParameter(1,b);
  //long_sh_shape->SetParameter(2,a);
  //long_sh_profile->Draw();
  //long_sh_profile->Fit(long_sh_shape,"QN0");
  //  long_sh_profile->Fit(long_sh_shape);
  float chi2_long = 0;
  //if(long_sh_shape) {
  // float chi2=long_sh_shape->GetChisquare();
    //float ndf=long_sh_shape->GetNDF();
  // chi2_long=chi2;//ndf;
    
  //}
  clusPar->Chi2_long=chi2_long;
  if(_algoParams->GetDebug()>2) {
    cout << "Cluster Paramerers of cluster " << clus_ID << " : " << endl
	 << "ID: " << clusPar->ID << endl
	 << "Start: " << clusPar->start << endl
	 << "End: " << clusPar->end << endl
	 << "Etot: " << clusPar->Etot << endl
	 << "Etot_g: " << clusPar->Etot_g << endl
	 << "Etot_g_pos: " << clusPar->Etot_g_pos << endl
	 << "E_GeV: " << clusPar->E_GeV << endl
      //<< "Es1ovEtot: " << clusPar->Es1ovEtot << endl
	 << "E1C: " << clusPar->E1C << endl
	 << "E4C: " << clusPar->E4C << endl
	 << "NHits: " << clusPar->nHits << endl
	 << "Theta: " << clusPar->theta << endl
	 << "Phi: " << clusPar->phi << endl
	 << "<depth>: " << clusPar->depth << endl
	 << "COG at: " << clusPar->COGx << ", " << clusPar->COGy << ", " << clusPar->COGz << endl
	 << "Distance to biggest: " << clusPar->distToBiggest << endl
	 << "Smallest Distance to biggest: " << clusPar->smallestDistToBiggest << endl
	 << "Distance to track: " << clusPar->distToTrack << endl
	 << "Smallest Distance to track: " << clusPar->smallestDistToTrack << endl
	 << "Chi2_long: " << clusPar->Chi2_long << endl 
	 << "Desnsity: " << clusPar->hitDensity << endl 
	 << "Direction Error: " << clusPar->dirErr << endl 
	 << "Seed Direction Error: " << clusPar->seedDirErr << endl 
	 << "Zone: " << clusPar->zone << endl;
  }

  //tree_file->Write();
  //delete long_sh_profile;
  //delete long_sh_shape;
}



void ECALGarlicClusterHelpers::CalculatePhotonProbability(LCEvent *evt, 
							  ExtendedCluster *myCluster, 
							  vector<ExtendedCluster*> *clusters, 
							  int clus_ID, 
							  ClusterParameters *clusPar,
							  vector<ExtendedTrack*> tracks, 
							  map<ExtendedTrack*,map<int,vector<ExtendedHit*> > > *allRemovedHits)
{


  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  clusPar->ID=clus_ID;
  LCCollection *preClusColl = 0; 
  preClusColl = evt->getCollection(_algoParams->GetEcalPreClusterCollectionName() );
  CellIDDecoder<CalorimeterHit> decoder(preClusColl);  
  
  // 1.) energy distribution in structures
  double en_s1=0;
  double en_s2=0;
  double en_s3=0;
  double enPs=0;
  double en1odd=0;
  double en1even=0;
  double en2odd=0;
  double en2even=0;
  double nPs=0;
  double n1odd=0;
  double n1even=0;
  double n2odd=0;
  double n2even=0;
  // 2.) cluster center of gravity
  vec3 cluster_cog;
  cluster_cog.x=0;
  cluster_cog.y=0;
  cluster_cog.z=0;
  double en_norm=0;
  // 3.) longitudinal shower profile
  //TH1F *long_sh_profile=new TH1F("LongitudinalShowerProfile","LongitudinalShowerProfile",80,0,40);
  // 4.) start/end shower
  double min_X0=50;
  double max_X0=0;
  // 5.) E4C/Etot
  //  vec3 *seed_dir=&(myCluster->seededFrom);
  // 6.) E1C/Etot
  double E1C_en=0;
  double E4C_en=0;
  double E9C_en=0;
  // 7.) nearest track
  double dist_track=9999;
  double shortest_dist_track=9999;
  // 8.) Etot/Multiplicity
  // 9.) Mean Shower Depth
  // 10.) Cluster Zone
  clusPar->zone=0;

  vec3 firstHitPos;
  bool hasBarrelHits=0;
  bool hasEndcapHits=0;
  vector<ExtendedHit* > *clusterHits = &(myCluster->hitVec);
  int NClusteredHits = clusterHits->size();
  // calculate COG
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    en_norm+=hit_en;
    cluster_cog.x+=(a_hit->getPosition()[0])*hit_en;
    cluster_cog.y+=(a_hit->getPosition()[1])*hit_en;
    cluster_cog.z+=(a_hit->getPosition()[2])*hit_en;
   int module = decoder(a_hit)["M"];
   if(0 < module && module < 6)
     hasBarrelHits=1;
   else
     hasEndcapHits=1;
  }

  if(myCluster->Ghosts!=0) {
    GhostCluster *ghostHits = myCluster->Ghosts;
    if(ghostHits!=0) {
      vector<GhostHit*> ghostHitVec = (ghostHits->ghostHitVec);
      int NGhostHits = ghostHitVec.size();
      if(NGhostHits>0) {
	//	cout << "Ghost cluster hits: " << NGhostHits << endl;
	for(int g_i=0;g_i<NGhostHits;g_i++) {
	  GhostHit *g_hit= dynamic_cast<GhostHit*>(ghostHitVec[g_i]);
	  float hit_en=g_hit->Energy;
	  cluster_cog.x+=(g_hit->Position).x*hit_en;
	  cluster_cog.y+=(g_hit->Position).y*hit_en;
	  cluster_cog.z+=(g_hit->Position).z*hit_en;
	  en_norm+=hit_en;
	}
      }
    }
  }
  cluster_cog.x=cluster_cog.x/en_norm;
  cluster_cog.y=cluster_cog.y/en_norm;
  cluster_cog.z=cluster_cog.z/en_norm;
  clusPar->COGx=cluster_cog.x;
  clusPar->COGy=cluster_cog.y;
  clusPar->COGz=cluster_cog.z;
  if(_algoParams->GetDebug()>2)
    cout << "Cluster COG is at : " << cluster_cog.x << ", " << cluster_cog.y << ", " << cluster_cog.z << endl;
  if(hasBarrelHits==1 && hasEndcapHits==0)
    clusPar->zone=1;
  if(hasBarrelHits==0 && hasEndcapHits==1)
    clusPar->zone=2;
  if(hasBarrelHits==1 && hasEndcapHits==1)
    clusPar->zone=3;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    //if(myHit->preShower==1)
    //  continue;
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    vec3 hitPos;
    hitPos.x=a_hit->getPosition()[0];
    hitPos.y=a_hit->getPosition()[1];
    hitPos.z=a_hit->getPosition()[2];
    //int cell_i = decoder(a_hit)["I"];
    //int cell_j = decoder(a_hit)["J"];
    int layer = decoder(a_hit)["K-1"];
    int psLayer = myHit->pseudoLayer;
	
    if(myHit->preShower==1) {
      enPs += hit_en;
      nPs++;
    }
    else {
      if(layer<20) {
      	if(layer%2 == 0) {
	  en1even += hit_en;
	  n1even++;
	}
	else {
	  en1odd += hit_en;
	  n1odd++;
	}
      }
      else {
	if(layer%2 == 0) {
	  en2even += hit_en;
	  n2even++;
	}
	else {
	  en2odd += hit_en;
	  n2odd++;
	}
      }
    }

    double lambda=(cluster_cog.x*hitPos.x)+(cluster_cog.y*hitPos.y)+(cluster_cog.z*hitPos.z);// calculate distance to line from IP through seed: "direction"
    double lambda_prime=fabs(lambda/((cluster_cog.x*cluster_cog.x)+(cluster_cog.y*cluster_cog.y)+(cluster_cog.z*cluster_cog.z)));
    vec3 lotpunkt;
    lotpunkt.x=cluster_cog.x*lambda_prime;
    lotpunkt.y=cluster_cog.y*lambda_prime;
    lotpunkt.z=cluster_cog.z*lambda_prime;
    double a_dist=_geomHelper->Get3dDistance(&lotpunkt,&hitPos);
    if(a_dist<(_geomParams->Get_padSizeEcal() [1]/2))
      E1C_en+=hit_en;
    if(a_dist<(sqrt(2.)*_geomParams->Get_padSizeEcal() [1]))
      E4C_en+=hit_en;
    if(a_dist<(sqrt(8.)*_geomParams->Get_padSizeEcal() [1]))
      E9C_en+=hit_en;
    double X0_passed=0; // calculate X0 in front of hit
    if(hitPos.z<_geomParams->Get_zOfBarrel() && hitPos.z>-_geomParams->Get_zOfBarrel()) { //Barrel Hit
      int stave = decoder(a_hit)["S-1"];
      //float stave_angle=acos(_geomParams->Get_barrelStaveDir()[stave].x);
      double stave_x=_geomParams->Get_barrelStaveDir()[stave].x;
      double stave_y=_geomParams->Get_barrelStaveDir()[stave].y;
      double r_frontplane=_geomParams->Get_rOfBarrel();
      double lambda=r_frontplane*r_frontplane;
      if(_algoParams->GetDebug()>2) {
	cout << "Hit at " << hitPos.x << ", " <<  hitPos.y << ", " << hitPos.z << endl;
	//cout << "Stave angle: " << stave_angle << ", r_frontplane: " << r_frontplane << endl;
	cout << "Stave: " << stave << " , r_frontplane* stave.x/.y: " << r_frontplane*stave_x << ", " << r_frontplane*stave_y << endl;
      }
      double lambda_prime=fabs(lambda/(hitPos.x*r_frontplane*stave_x+hitPos.y*r_frontplane*stave_y));
      if(_algoParams->GetDebug()>2)
	cout << "lambda_prime: " << lambda_prime << endl;
      vec3 frontplaneProj;
      frontplaneProj.x=hitPos.x*lambda_prime;
      frontplaneProj.y=hitPos.y*lambda_prime;
      frontplaneProj.z=hitPos.z*lambda_prime;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
      if(_algoParams->GetDebug()>2)
	cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
      double abs_th=0;
      int lay_i=psLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	lay_i--;
      }
      double X0_ratio = 0;
      if(psLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[psLayer])-_geomParams->Get_rOfBarrel());
      else 
	X0_ratio = 0;
      X0_passed=(X0_ratio*dist)/3.5;
      if(_algoParams->GetDebug()>2)
	cout << "X0 passed: " << X0_passed << endl;
    }
    else { //Endcap hit...
      double phi_pos = atan2((double)hitPos.x,(double)hitPos.y);
      int gamma;
      if(phi_pos<0) 
	phi_pos=twopi+phi_pos;
      gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
      double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
      double xPos=hitPos.x*cos(GAMMA)-hitPos.y*sin(GAMMA);
      double yPos=hitPos.y*cos(GAMMA)+hitPos.x*sin(GAMMA);
      double zPos=hitPos.z;
      // float cosHit=zPos/sqrt(yPos*yPos+zPos*zPos);
      //if(!(cosHit<_cosOfBarrel && cosHit>-_cosOfBarrel && clusPar->zone==2)) { 
      // if(cosHit<_cosOfBarrel) { // ...without Barrel Projection and cluster in overlap zone
      bool barrel_proj = 1;
      
      if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) {
	if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
	  barrel_proj=0;
      }
      else {
	double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
	if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
	  barrel_proj=0;
      }
      
      if(barrel_proj==0 && clusPar->zone==2) { // Endcap hit without Barrel projection
    
	double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	//	if(z_ratio>1)
	//	  cout << "Warning: z_ratio >1: " << z_ratio << endl;
	vec3 frontplaneProj;
	frontplaneProj.x=hitPos.x*z_ratio;
	frontplaneProj.y=hitPos.y*z_ratio;
	frontplaneProj.z=hitPos.z*z_ratio;
	double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	double abs_th=0;
	int lay_i=psLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	  lay_i--;
	}
	double X0_ratio = 0;
	if(psLayer>0)
	  X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	else
	  X0_ratio = 0;
	X0_passed=(X0_ratio*dist)/3.5;
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	if(_algoParams->GetDebug()>2)
	  cout << "X0 passed: " << X0_passed << endl;
      }
      else { //...with Barrel Projection
	if(_algoParams->GetDebug()>2)
	  cout << "hit with barrel projection" << endl;
	vec3 entryBarrel;   //  first get material passed in barrel
	entryBarrel.y=_geomParams->Get_rOfBarrel();
	entryBarrel.x=xPos*entryBarrel.y/yPos;
	entryBarrel.z=zPos*entryBarrel.y/yPos;
	vec3 exitBarrel;
	if(zPos>0)
	  exitBarrel.z=_geomParams->Get_zOfBarrel();
	else
	  exitBarrel.z=-_geomParams->Get_zOfBarrel();
	exitBarrel.x=xPos*exitBarrel.z/zPos;
	exitBarrel.y=yPos*exitBarrel.z/zPos;
	double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
	int   bestBarrelLayer=0;
	double barrelLayerDist = 9999.;
	double prodRadiusMax=0.;
	for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
	  double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
	  if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
	}
	for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	  if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	    barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	    bestBarrelLayer = ilayer;
	  }
	}
	if(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer]<0 && bestBarrelLayer>0) {
	  bestBarrelLayer--;   
	}
	double leftBarrelDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
	
	double abs_th=0;
	int lay_i=bestBarrelLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	  lay_i--;
	}
	abs_th+=leftBarrelDist;
	double X0_ratio = 0;
	if(bestBarrelLayer>0)
	  X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
	else
	  X0_ratio = 0;
	double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	// then adding "normal"endcap material passed
	double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	vec3 frontplaneProj;
	frontplaneProj.x=hitPos.x*z_ratio;
	frontplaneProj.y=hitPos.y*z_ratio;
	frontplaneProj.z=hitPos.z*z_ratio;
	double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	abs_th=0;
	lay_i=psLayer;
	while(lay_i>0) {
	  abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	  lay_i--;
	}
	if(psLayer>0)
	  X0_ratio=abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	else 
	  X0_ratio = 0;
	if(_algoParams->GetDebug()>2)
	  cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	X0_passed=((X0_ratio*dist)/3.5)+barrelX0Passed;
	if(_algoParams->GetDebug()>2)
	  cout << "X0 passed: " << X0_passed << endl;
      }
    }

    //cout << "X0 passed: " << X0_passed << endl;
    if(X0_passed>max_X0)
      max_X0=X0_passed;
    if(X0_passed<min_X0) {
      min_X0=fabs(X0_passed);
      firstHitPos=hitPos;
    }
    //    if(X0_passed>30)
    //  cout << "Still sth wrong in X0 passed, " << X0_passed << endl;
    //long_sh_profile->Fill(X0_passed,hit_en);
    if(layer<10)
      en_s1+=hit_en;
    else {
      if(layer<25)
	en_s2+=hit_en;
      else
	en_s3+=hit_en;
    }
  }

  double Etot=en_s1+en_s2+en_s3;

  if(myCluster->Ghosts!=0) {
    GhostCluster *ghostHits = myCluster->Ghosts;
    if(ghostHits!=0) {
      vector<GhostHit*> ghostHitVec = (ghostHits->ghostHitVec);
      int NGhostHits = ghostHitVec.size();
      if(NGhostHits>0) {
	//	cout << "Ghost cluster hits: " << NGhostHits << endl;
	for(int g_i=0;g_i<NGhostHits;g_i++) {
	  GhostHit *g_hit= dynamic_cast<GhostHit*>(ghostHitVec[g_i]);
	  float hit_en=g_hit->Energy;
	  int layer = g_hit->Layer;
	  int psLayer = g_hit->PseudoLayer;
	  vec3 hitPos;
	  hitPos.x=(g_hit->Position).x;
	  hitPos.y=(g_hit->Position).y;
	  hitPos.z=(g_hit->Position).z;
	  double lambda=(cluster_cog.x*hitPos.x)+(cluster_cog.y*hitPos.y)+(cluster_cog.z*hitPos.z);// calculate distance to line from IP through seed: "direction"
	  double lambda_prime=fabs(lambda/((cluster_cog.x*cluster_cog.x)+(cluster_cog.y*cluster_cog.y)+(cluster_cog.z*cluster_cog.z)));
	  vec3 lotpunkt;
	  lotpunkt.x=cluster_cog.x*lambda_prime;
	  lotpunkt.y=cluster_cog.y*lambda_prime;
	  lotpunkt.z=cluster_cog.z*lambda_prime;
	  double a_dist=_geomHelper->Get3dDistance(&lotpunkt,&hitPos);
	  if(a_dist<(_geomParams->Get_padSizeEcal() [1]/2))
	    E1C_en+=hit_en;
	  if(a_dist<(sqrt(2.)*_geomParams->Get_padSizeEcal() [1]))
	    E4C_en+=hit_en;
	  if(a_dist<(sqrt(8.)*_geomParams->Get_padSizeEcal() [1]))
	    E9C_en+=hit_en;
	  double X0_passed=0; // calculate X0 in front of hit
	  if(hitPos.z<_geomParams->Get_zOfBarrel() && hitPos.z>-_geomParams->Get_zOfBarrel()) { //Barrel Hit
	    int stave = g_hit->Stave;
	    double stave_x=_geomParams->Get_barrelStaveDir()[stave].x;
	    double stave_y=_geomParams->Get_barrelStaveDir()[stave].y;
	    double r_frontplane=_geomParams->Get_rOfBarrel();
	    double lambda=r_frontplane*r_frontplane;
	    if(_algoParams->GetDebug()>2) {
	      cout << "Hit at " << hitPos.x << ", " <<  hitPos.y << ", " << hitPos.z << endl;
	      //cout << "Stave angle: " << stave_angle << ", r_frontplane: " << r_frontplane << endl;
	      cout << "Stave: " << stave << " , r_frontplane* stave.x/.y: " << r_frontplane*stave_x << ", " << r_frontplane*stave_y << endl;
	    }
	    double lambda_prime=fabs(lambda/(hitPos.x*r_frontplane*stave_x+hitPos.y*r_frontplane*stave_y));
	    if(_algoParams->GetDebug()>2)
	      cout << "lambda_prime: " << lambda_prime << endl;
	    vec3 frontplaneProj;
	    frontplaneProj.x=hitPos.x*lambda_prime;
	    frontplaneProj.y=hitPos.y*lambda_prime;
	    frontplaneProj.z=hitPos.z*lambda_prime;
	    double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	    if(_algoParams->GetDebug()>2)
	      cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
	    double abs_th=0;
	    int lay_i=psLayer;
	    while(lay_i>0) {
	      abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
	      lay_i--;
	    }
	    double X0_ratio = 0;
	    if(psLayer>0)
	      X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[psLayer])-_geomParams->Get_rOfBarrel());
	    else 
	      X0_ratio = 0;
	    X0_passed=(X0_ratio*dist)/3.5;
	    if(_algoParams->GetDebug()>2)
	      cout << "X0 passed: " << X0_passed << endl;
	  }
	  else { //Endcap hit...
	    double phi_pos = atan2((double)hitPos.x,(double)hitPos.y);
	    int gamma;
	    if(phi_pos<0) 
	      phi_pos=twopi+phi_pos;
	    gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
	    double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
	    double xPos=hitPos.x*cos(GAMMA)-hitPos.y*sin(GAMMA);
	    double yPos=hitPos.y*cos(GAMMA)+hitPos.x*sin(GAMMA);
	    double zPos=hitPos.z;
	    // float cosHit=zPos/sqrt(yPos*yPos+zPos*zPos);
	    //if(!(cosHit<_cosOfBarrel && cosHit>-_cosOfBarrel && clusPar->zone==2)) { 
	    // if(cosHit<_cosOfBarrel) { // ...without Barrel Projection and cluster in overlap zone
	    bool barrel_proj = 1;
	    
	    if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) {
	      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
		barrel_proj=0;
	    }
	    else {
	      double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
	      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
		barrel_proj=0;
	    }
	    
	    if(barrel_proj==0 && clusPar->zone==2) { // Endcap hit without Barrel projection
	      
	      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	      //	if(z_ratio>1)
	      //	  cout << "Warning: z_ratio >1: " << z_ratio << endl;
	      vec3 frontplaneProj;
	      frontplaneProj.x=hitPos.x*z_ratio;
	      frontplaneProj.y=hitPos.y*z_ratio;
	      frontplaneProj.z=hitPos.z*z_ratio;
	      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	      double abs_th=0;
	      int lay_i=psLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
		lay_i--;
	      }
	      double X0_ratio = 0;
	      if(psLayer>0)
		X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	      else
		X0_ratio = 0;
	      X0_passed=(X0_ratio*dist)/3.5;
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      if(_algoParams->GetDebug()>2)
		cout << "X0 passed: " << X0_passed << endl;
	    }
	    else { //...with Barrel Projection
	      if(_algoParams->GetDebug()>2)
		cout << "hit with barrel projection" << endl;
	      vec3 entryBarrel;   //  first get material passed in barrel
	      entryBarrel.y=_geomParams->Get_rOfBarrel();
	      entryBarrel.x=xPos*entryBarrel.y/yPos;
	      entryBarrel.z=zPos*entryBarrel.y/yPos;
	      vec3 exitBarrel;
	      if(zPos>0)
		exitBarrel.z=_geomParams->Get_zOfBarrel();
	      else
		exitBarrel.z=-_geomParams->Get_zOfBarrel();
	      exitBarrel.x=xPos*exitBarrel.z/zPos;
	      exitBarrel.y=yPos*exitBarrel.z/zPos;
	      double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
	      int   bestBarrelLayer=0;
	      double barrelLayerDist = 9999.;
	      double prodRadiusMax=0.;
	      for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
		double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
		if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
	      }
	      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
		if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
		  barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
		  bestBarrelLayer = ilayer;
		}
	      }
	      if(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer]<0 && bestBarrelLayer>0) {
		bestBarrelLayer--;   
	      }
	      double leftBarrelDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
	      
	      double abs_th=0;
	      int lay_i=bestBarrelLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
		lay_i--;
	      }
	      abs_th+=leftBarrelDist;
	      double X0_ratio = 0;
	      if(bestBarrelLayer>0)
		X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
	      else
		X0_ratio = 0;
	      double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      // then adding "normal"endcap material passed
	      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/hitPos.z);
	      vec3 frontplaneProj;
	      frontplaneProj.x=hitPos.x*z_ratio;
	      frontplaneProj.y=hitPos.y*z_ratio;
	      frontplaneProj.z=hitPos.z*z_ratio;
	      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&hitPos);
	      abs_th=0;
	      lay_i=psLayer;
	      while(lay_i>0) {
		abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
		lay_i--;
	      }
	      if(psLayer>0)
		X0_ratio=abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[psLayer])-_geomParams->Get_zOfEndcap());
	      else 
		X0_ratio = 0;
	      if(_algoParams->GetDebug()>2)
		cout << "Passed " << psLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
	      X0_passed=((X0_ratio*dist)/3.5)+barrelX0Passed;
	      if(_algoParams->GetDebug()>2)
		cout << "X0 passed: " << X0_passed << endl;
	    }
	  }
	  //long_sh_profile->Fill(X0_passed,hit_en);
	  if(layer<10)
	    en_s1+=hit_en;
	  else {
	    if(layer<25)
	      en_s2+=hit_en;
	    else
	      en_s3+=hit_en;
	  }
	}
      }
    }
  }

  clusPar->start=min_X0;
  clusPar->end=max_X0;
  clusPar->Etot=Etot;
  clusPar->Es1=en_s1;
  clusPar->Es2=en_s2;
  clusPar->Es3=en_s3;

  double etot_g = Etot;

  //cout << "Sums: " << Etot << ", " << en1even+en1odd+en2even+en2odd << endl;

  int NCountedGhostHits=0;
  if(myCluster->Ghosts!=0) {
    GhostCluster *ghostHits = myCluster->Ghosts;
    if(ghostHits!=0) {
      vector<GhostHit*> ghostHitVec = (ghostHits->ghostHitVec);
      int NGhostHits = ghostHitVec.size();
      if(NGhostHits>0) {
	//cout << "Ghost cluster hits: " << NGhostHits << endl;
	for(int g_i=0;g_i<NGhostHits;g_i++) {
	  GhostHit *g_hit= dynamic_cast<GhostHit*>(ghostHitVec[g_i]);
	  float g_hit_en = g_hit->Energy;
	  int layer = g_hit->Layer;
	  etot_g+=g_hit_en;
	  if(layer<20) {
	    if(layer%2 == 0) {
	      en1even += g_hit_en;	     
	      //cout << "Added Ghost hit En: " << g_hit_en << endl;
	      if(g_hit->Count==1)
		n1even++;
	    }
	    else {
	      en1odd += g_hit_en;
	      if(g_hit->Count==1)
		n1odd++;
	    }
	  }
	  else {
	    if(layer%2 == 0) {
	      en2even += g_hit_en;
	      if(g_hit->Count==1)
		n2even++;
	    }
	    else {
	      en2odd += g_hit_en;
	      if(g_hit->Count==1)
		n2odd++;
	    }
	  }
	  if(g_hit->Count==1)
	    NCountedGhostHits++;
	  //cout << etot_g << endl;
	  //cout << en1even+en1odd+en2even+en2odd << endl;
	}
      }
    }
  }

  clusPar->Etot_g = etot_g;
  clusPar->EnPs=enPs;
  clusPar->En1even=en1even;
  clusPar->En1odd=en1odd;
  clusPar->En2even=en2even;
  clusPar->En2odd=en2odd;
  clusPar->NPs=nPs;
  clusPar->N1even=n1even;
  clusPar->N1odd=n1odd;
  clusPar->N2even=n2even;
  clusPar->N2odd=n2odd;
  clusPar->nGhostHits=NCountedGhostHits;
  
  clusPar->E1C=E1C_en;
  clusPar->E4C=E4C_en;
  clusPar->E9C=E9C_en;

  double firstHitDistance=9999;
  int firstPSLayer=99;
  // calculate distance from biggest cluster in same roi
  double projCOGDist=9999;
  double smallest_dist=2*_geomParams->Get_zOfBarrel();
  if(clusPar->ID>0) {
    ExtendedCluster *myBigCluster=(*clusters)[0];
    ClusterParameters *cp0=myBigCluster->Parameters;
    vec3 another_cog;
    another_cog.x=cp0->COGx;
    another_cog.y=cp0->COGy;
    another_cog.z=cp0->COGz;
    if(_algoParams->GetDebug()>2)
      cout << "Biggest cluster COG: " << another_cog.x << ", " << another_cog.y << ", " << another_cog.z << endl
	   << "Local COG: " << cluster_cog.x << ", " <<  cluster_cog.y << ", " <<  cluster_cog.z << endl;
    double lambda=(cluster_cog.x*another_cog.x)+(cluster_cog.y*another_cog.y)+(cluster_cog.z*another_cog.z);// calculate distance to line from IP through seed: "direction"
    double lambda_prime=fabs(lambda/((another_cog.x*another_cog.x)+(another_cog.y*another_cog.y)+(another_cog.z*another_cog.z)));
    vec3 lotpunkt;
    lotpunkt.x=another_cog.x*lambda_prime;
    lotpunkt.y=another_cog.y*lambda_prime;
    lotpunkt.z=another_cog.z*lambda_prime;
    double a_dist=_geomHelper->Get3dDistance(&lotpunkt,&cluster_cog);
    projCOGDist=a_dist;
  
    vec3 firstBigClHit;
    // now get closest distance between two clusters
    vector<ExtendedHit* > *bigClusterHits = &(myBigCluster->hitVec);
    int NBigClusteredHits = bigClusterHits->size();
    for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
      vec3 a_pos;
      a_pos.x=(myHit->hit)->getPosition()[0];
      a_pos.y=(myHit->hit)->getPosition()[1];
      a_pos.z=(myHit->hit)->getPosition()[2];
      for(int big_cl_hit_i=0;big_cl_hit_i<NBigClusteredHits;big_cl_hit_i++) {
	ExtendedHit *myBigHit = dynamic_cast<ExtendedHit* > ((*bigClusterHits)[big_cl_hit_i]);
	vec3 a_big_pos;
	a_big_pos.x=(myBigHit->hit)->getPosition()[0];
	a_big_pos.y=(myBigHit->hit)->getPosition()[1];
	a_big_pos.z=(myBigHit->hit)->getPosition()[2];
	double dist=_geomHelper->Get3dDistance(&a_pos,&a_big_pos);
	if(dist<smallest_dist)
	  smallest_dist=dist;
	if(myBigHit->pseudoLayer<=firstPSLayer) {
	  firstPSLayer=myBigHit->pseudoLayer;
	  double a_firstDistance = _geomHelper->Get2dProjDistance(&a_big_pos,&firstHitPos);
	  if(a_firstDistance<firstHitDistance) {
	    firstBigClHit=a_big_pos;
	    firstHitDistance=a_firstDistance;
	  }
	}
      }
    }
  }
  clusPar->distToBiggest=projCOGDist;
  clusPar->smallestDistToBiggest=smallest_dist;
  clusPar->distFirstCell=firstHitDistance;

  double r_phi=sqrt(cluster_cog.x*cluster_cog.x+cluster_cog.y*cluster_cog.y);
  double theta_cog = atan2((double)r_phi,(double)cluster_cog.z);
  clusPar->theta=theta_cog;
  // get mean shower depth from COG: same procedure as for hits, FIXME: be more precise?
  double phi_cog = atan2((double)cluster_cog.x,(double)cluster_cog.y);
  int gamma_cog;
  if(phi_cog<0) 
    phi_cog=twopi+phi_cog;
  clusPar->phi=phi_cog;
  gamma_cog=(int)(((phi_cog-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
  double GAMMA_cog=gamma_cog*(twopi/_geomParams->Get_symmetry());
  double X0_cog=0;
  vec3 frontplaneProj;
  frontplaneProj.x = 0;
  frontplaneProj.y = 0;
  frontplaneProj.z = 0;
  if(cluster_cog.z>-_geomParams->Get_zOfBarrel() && cluster_cog.z<_geomParams->Get_zOfBarrel()) { // in Barrel
    // find best barrel layer
    int   bestBarrelLayer=0;
    double barrelLayerDist = 9999.;
    double prodRadiusMax=0.;
    for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
      double prodRadius = cluster_cog.x*_geomParams->Get_barrelStaveDir()[istave].x+cluster_cog.y*_geomParams->Get_barrelStaveDir()[istave].y;
      if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
    }
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	bestBarrelLayer = ilayer;
      }
    }
    if((prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])<0 && bestBarrelLayer>0) {
      bestBarrelLayer--;   
    }
    double leftDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
    if(_algoParams->GetDebug()>2)
      cout << "Accounted layer " << bestBarrelLayer << " for COG, leaving distance: " << leftDist << endl;
    int layer=bestBarrelLayer;
    double r_frontplane=_geomParams->Get_rOfBarrel();
    double lambda=r_frontplane*r_frontplane;
    if(_algoParams->GetDebug()>2) {
      cout << "GAMMA_COG: " << GAMMA_cog << ", r_frontplane: " << r_frontplane << endl;
      cout << "r_frontplane* (sin,cos): " << r_frontplane*sin(GAMMA_cog) << ", " << r_frontplane*cos(GAMMA_cog) << endl;
    }
    double lambda_prime=fabs(lambda/(cluster_cog.x*r_frontplane*sin(GAMMA_cog)+cluster_cog.y*r_frontplane*cos(GAMMA_cog)));
    if(_algoParams->GetDebug()>2)
      cout << "lambda_prime: " << lambda_prime << endl;
    frontplaneProj.x=cluster_cog.x*lambda_prime;
    frontplaneProj.y=cluster_cog.y*lambda_prime;
    frontplaneProj.z=cluster_cog.z*lambda_prime;
    double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
    if(_algoParams->GetDebug()>2)
      cout << "Front plane projection at " << frontplaneProj.x << ", " << frontplaneProj.y << ", " << frontplaneProj.z << " ,thus distance = " << dist  << endl;
    double abs_th=0;
    int lay_i=layer;
    while(lay_i>0) {
      abs_th+=_geomParams->Get_absThicknessBarrelLayer()[lay_i];
      lay_i--;
    }
    abs_th+=leftDist;
    double X0_ratio=0;
    if(bestBarrelLayer>0)
      X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[layer])-_geomParams->Get_rOfBarrel());
    else
      X0_ratio = abs_th; 
    X0_cog=(X0_ratio*dist)/3.5;
    if(_algoParams->GetDebug()>2)
      cout << "Passed " << layer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
  }
  else { //in Endcap ...
    // determine ENDCAP layer
    int   bestEndcapLayer=0;
    double endcapLayerDist = 9999.;
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    if((fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])<0 && bestEndcapLayer>0) {
      bestEndcapLayer--;   
    }
    double leftDist = 0;
    if((fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])<0 && bestEndcapLayer==0)
      leftDist = 0;
    else
      leftDist = fabs(cluster_cog.z)-_geomParams->Get_positionEndcapLayer()[bestEndcapLayer];
    if(_algoParams->GetDebug()>2)
      cout << "Accounted Endcap layer " << bestEndcapLayer << " for COG, leaving distance: " << leftDist << endl;
    double xPos=cluster_cog.x*cos(GAMMA_cog)-cluster_cog.y*sin(GAMMA_cog);
    double yPos=cluster_cog.y*cos(GAMMA_cog)+cluster_cog.x*sin(GAMMA_cog);
    double zPos=cluster_cog.z;
    //    float cosCOG=zPos/sqrt(yPos*yPos+zPos*zPos);
    //    if(cosCOG<_cosOfBarrel) { // ...without Barrel Projection
    //if(!(cosCOG<_cosOfBarrel && cosCOG>-_cosOfBarrel)) { 
    bool barrel_proj = 1;
    
    if((fabs(xPos)*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) {
      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < _geomParams->Get_rOfBarrel())
	barrel_proj=0;
    }
    else {
      double y_offset = ((fabs(xPos) - (_geomParams->Get_rOfBarrel()*tan(twopi/(2*_geomParams->Get_symmetry())))) * tan(twopi/_geomParams->Get_symmetry()) );
      if((yPos*_geomParams->Get_zOfBarrel()/zPos) < (_geomParams->Get_rOfBarrel()-y_offset))
	barrel_proj=0;
    }
    
    if(barrel_proj==0) { // Endcap hit without Barrel projection
    
      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/cluster_cog.z);
      frontplaneProj.x=cluster_cog.x*z_ratio;
      frontplaneProj.y=cluster_cog.y*z_ratio;
      frontplaneProj.z=cluster_cog.z*z_ratio;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
      double abs_th=0;
      int lay_i=bestEndcapLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftDist;
      double X0_ratio = 0;
      if(bestEndcapLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionEndcapLayer()[bestEndcapLayer])-_geomParams->Get_zOfEndcap());
      else
	X0_ratio = abs_th;
      X0_cog=(X0_ratio*dist)/3.5;
      if(_algoParams->GetDebug()>2)
	cout << "Passed " << bestEndcapLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 " << endl;
    }
    else { //...with Barrel Projection
      vec3 entryBarrel;   //  first get material passed in barrel
      entryBarrel.y=_geomParams->Get_rOfBarrel();
      entryBarrel.x=xPos*entryBarrel.y/yPos;
      entryBarrel.z=zPos*entryBarrel.y/yPos;
      vec3 exitBarrel;
      if(zPos>0)
	exitBarrel.z=_geomParams->Get_zOfBarrel();
      else
	exitBarrel.z=-_geomParams->Get_zOfBarrel();
      exitBarrel.x=xPos*exitBarrel.z/zPos;
      exitBarrel.y=yPos*exitBarrel.z/zPos;
      double barrelMaterialPassed=_geomHelper->Get3dDistance(&entryBarrel,&exitBarrel);
      int   bestBarrelLayer=0;
      double barrelLayerDist = 9999.;
      double prodRadiusMax=0.;
      for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
	double prodRadius = exitBarrel.x*_geomParams->Get_barrelStaveDir()[istave].x+exitBarrel.y*_geomParams->Get_barrelStaveDir()[istave].y;
	if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
      }
      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	  barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	  bestBarrelLayer = ilayer;
	}
      }
      
      if((prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])<0 && bestBarrelLayer>0) {
	bestBarrelLayer--;   
      }
      double leftBarrelDist = prodRadiusMax-_geomParams->Get_positionBarrelLayer()[bestBarrelLayer];
      
      double abs_th=0;
      int lay_i=bestBarrelLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftBarrelDist;
      if(_algoParams->GetDebug()>2)
	cout << "From projection: passed " << bestBarrelLayer << " barrel layers, while leaving distance: " << leftBarrelDist  << endl;
      double X0_ratio = 0;
      if(bestBarrelLayer>0)
	X0_ratio = abs_th/(fabs(_geomParams->Get_positionBarrelLayer()[bestBarrelLayer])-_geomParams->Get_rOfBarrel());
      else
	X0_ratio = abs_th;
      double barrelX0Passed=(X0_ratio*barrelMaterialPassed)/3.5;  
      if(_algoParams->GetDebug()>2)
	cout << "Passed " << bestBarrelLayer << " Barrel layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 ..." << endl;
      // then adding "normal"endcap material passed
      double z_ratio=fabs(_geomParams->Get_zOfEndcap()/cluster_cog.z);
      frontplaneProj.x=cluster_cog.x*z_ratio;
      frontplaneProj.y=cluster_cog.y*z_ratio;
      frontplaneProj.z=cluster_cog.z*z_ratio;
      double dist=_geomHelper->Get3dDistance(&frontplaneProj,&cluster_cog);
      abs_th=0;
      lay_i=bestEndcapLayer;
      while(lay_i>0) {
	abs_th+=_geomParams->Get_absThicknessEndcapLayer()[lay_i];
	lay_i--;
      }
      abs_th+=leftDist;
      X0_ratio=abs_th/(fabs(cluster_cog.z)-_geomParams->Get_zOfEndcap());
      if(_algoParams->GetDebug()>2)
	cout << "...and " << bestEndcapLayer << " Endcap layers, i.e. " << abs_th << " mm = " << abs_th/3.5 << " X0 ..." << endl;
      X0_cog=((X0_ratio*dist)/3.5)+barrelX0Passed;
    }
  }
  clusPar->depth=X0_cog;
  clusPar->POSx = frontplaneProj.x;
  clusPar->POSy = frontplaneProj.y;
  clusPar->POSz = frontplaneProj.z;
  
  ExtendedTrack *nearestTrack = 0;

  int nTracks = tracks.size();
  for(int t_i=0;t_i<nTracks;t_i++) {
    ExtendedTrack *a_track = dynamic_cast<ExtendedTrack*> (tracks[t_i]);
    HelixClass *a_helix = a_track->helix;
    float dist[3];
    float cl_cog[3];
    cl_cog[0]=clusPar->COGx;
    cl_cog[1]=clusPar->COGy;
    cl_cog[2]=clusPar->COGz;
    float dummy = a_helix->getDistanceToPoint(cl_cog,dist);
    if(dist[2]<dist_track)
      dist_track=dist[2];
    for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
      float hitPos[3];
      hitPos[0]=myHit->hit->getPosition()[0];
      hitPos[1]=myHit->hit->getPosition()[1];
      hitPos[2]=myHit->hit->getPosition()[2];
      dummy = a_helix->getDistanceToPoint(hitPos,dist);
      if(dist[2]<shortest_dist_track) {
	shortest_dist_track=dist[2];
	nearestTrack = a_track;
      }
    }
  }
  clusPar->distToTrack=dist_track;
  clusPar->smallestDistToTrack=shortest_dist_track;

  //check if cluster envelopes a large number of removed hits completely (only if cluster is near enough)
  
  if(shortest_dist_track < 3.5*_geomParams->Get_padSizeEcal() [1]) {
    
    map<ExtendedTrack*,map<int,vector<ExtendedHit*> > >::iterator track_it;
    vector<ExtendedHit*>::iterator hit_it;
    
    for(track_it=allRemovedHits->begin();track_it!=allRemovedHits->end();track_it++) {
      if( ((track_it->first)) == (nearestTrack) ) {
	if (_algoParams->GetDebug()>1) cout << "Found nearest Track" << endl;
	bool surrounded[29];
	for(int lay=0;lay<29;lay++) {
	  if(_algoParams->GetDebug()>2)
	    cout << "Layer " << lay << endl;
	  vector<ExtendedHit*> *rHits = &((track_it->second)[lay]);
	  int NRemovedHitsInLayer = rHits->size();
	  if(_algoParams->GetDebug()>2)
	    cout << "Removed Hits: " << NRemovedHitsInLayer << endl;
	  if(NRemovedHitsInLayer==0) {
	    surrounded[lay]=false;
	    continue;
	  }
	  surrounded[lay]=true;
	  for(hit_it=rHits->begin();hit_it!=rHits->end();hit_it++) {
	    int search_layer=(*(_geomParams->Get_defaultDecoder()))((*(hit_it))->hit)["K-1"];
	    int search_i=(*(_geomParams->Get_defaultDecoder()))((*(hit_it))->hit)["I"];
	    int search_j=(*(_geomParams->Get_defaultDecoder()))((*(hit_it))->hit)["J"];
	    int search_m=(*(_geomParams->Get_defaultDecoder()))((*(hit_it))->hit)["M"];
	    //int search_s=(*(_geomParams->Get_defaultDecoder()))((*(hit_it))->hit)["S-1"];
	    if(_algoParams->GetDebug()>2)
	      cout << "Checking removed hit " << search_i << ", " << search_j << endl; 
	    // check if its neighbours are either cluster hits or removed hits
	    // break the loop if: - hit is on the edge of module 
	    if(search_m>0 && search_m<6) {// barrel
	      if(search_i==0) { // IGNORE HIGH EDGE OF MODULE SINCE PAD IS NOT NECESSARILY 179
		surrounded[lay]=false;
		break;
	      }
	    }	    
	    else if(search_m==0 || search_m==6) { //endcap
	      if(search_j==0 || search_i==0) { // IGNORE HIGH EDGE OF MODULE SINCE PAD IS NOT NECESSARILY 179
		surrounded[lay]=false;
		break;
	      }
	    }
	    bool hitArr[3][3]; // to save state of neighbour hits
	    for(int i=0;i<3;i++) {
	      for(int j=0;j<3;j++) {
		hitArr[i][j]=false;
	      }
	    }
	    hitArr[1][1]=true; // the hit itself is true
	    for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) { // now looping once over cluster is enough
	      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
	      if(myHit == *hit_it)
		continue;
	      int hit_layer=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["K-1"];
	      if(hit_layer!=search_layer)
		continue;
	      int hit_i=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["I"];
	      int hit_j=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["J"];
	      int hit_m=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["M"];
	      //int hit_s=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["S-1"];
	      if(search_m>0 && search_m<6 && search_m!=hit_m) {// barrel
//		if(search_j==0 && hit_j==179) {
		if(search_j==0 && hit_j==(10*_geomParams->Get_nCellsPerWafer()-1)) {
		  hit_j=-1;
		}
//		if(search_j==179 && hit_j==0) {
		if(search_j==(10*_geomParams->Get_nCellsPerWafer()-1) && hit_j==0) {
		  hit_j=(10*_geomParams->Get_nCellsPerWafer());
		}
	      }

	      if( fabs((float)(search_i-hit_i))<2 && fabs((float)(search_j-hit_j))<2 ) {
		hitArr[(search_i-hit_i+1)][(search_j-hit_j+1)] = true;
		if(_algoParams->GetDebug()>2)
		  cout << "Found cluster hit " << hit_i << ", " << hit_j << endl; 
	      }
	    }
	    for(int r_hit_i=0;r_hit_i<NRemovedHitsInLayer;r_hit_i++) { // looping over removed hits
	      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*rHits)[r_hit_i]);
	      if(myHit == *hit_it)
		continue;
	      int hit_layer=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["K-1"];
	      if(hit_layer!=search_layer)
		continue;
	      int hit_i=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["I"];
	      int hit_j=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["J"];
	      int hit_m=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["M"];
	      //int hit_s=(*(_geomParams->Get_defaultDecoder()))(myHit->hit)["S-1"];
	      if(search_m>0 && search_m<6 && search_m!=hit_m) {// barrel
//		if(search_j==0 && hit_j==179) {
		if(search_j==0 && hit_j==(10*_geomParams->Get_nCellsPerWafer()-1)) {
		  hit_j=-1;
		}
//		if(search_j==179 && hit_j==0) {
		if(search_j==(10*_geomParams->Get_nCellsPerWafer()-1) && hit_j==0) {
		  hit_j=(10*_geomParams->Get_nCellsPerWafer());
		}
	      }
	      if( fabs((float)(search_i-hit_i))<2 && fabs((float)(search_j-hit_j))<2 ) {
		hitArr[(search_i-hit_i+1)][(search_j-hit_j+1)] = true;
		if(_algoParams->GetDebug()>2)
		  cout << "Found removed hit " << hit_i << ", " << hit_j << endl; 
	      }
	    }
	    for(int i=0;i<3;i++) {
	      for(int j=0;j<3;j++) {
		if(hitArr[i][j]==false) {
		  surrounded[lay]=false;
		  break;
		}
	      }
	      if(surrounded[lay]==false) 
		break;
	    }
	    if(surrounded[lay]==false) {
	      break;
	    }
	  }

	}
	unsigned int surrounded_layers=0;
	for(int lay=0;lay<29;lay++) {
	  if(surrounded[lay]==true) {
	    surrounded_layers++;
	    if(_algoParams->GetDebug()>2)
	      cout << "Layer " << lay << " surrounds Track" << endl;
	  }
	}
	clusPar->surroundingLayers=surrounded_layers;
	break;
      }
    }
  }
  //else
    //cout << "Cluster too far off track" << endl;


  int nPreshowerHits=0;
  int nPL0Hits=0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->preShower)
      nPreshowerHits++;
    if(myHit->pseudoLayer==0)
      nPL0Hits++;
  }
  clusPar->nHits=(NClusteredHits+NCountedGhostHits);//-nPreshowerHits);
  clusPar->psHits=nPreshowerHits;
  clusPar->pl0Hits=nPL0Hits;
    
  int NWeightableHits = NClusteredHits;//-nPreshowerHits; 
  //postion correction via phi and theta

  double clus_phi = clusPar->phi;
  int a_gamma;
  if(clus_phi<0) 
    clus_phi=twopi+clus_phi;
  a_gamma=(int)(((clus_phi-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
  double GAMMA=a_gamma*(twopi/_geomParams->Get_symmetry());
  clus_phi = fabs(clus_phi-GAMMA);
      
  if(clusPar->zone==1) {  // barrel
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    clusPar->Etot_g_theta = clusPar->Etot_g*sin(clusPar->theta);
    clusPar->Etot_g_pos = clusPar->Etot_g*sin(clusPar->theta)*cos(clus_phi);
  }

  if(clusPar->zone==2) {  // endcap
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    clusPar->Etot_g_theta = clusPar->Etot_g*fabs(cos(clusPar->theta));
    clusPar->Etot_g_pos = clusPar->Etot_g*fabs(cos(clusPar->theta))*cos(clus_phi);
  }

  if(clusPar->zone==3) {  // overlap: decide from position of cog 
    clusPar->Etot_g_phi = clusPar->Etot_g*cos(clus_phi);
    if(cluster_cog.z<_geomParams->Get_zOfBarrel() && cluster_cog.z>-_geomParams->Get_zOfBarrel()) { //treat as in barrel
      clusPar->Etot_g_theta = clusPar->Etot_g*sin(clusPar->theta);
      clusPar->Etot_g_pos = clusPar->Etot_g*sin(clusPar->theta)*cos(clus_phi);
    }
    else { // treat as in endcap
      clusPar->Etot_g_theta = clusPar->Etot_g*fabs(cos(clusPar->theta));
      clusPar->Etot_g_pos = clusPar->Etot_g*fabs(cos(clusPar->theta))*cos(clus_phi);
    }
  }

  //if(clusPar->zone==1 || clusPar->zone==3) {
  //  //    clusPar->E_GeV_pre = toGeVFctn->Eval(clusPar->Etot_g);
  //  clusPar->E_GeV_pre = ECALGarlicEnergyEstimator::getGeV(clusPar->Etot_g);
  //} else {
  //  //    clusPar->E_GeV_pre = toGeVFctn_EC->Eval(clusPar->Etot_g);
  //  clusPar->E_GeV_pre = ECALGarlicEnergyEstimator::getGeV_EC(clusPar->Etot_g);
  //}

  bool isbarrel = clusPar->zone==1 || clusPar->zone==3;
  clusPar->E_GeV_pre = _energyEstimator->getGeV(clusPar->Etot_g, isbarrel);


  clusPar->E_GeV_noLC = clusPar->E_GeV;
  clusPar->Etot_g_noLC = clusPar->Etot_g;

  clusPar->E_GeV = clusPar->E_GeV_pre;
  
  if(_algoParams->GetCorrectPhi()) {
    if(clusPar->zone==1 || clusPar->zone==3) {
      //      double corrected_energy = CorrectPhi(clusPar);

      double corrected_energy = _energyEstimator->correctEnergyPhi(clusPar->E_GeV, clusPar->phi);


      clusPar->E_GeV_pre_noPhi = clusPar->E_GeV;
      clusPar->E_GeV_pre = corrected_energy ;      
      clusPar->E_GeV = clusPar->E_GeV_pre;
    }
  }


  //double estimated_energy_en = EstimateEnergyByEnergy(clusPar);
  //clusPar->E_GeV_en = estimated_energy_en;
  //double estimated_energy_hits = EstimateEnergyByHits(clusPar,clusPar->E_GeV_en);
  //clusPar->E_GeV_hits = estimated_energy_hits;
  if(_geomParams->Get_nCellsPerWafer()==18) {
    double estimated_energy_mix = _energyEstimator->EstimateEnergyByMix(clusPar,clusPar->E_GeV_pre);
    clusPar->E_GeV_mix = estimated_energy_mix;
  //if(_optimiseResolution) {
    //double estimated_energy = EstimateEnergy(clusPar,clusPar->E_GeV_en);
    //double estimated_energy = EstimateEnergy(clusPar,clusPar->E_GeV_pre);
  //clusPar->E_GeV_opt = estimated_energy;
  //}
    clusPar->E_GeV = clusPar->E_GeV_mix ;
    clusPar->E_GeV_noTheta = clusPar->E_GeV;
    clusPar->E_GeV_noPhi = clusPar->E_GeV;
    if(_algoParams->GetCorrectTheta()) {
      if(clusPar->zone==1 || clusPar->zone==3) {
	//	double corrected_energy = CorrectTheta(clusPar);

	double corrected_energy = _energyEstimator->correctEnergyTheta(clusPar->E_GeV, clusPar->theta);


	clusPar->E_GeV_noTheta = clusPar->E_GeV;
	clusPar->E_GeV = corrected_energy ;
      }
    }
    if(_algoParams->GetCorrectPhi()) {
      if(clusPar->zone==1 || clusPar->zone==3) {
	//	double corrected_energy = CorrectPhi(clusPar);

	double corrected_energy = _energyEstimator->correctEnergyPhi(clusPar->E_GeV, clusPar->phi);

	clusPar->E_GeV_noPhi = clusPar->E_GeV;
	clusPar->E_GeV = corrected_energy ;      
      }
    }
  }

  // Fill Float arrays of cluster hits
  float *aEn = new float[NWeightableHits];
  float *xPos = new float[NWeightableHits];
  float *yPos = new float[NWeightableHits];
  float *zPos = new float[NWeightableHits];
  int seed_hits=0;
  int hit_cnt=0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    //if(myHit->preShower==1)
    // continue;
    aEn[hit_cnt] = (myHit->hit)->getEnergy();
    xPos[hit_cnt] = (myHit->hit)->getPosition()[0];
    yPos[hit_cnt] = (myHit->hit)->getPosition()[1];
    zPos[hit_cnt] = (myHit->hit)->getPosition()[2];
    hit_cnt++;
    if(myHit->pseudoLayer>_algoParams->GetNLayersForSeeding()-1)
      continue;
    else
      seed_hits++;
  }
  ClusterShapes clusShape(NWeightableHits,aEn,xPos,yPos,zPos);
  float *cog_byShape = clusShape.getCentreOfGravity();
  if(_algoParams->GetDebug()>2)
    cout << "COG by shape : " << cog_byShape[0] << ", " << cog_byShape[1] << ", " << cog_byShape[2] << endl;  
  float *main_axes = clusShape.getEigenVecInertia();
  if(_algoParams->GetDebug()>2)
    cout << "Main axis : " << main_axes[0] << ", " << main_axes[1] << ", " << main_axes[2] << endl;  
  float width = clusShape.getWidth();
  if(_algoParams->GetDebug()>2)
    cout << "Mean width : " << width << endl;
  float e_vol = clusShape.getElipsoid_vol();
  float en_density = clusShape.getElipsoid_density();
  float hit_density = NWeightableHits/e_vol;
  if(_algoParams->GetDebug()>2)
    cout << "Hit density: " << hit_density << ", Energy density: " << en_density << endl;
  float eccentr = clusShape.getElipsoid_eccentricity();
  if(_algoParams->GetDebug()>2)
    cout << "Eccentricity: " << eccentr << endl;
  
  delete aEn;
  delete xPos;
  delete yPos;
  delete zPos;

  double cosDiff = ((cog_byShape[0]*main_axes[0])+(cog_byShape[1]*main_axes[1])+(cog_byShape[2]*main_axes[2]))/(sqrt((cog_byShape[0]*cog_byShape[0])+(cog_byShape[1]*cog_byShape[1])+(cog_byShape[2]*cog_byShape[2]))*sqrt((main_axes[0]*main_axes[0])+(main_axes[1]*main_axes[1])+(main_axes[2]*main_axes[2])));
  
  if( (!isnan(cosDiff)) && (std::isfinite(cosDiff)) )
    clusPar->dirErr = cosDiff;
  else
    clusPar->dirErr = 0;

  if( (!isnan(hit_density)) && (std::isfinite(hit_density)) )
    clusPar->hitDensity = hit_density;
  else
    clusPar->hitDensity = 0;
  if( (!isnan(en_density)) && (std::isfinite(en_density)) )
    clusPar->enDensity = en_density;
  else
    clusPar->enDensity = 0;
  if( (!isnan(eccentr)) && (std::isfinite(eccentr)) )
    clusPar->Eccentricity = eccentr;
  else
    clusPar->Eccentricity = 1;
  if( (!isnan(width)) && (std::isfinite(width)) )
    clusPar->Width = width;
  else
    clusPar->Width = 0;
  if( (!isnan(e_vol)) && (std::isfinite(e_vol)) )   
    clusPar->Volume = e_vol;
  else
    clusPar->Volume = 0;

   //same for seed
  float *seed_aEn = new float[seed_hits];
  float *seed_xPos = new float[seed_hits];
  float *seed_yPos = new float[seed_hits];
  float *seed_zPos = new float[seed_hits];
  hit_cnt=0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->pseudoLayer>_algoParams->GetNLayersForSeeding()-1)
      continue;
    //if(myHit->preShower==1)
    //  continue;
    seed_aEn[hit_cnt]=(myHit->hit)->getEnergy();
    seed_xPos[hit_cnt]=(myHit->hit)->getPosition()[0];
    seed_yPos[hit_cnt]=(myHit->hit)->getPosition()[1];
    seed_zPos[hit_cnt]=(myHit->hit)->getPosition()[2];
    hit_cnt++;
  }
  ClusterShapes seedShape(seed_hits,seed_aEn,seed_xPos,seed_yPos,seed_zPos);
  float *seed_cog_byShape = seedShape.getCentreOfGravity();
  if(_algoParams->GetDebug()>2)
    cout << "Seed COG by shape : " << seed_cog_byShape[0] << ", " << seed_cog_byShape[1] << ", " << seed_cog_byShape[2] << endl;  
  float *seed_main_axes = seedShape.getEigenVecInertia();
  if(_algoParams->GetDebug()>2)
    cout << "Seed Main axis : " << seed_main_axes[0] << ", " << seed_main_axes[1] << ", " << seed_main_axes[2] << endl;  
  float seed_width = seedShape.getWidth();
  if(_algoParams->GetDebug()>2)
    cout << "Seed Mean width : " << seed_width << endl;
  float seed_e_vol = seedShape.getElipsoid_vol();
  float seed_en_density = seedShape.getElipsoid_density();
  float seed_hit_density = seed_hits/seed_e_vol;
  if(_algoParams->GetDebug()>2)
    cout << "Seed Hit density: " << seed_hit_density << ", Energy density: " << seed_en_density << endl;
  float seed_eccentr = seedShape.getElipsoid_eccentricity();
  if(_algoParams->GetDebug()>2)
    cout << "Seed Eccentricity: " << seed_eccentr << endl;
   
  delete seed_aEn;
  delete seed_xPos;
  delete seed_yPos;
  delete seed_zPos;

  double seed_cosDiff = ((seed_cog_byShape[0]*seed_main_axes[0])+(seed_cog_byShape[1]*seed_main_axes[1])+(seed_cog_byShape[2]*seed_main_axes[2]))/(sqrt((seed_cog_byShape[0]*seed_cog_byShape[0])+(seed_cog_byShape[1]*seed_cog_byShape[1])+(seed_cog_byShape[2]*seed_cog_byShape[2]))*sqrt((seed_main_axes[0]*seed_main_axes[0])+(seed_main_axes[1]*seed_main_axes[1])+(seed_main_axes[2]*seed_main_axes[2])));


  if( (!isnan(seed_cosDiff)) && (std::isfinite(seed_cosDiff)) )
    clusPar->seedDirErr = seed_cosDiff;
  else
    clusPar->seedDirErr = 0;

  /*
  CalculateSeedDirection(myCluster);
  CalculateClusterDirection(myCluster);
  FitEllipsoid(myCluster);
  
  clusPar->Density = ((myCluster->Shape)->Volume)/(clusPar->nHits);

  float r_phi_dir=sqrt((myCluster->dir).x*(myCluster->dir).x+(myCluster->dir).y*(myCluster->dir).y);
  float theta_dir = atan2(r_phi_dir,(myCluster->dir).z);
  if(theta_dir<0) 
    theta_dir=twopi+theta_dir;
  float phi_dir = atan2((myCluster->dir).x,(myCluster->dir).y);
  if(phi_dir<0) 
    phi_dir=twopi+phi_dir;
  float dir_err = sqrt(pow(phi_dir-clus_phi,2)+pow((theta_dir-(clusPar->theta)),2));
  clusPar->dirErr = dir_err;

  r_phi_dir=sqrt((myCluster->seed_dir).x*(myCluster->seed_dir).x+(myCluster->seed_dir).y*(myCluster->seed_dir).y);
  theta_dir = atan2(r_phi_dir,(myCluster->seed_dir).z);
  if(theta_dir<0) 
    theta_dir=twopi+theta_dir;
  phi_dir = atan2((myCluster->seed_dir).x,(myCluster->seed_dir).y);
  if(phi_dir<0) 
    phi_dir=twopi+phi_dir;
  float seed_dir_err = sqrt(pow(phi_dir-clus_phi,2)+pow((theta_dir-(clusPar->theta)),2));
  clusPar->seedDirErr = seed_dir_err;
  */

  // Fit longitudinal shower shape
  float chi2_long = 0;
  clusPar->Chi2_long=chi2_long;
  if(_algoParams->GetDebug()>2) {
    cout << "Cluster Paramerers of cluster " << clus_ID << " : " << endl
	 << "ID: " << clusPar->ID << endl
	 << "Start: " << clusPar->start << endl
	 << "End: " << clusPar->end << endl
	 << "Etot: " << clusPar->Etot << endl
	 << "Etot_g: " << clusPar->Etot_g << endl
	 << "Etot_g_pos: " << clusPar->Etot_g_pos << endl
	 << "E_GeV: " << clusPar->E_GeV << endl
      //	 << "Es1ovEtot: " << clusPar->Es1ovEtot << endl
	 << "E1C: " << clusPar->E1C << endl
	 << "E4C: " << clusPar->E4C << endl
	 << "NHits: " << clusPar->nHits << endl
	 << "Theta: " << clusPar->theta << endl
	 << "Phi: " << clusPar->phi << endl
	 << "<depth>: " << clusPar->depth << endl
	 << "COG at: " << clusPar->COGx << ", " << clusPar->COGy << ", " << clusPar->COGz << endl
	 << "Distance to biggest: " << clusPar->distToBiggest << endl
	 << "Smallest Distance to biggest: " << clusPar->smallestDistToBiggest << endl
      	 << "Distance to track: " << clusPar->distToTrack << endl
	 << "Smallest Distance to track: " << clusPar->smallestDistToTrack << endl
	 << "Chi2_long: " << clusPar->Chi2_long << endl 
      	 << "Hit Desnsity: " << clusPar->hitDensity << endl 
      	 << "Energy Desnsity: " << clusPar->enDensity << endl 
	 << "Direction Error: " << clusPar->dirErr << endl 
	 << "Zone: " << clusPar->zone << endl;
  }

  clusPar->signalProbMult = 0;
  clusPar->bckgrndProbMult = 0;
  clusPar->photonProbMult = 0;

  // daniel removed MC information from here

  clusPar->isReal=1;

  double MLPout = -999;

  vector<double> inputVec;

  if(clusPar->E_GeV<=3.) {
    
    inputVec.push_back(clusPar->start);
    inputVec.push_back(clusPar->end);
    inputVec.push_back(clusPar->depth);
    inputVec.push_back((clusPar->E1C)/(clusPar->Etot_g));
    inputVec.push_back((clusPar->E4C)/(clusPar->Etot_g));
    inputVec.push_back((clusPar->E9C)/(clusPar->Etot_g));
    inputVec.push_back(clusPar->dirErr);
    inputVec.push_back(clusPar->seedDirErr);
    inputVec.push_back(clusPar->Width);
    inputVec.push_back(clusPar->Eccentricity);
    inputVec.push_back(clusPar->nHits);
    inputVec.push_back(clusPar->Volume);
    inputVec.push_back(clusPar->depth-clusPar->start);

    if(clusPar->zone==1 || clusPar->zone==3) {

      if(clusPar->E_GeV<=0.25) {
	MLPout = MLPResponse0_25_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.25 && clusPar->E_GeV<=0.35) {
	MLPout = MLPResponse0_35_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.35 && clusPar->E_GeV<=0.5) {
	MLPout = MLPResponse0_5_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.5 && clusPar->E_GeV<=0.75) {
      MLPout = MLPResponse0_75_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.75 && clusPar->E_GeV<=1.) {
	MLPout = MLPResponse1_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1. && clusPar->E_GeV<=1.25) {
	MLPout = MLPResponse1_25_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.25 && clusPar->E_GeV<=1.5) {
	MLPout = MLPResponse1_5_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.5 && clusPar->E_GeV<=1.75) {
	MLPout = MLPResponse1_75_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.75 && clusPar->E_GeV<=2.25) {
	MLPout = MLPResponse2_25_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>2.25 && clusPar->E_GeV<=3.) {
	MLPout = MLPResponse3_B->GetMvaValue(inputVec);
      }
    }

    if(clusPar->zone==2) {
      if(clusPar->E_GeV<=0.25) {
	MLPout = MLPResponse0_25_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.25 && clusPar->E_GeV<=0.35) {
	MLPout = MLPResponse0_35_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.35 && clusPar->E_GeV<=0.5) {
	MLPout = MLPResponse0_5_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.5 && clusPar->E_GeV<=0.75) {
      MLPout = MLPResponse0_75_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>0.75 && clusPar->E_GeV<=1.) {
	MLPout = MLPResponse1_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1. && clusPar->E_GeV<=1.25) {
	MLPout = MLPResponse1_25_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.25 && clusPar->E_GeV<=1.5) {
	MLPout = MLPResponse1_5_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.5 && clusPar->E_GeV<=1.75) {
	MLPout = MLPResponse1_75_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>1.75 && clusPar->E_GeV<=2.25) {
	MLPout = MLPResponse2_25_EC->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>2.25 && clusPar->E_GeV<=3.) {
	MLPout = MLPResponse3_EC->GetMvaValue(inputVec);
      }
    }
  } else { // higher energy >3GeV

    inputVec.push_back(clusPar->start);
    inputVec.push_back(clusPar->end);
    inputVec.push_back(clusPar->depth);
    inputVec.push_back((clusPar->E1C)/(clusPar->Etot_g));
    inputVec.push_back((clusPar->E4C)/(clusPar->Etot_g));
    inputVec.push_back((clusPar->E9C)/(clusPar->Etot_g));
    inputVec.push_back(clusPar->dirErr);
    inputVec.push_back(clusPar->seedDirErr);
    inputVec.push_back((clusPar->Es3)/(clusPar->Etot_g));
    inputVec.push_back(clusPar->Width);
    inputVec.push_back(clusPar->Eccentricity);
    inputVec.push_back(clusPar->nHits);
    inputVec.push_back(clusPar->Volume);
    inputVec.push_back(clusPar->depth-clusPar->start);

    if(clusPar->zone==1 || clusPar->zone==3) {
      if(clusPar->E_GeV>3. && clusPar->E_GeV<=5.) {
	if(_geomParams->Get_nCellsPerWafer()==18)
	  MLPout = MLPResponse5_B->GetMvaValue(inputVec);
	else
	  MLPout = MLPResponse10_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>5. && clusPar->E_GeV<=10.) {
	MLPout = MLPResponse10_B->GetMvaValue(inputVec);
      }
      if(clusPar->E_GeV>10. && clusPar->E_GeV<=20.) {
	if(_geomParams->Get_nCellsPerWafer()==18)
	  MLPout = MLPResponse20_B->GetMvaValue(inputVec);
	else
	  MLPout = 1;
      }
      if(clusPar->E_GeV>20.) {
	MLPout = 1;
      }
    } else if (clusPar->zone==2) {
      if(clusPar->E_GeV>3. && clusPar->E_GeV<=5.) {
	if(_geomParams->Get_nCellsPerWafer()==18)
	  MLPout = MLPResponse5_EC->GetMvaValue(inputVec);
	else
	  MLPout = MLPResponse10_EC->GetMvaValue(inputVec);
      } else if (clusPar->E_GeV>5. && clusPar->E_GeV<=10.) {
	MLPout = MLPResponse10_EC->GetMvaValue(inputVec);
      } else if (clusPar->E_GeV>10. && clusPar->E_GeV<=20.) {
	if(_geomParams->Get_nCellsPerWafer()==18) {
	  MLPout = MLPResponse20_EC->GetMvaValue(inputVec);
	} else {
	  MLPout = 1;
	}
      }
      if(clusPar->E_GeV>20.) {
	MLPout = 1;
      }
    }
    //if(clusPar->E_GeV>50) {
    //MLPout = 1;
    //}    
  }

  clusPar->MLP = MLPout;
  if (_algoParams->GetDebug()>1) cout << "TMVA classifier outputs for cluster " << clusPar->ID << ": MLP: " <<  clusPar->MLP << endl;

  /*
  vector<double> input_three;

  //input_three.push_back(clusPar->smallestDistToTrack);
  //input_three.push_back(clusPar->MLP);
  //input_three.push_back(clusPar->E_GeV);

  //float finalCut = ThreeDimResponse->GetMvaValue(input_three);
  clusPar->finalCut=finalCut;
  */

}


void ECALGarlicClusterHelpers::CalculateSeedDirection(ExtendedCluster* a_cluster)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  //calculate COG of seed
  vec3 COG;
  vector<ExtendedHit* > *clusterHits = &(a_cluster->hitVec);
  int NClusteredHits = clusterHits->size();
  double SumE = 0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->pseudoLayer>_algoParams->GetNLayersForSeeding()-1)
      continue;
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    SumE+=hit_en;
    COG.x+=(a_hit->getPosition()[0])*hit_en;
    COG.y+=(a_hit->getPosition()[1])*hit_en;
    COG.z+=(a_hit->getPosition()[2])*hit_en;
  }
  COG.x/=SumE;
  COG.y/=SumE;
  COG.z/=SumE;

  cout << "Seed COG: " << COG.x << ", " << COG.y << ", " << COG.z << endl;

  TVector3 uu(COG.x,COG.y,COG.z);
  TVector3 TCOG(COG.x,COG.y,COG.z);
  TMatrixD Sums(3,3);
  SumE = 0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->pseudoLayer>_algoParams->GetNLayersForSeeding()-1)
      continue;
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    SumE+=hit_en;
    vec3 hitPos;
    hitPos.x=(a_hit->getPosition()[0]);
    hitPos.y=(a_hit->getPosition()[1]);
    hitPos.z=(a_hit->getPosition()[2]);
    TVector3 xc(hitPos.x,hitPos.y,hitPos.z); //global coordinates
    xc -= TCOG; //local coordinates
    TVector3 X(xc.Dot(uu[0]), xc.Dot(uu[1]), xc.Dot(uu[2]));
    for(int d1=0; d1<3; d1++)
      for(int d2=0; d2<d1+1; d2++)
	Sums(d1, d2) += hit_en * X[d1] * X[d2];
  }
  
  // Angle of the main components
  double tanPhi = Sums(1,0)/Sums(0,0);
  double tanTh = Sums(2,0)/Sums(0,0);
  
  double cosPhi = 1.0/sqrt(1.0+tanPhi*tanPhi);
  double sinPhi = tanPhi*cosPhi;
  double cosTh = 1.0/sqrt(1.0+tanTh*tanTh);
  double sinTh = tanTh*cosTh;
  
  //  if(_algoParams->GetDebug()>2)
    cout << "Seed direction: cosPhi, cosTh = " << cosPhi << ", " << cosTh << endl;

  // New main direction vector
  TVector3 nu = uu[0]*cosPhi*cosTh + uu[1]*sinPhi*cosTh + uu[2]*cosPhi*sinTh;
  vec3 main_axis;
  main_axis.x = nu.X();
  main_axis.y = nu.Y();
  main_axis.z = nu.Z();

  cout << "Seed MainAxis: " << main_axis.x << ", " << main_axis.y << ", " << main_axis.z << endl;

  a_cluster->seed_dir=main_axis;
}


void ECALGarlicClusterHelpers::CalculateClusterDirection(ExtendedCluster* a_cluster)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  //initial direction is the COG vector
  vec3 COG;
  COG.x = (a_cluster->Parameters)->COGx;
  COG.y = (a_cluster->Parameters)->COGy;
  COG.z = (a_cluster->Parameters)->COGz;
  TVector3 uu(COG.x,COG.y,COG.z);
  TVector3 TCOG(COG.x,COG.y,COG.z);
  TMatrixD Sums(3,3);
  double SumE = 0;

  cout << "Cluster COG: " << COG.x << ", " << COG.y << ", " << COG.z << endl;

  vector<ExtendedHit* > *clusterHits = &(a_cluster->hitVec);
  int NClusteredHits = clusterHits->size();
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    SumE+=hit_en;
    vec3 hitPos;
    hitPos.x=(a_hit->getPosition()[0]);
    hitPos.y=(a_hit->getPosition()[1]);
    hitPos.z=(a_hit->getPosition()[2]);
    TVector3 xc(hitPos.x,hitPos.y,hitPos.z); //global coordinates
    xc -= TCOG; //local coordinates
    TVector3 X(xc.Dot(uu[0]), xc.Dot(uu[1]), xc.Dot(uu[2]));
    for(int d1=0; d1<3; d1++)
      for(int d2=0; d2<d1+1; d2++)
	Sums(d1, d2) += hit_en * X[d1] * X[d2];
  }
  
  // Angle of the main components
  double tanPhi = Sums(1,0)/Sums(0,0);
  double tanTh = Sums(2,0)/Sums(0,0);
  
  double cosPhi = 1.0/sqrt(1.0+tanPhi*tanPhi);
  double sinPhi = tanPhi*cosPhi;
  double cosTh = 1.0/sqrt(1.0+tanTh*tanTh);
  double sinTh = tanTh*cosTh;
  
  if(_algoParams->GetDebug()>2)
    cout << "Cluster direction: cosPhi, cosTh = " << cosPhi << ", " << cosTh << endl;

  // New main direction vector
  TVector3 nu = uu[0]*cosPhi*cosTh + uu[1]*sinPhi*cosTh + uu[2]*cosPhi*sinTh;
  vec3 main_axis;
  main_axis.x = nu.X();
  main_axis.y = nu.Y();
  main_axis.z = nu.Z();

  cout << "Cluster MainAxis: " << main_axis.x << ", " << main_axis.y << ", " << main_axis.z << endl;

  a_cluster->dir=main_axis;
}


void ECALGarlicClusterHelpers::FitEllipsoid(ExtendedCluster* a_cluster)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  Ellipsoid *shape = new Ellipsoid();
  vec3 main_axis = a_cluster->dir;
  TVector3 uu(main_axis.x,main_axis.y,main_axis.z);
  vec3 COG;
  COG.x = (a_cluster->Parameters)->COGx;
  COG.y = (a_cluster->Parameters)->COGy;
  COG.z = (a_cluster->Parameters)->COGz;
  TVector3 TCOG(COG.x,COG.y,COG.z);
  TMatrixD Sums(3,3);
  double SumE = 0;
  vector<ExtendedHit* > *clusterHits = &(a_cluster->hitVec);
  int NClusteredHits = clusterHits->size();
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    SumE+=hit_en;
    vec3 hitPos;
    hitPos.x=(a_hit->getPosition()[0]);
    hitPos.y=(a_hit->getPosition()[1]);
    hitPos.z=(a_hit->getPosition()[2]);
    TVector3 xc(hitPos.x,hitPos.y,hitPos.z); //global coordinates
    xc -= TCOG; //local coordinates
    TVector3 X(xc.Dot(uu[0]), xc.Dot(uu[1]), xc.Dot(uu[2]));
    for(int d1=0; d1<3; d1++)
      for(int d2=0; d2<d1+1; d2++)
	Sums(d1, d2) += hit_en * X[d1] * X[d2];
  }

  shape->WidthR = sqrt( (Sums(1,1)+Sums(2,2))/SumE + pow(_geomParams->Get_padSizeEcal() [0],2)/12 );
  shape->WidthL = sqrt( Sums(0,0)/SumE + pow(_geomParams->Get_padSizeEcal() [0],2)/12);
  shape->Volume = 4/3*pi*(shape->WidthL)*(shape->WidthR)*(shape->WidthR);
  a_cluster->Shape=shape;
}

vec3 ECALGarlicClusterHelpers::FollowClusterDirection(vec3 *myPos, int layersFromStart, 
						      vec3 *clusterDirection, ClusterLocation location, 
						      bool &gapHoleVeto, int &gapJumpTo)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  if(_algoParams->GetDebug()>2) {
    cout << "Our location is " << location << endl;
    cout << "Cluster direction: " << clusterDirection->x << ", " << clusterDirection->y << ", " << clusterDirection->z << ", " << endl;
  }
  if(layersFromStart==0)
    return *myPos;
  vec3 endpoint;
  float zhit = fabs(myPos->z);
  // here we determine the pseudo layer for myPos: for hits this should give the same as for the real hit, for the ECAL front plane projections we should end up in layer 0 FIXME: Do we???
  double prodRadiusMax=0.;
  
  for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
    double prodRadius = myPos->x*_geomParams->Get_barrelStaveDir()[istave].x+myPos->y*_geomParams->Get_barrelStaveDir()[istave].y;
    if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
  }

  //int   bestLayer=0;
  int   bestBarrelLayer=0;
  double barrelLayerDist = 9999.;
  int   bestEndcapLayer=0;
  double endcapLayerDist = 9999.;
  
  int startPseudoLayer=200;
  double scale=1;
  gapJumpTo=0;

  switch(location){
  case CLUS_LOCATION_BARREL: {
    if(_algoParams->GetDebug()>2)
      cout << "Following cluster direction in barrel" << endl;
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	bestBarrelLayer = ilayer;
      }
    }
    startPseudoLayer=bestBarrelLayer;
    if(_algoParams->GetDebug()>2)
      cout << "...starting from ps layer " << startPseudoLayer << endl;
    double phi_pos = atan2((double)clusterDirection->x,(double)clusterDirection->y);
    if(_algoParams->GetDebug()>2)
      cout << "PHI pos: " << phi_pos << endl;
    if(phi_pos<(twopi/_geomParams->Get_symmetry())&&phi_pos>(-twopi/_geomParams->Get_symmetry())) {
      if(_algoParams->GetDebug()>2)
	cout << "Stave angle: 0" << endl;
      vec3 direction;	    
      direction.x = clusterDirection->x;
      direction.y = clusterDirection->y;
      direction.z = clusterDirection->z;
      double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[startPseudoLayer];
      if(_algoParams->GetDebug()>2)
	cout << "The distance between ps layer " << startPseudoLayer+layersFromStart << " and " << startPseudoLayer << " is " << y_dist << " , so the scale is " << fabs(y_dist/direction.y) << endl;
      scale = fabs(y_dist/direction.y);
      direction.x=direction.x*scale;
      direction.y=direction.y*scale;
      direction.z=direction.z*scale;
      endpoint.x = myPos->x+direction.x;
      endpoint.y = myPos->y+direction.y;
      endpoint.z = myPos->z+direction.z;
      }
    else {
      int gamma;
      if(phi_pos<0) 
	phi_pos=twopi+phi_pos;
      gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
      if(_algoParams->GetDebug()>2)
	cout << "Symmetry " << _geomParams->Get_symmetry() << endl;
      double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
      if(_algoParams->GetDebug()>2)
	cout << "Stave angle: " << GAMMA << " and gamma: "<< gamma << endl;
      float X_DIR=clusterDirection->x;
      float Y_DIR=clusterDirection->y;
      vec3 direction;	    
      direction.x = X_DIR*cos(GAMMA)-Y_DIR*sin(GAMMA);
      direction.y = Y_DIR*cos(GAMMA)+X_DIR*sin(GAMMA);
      direction.z = clusterDirection->z;
      float X_POS=myPos->x;
      float Y_POS=myPos->y;
      vec3 start_pos;
      start_pos.x = X_POS*cos(GAMMA)-Y_POS*sin(GAMMA);
      start_pos.y = Y_POS*cos(GAMMA)+X_POS*sin(GAMMA);
      start_pos.z = myPos->z;
      // now in local frame: find position n layers deeper
      // find y distance
      double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[startPseudoLayer];
      if(_algoParams->GetDebug()>2)
	cout << "The distance between ps layer " << startPseudoLayer+layersFromStart << " and " << startPseudoLayer << " is " << y_dist << " , so the scale is " << fabs(y_dist/direction.y) << endl;
      scale = fabs(y_dist/direction.y);
      direction.x=direction.x*scale;
      direction.y=direction.y*scale;
      direction.z=direction.z*scale;
      X_POS=start_pos.x+direction.x;
      Y_POS=start_pos.y+direction.y;
      if(_algoParams->GetDebug()>2)
      cout << "In local frame new point is: " << X_POS << ", " << Y_POS << endl;
      //no go back to global frame
      endpoint.x = X_POS*cos(GAMMA)+Y_POS*sin(GAMMA);
      endpoint.y = Y_POS*cos(GAMMA)-X_POS*sin(GAMMA);
      endpoint.z = myPos->z+direction.z;
    }
    break;
  }
  case CLUS_LOCATION_ENDCAP:  {
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(fabs(zhit)-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    startPseudoLayer=bestEndcapLayer;
    double z_dist=_geomParams->Get_positionEndcapLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionEndcapLayer()[startPseudoLayer];
    scale = fabs(z_dist/clusterDirection->z);
    endpoint.x = myPos->x+(clusterDirection->x*scale);
    endpoint.y = myPos->y+(clusterDirection->y*scale);
    endpoint.z = myPos->z+(clusterDirection->z*scale);
    break;
  }
  case CLUS_LOCATION_OVERLAP:  {
    if(_algoParams->GetDebug()>2)
      cout << "Following cluster direction in Overlap Zone: " << endl;
    
    if(myPos->z < _geomParams->Get_zOfBarrel() && myPos->z > -_geomParams->Get_zOfBarrel()) { //Barrel hit, treat carefully
      if(_algoParams->GetDebug()>2) cout << "From a barrel hit" << endl;
      double dist_from_edge = fabs(_geomParams->Get_zOfBarrel()-fabs(myPos->z));
      if(_algoParams->GetDebug()>2) cout << "Distance from the edge is " << dist_from_edge << endl;
      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	  barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	  bestBarrelLayer = ilayer;
	}
      }
      startPseudoLayer=bestBarrelLayer;
      if(_algoParams->GetDebug()>2) cout << "...starting from ps layer " << startPseudoLayer << endl;
      double phi_pos = atan2((double)clusterDirection->x,(double)clusterDirection->y);
      if(_algoParams->GetDebug()>2) cout << "PHI pos: " << phi_pos << endl;
      if(phi_pos<(twopi/_geomParams->Get_symmetry())&&phi_pos>(-twopi/_geomParams->Get_symmetry())) {
	if(_algoParams->GetDebug()>2) cout << "Stave angle: 0" << endl;
	vec3 direction;	    
	direction.x = clusterDirection->x;
	direction.y = clusterDirection->y;
	direction.z = clusterDirection->z;
	double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[startPseudoLayer];
	if(_algoParams->GetDebug()>2) cout << "The distance between ps layer " << startPseudoLayer+layersFromStart << " and " << startPseudoLayer << " is " << y_dist << " , so the scale is " << fabs(y_dist/direction.y) << endl;
	scale = fabs(y_dist/direction.y);
	if((fabs(direction.z)*scale)<dist_from_edge) {
	  direction.x=direction.x*scale;
	  direction.y=direction.y*scale;
	  direction.z=direction.z*scale;
	  endpoint.x = myPos->x+direction.x;
	  endpoint.y = myPos->y+direction.y;
	  endpoint.z = myPos->z+direction.z;
	}
	else { // cant go all the way over the edge...
	  if(_algoParams->GetDebug()>2) cout << "Stepping over the edge..." << endl;
	  double z_scale = fabs(dist_from_edge/direction.z);
	  // now bridge cable gap
	  vec3 exitBarrel;
	  exitBarrel.x = myPos->x+(direction.x*z_scale);
	  exitBarrel.y = myPos->y+(direction.y*z_scale);
	  exitBarrel.z = myPos->z+(direction.z*z_scale);
	  double gap_scale=fabs((_geomParams->Get_positionEndcapLayer()[0]-fabs(exitBarrel.z))/direction.z);
	  vec3 enterEndcap;
	  enterEndcap.x=exitBarrel.x+(direction.x*gap_scale)*_algoParams->GetXY_gap_transit_factor();
	  enterEndcap.y=exitBarrel.y+(direction.y*gap_scale)*_algoParams->GetXY_gap_transit_factor();
	  if(exitBarrel.z>0)
	    enterEndcap.z=_geomParams->Get_positionEndcapLayer()[0];
	  else
	    enterEndcap.z=-_geomParams->Get_positionEndcapLayer()[0];
	  int enterEndcapInPsLayer=0;
	  double psDist=9999;
	  for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	    if(fabs(enterEndcap.y-_geomParams->Get_positionBarrelLayer()[ilayer])<psDist){
	      psDist=fabs(enterEndcap.y-_geomParams->Get_positionBarrelLayer()[ilayer]);
	      enterEndcapInPsLayer = ilayer;
	    }
	  }
	  if(_algoParams->GetDebug()>2) cout << "...and would enter Endcap in PSlayer " << enterEndcapInPsLayer << " at " << 
			 enterEndcap.x << ", " << enterEndcap.y << ", " << enterEndcap.z << endl;
	  if(enterEndcapInPsLayer>(startPseudoLayer+layersFromStart)) {
	    gapHoleVeto=true;
	    endpoint.x=0;	    
	    endpoint.y=0;	    
	    endpoint.z=0;
	    break;
	  }	    
	  else {
	    if(enterEndcapInPsLayer==(startPseudoLayer+layersFromStart))
	       gapJumpTo=enterEndcapInPsLayer;
	    double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[enterEndcapInPsLayer];
	    scale = fabs(y_dist/clusterDirection->y);
	    endpoint.x = enterEndcap.x+(clusterDirection->x*scale);
	    endpoint.y = enterEndcap.y+(clusterDirection->y*scale);
	    endpoint.z = enterEndcap.z+(clusterDirection->z*scale);
	  }
	}
      }
      else {
	int gamma;
	if(phi_pos<0) 
	  phi_pos=twopi+phi_pos;
	gamma=(int)(((phi_pos-(twopi/(2*_geomParams->Get_symmetry())))/(twopi/_geomParams->Get_symmetry()))+1);
	if(_algoParams->GetDebug()>2)
	cout << "Symmetry " << _geomParams->Get_symmetry() << endl;
	double GAMMA=gamma*(twopi/_geomParams->Get_symmetry());
	if(_algoParams->GetDebug()>2)
	  cout << "Stave angle: " << GAMMA << " and gamma: "<< gamma << endl;
	float X_DIR=clusterDirection->x;
	float Y_DIR=clusterDirection->y;
	vec3 direction;	    
	direction.x = X_DIR*cos(GAMMA)-Y_DIR*sin(GAMMA);
	direction.y = Y_DIR*cos(GAMMA)+X_DIR*sin(GAMMA);
	direction.z = clusterDirection->z;
	float X_POS=myPos->x;
	float Y_POS=myPos->y;
	vec3 start_pos;
	start_pos.x = X_POS*cos(GAMMA)-Y_POS*sin(GAMMA);
	start_pos.y = Y_POS*cos(GAMMA)+X_POS*sin(GAMMA);
	start_pos.z = myPos->z;
	// now in local frame: find position n layers deeper
	// find y distance
	double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[startPseudoLayer];
	if(_algoParams->GetDebug()>2)
	  cout << "The distance between ps layer " << startPseudoLayer+layersFromStart << " and " << startPseudoLayer << " is " << y_dist << " , so the scale is " << fabs(y_dist/direction.y) << endl;
	scale = fabs(y_dist/direction.y);
	if((fabs(direction.z)*scale)<dist_from_edge) {
	  direction.x=direction.x*scale;
	  direction.y=direction.y*scale;
	  direction.z=direction.z*scale;
	  X_POS=start_pos.x+direction.x;
	  Y_POS=start_pos.y+direction.y;
	  if(_algoParams->GetDebug()>2)
	    cout << "In local frame new point is: " << X_POS << ", " << Y_POS << endl;
	  //no go back to global frame
	  endpoint.x = X_POS*cos(GAMMA)+Y_POS*sin(GAMMA);
	  endpoint.y = Y_POS*cos(GAMMA)-X_POS*sin(GAMMA);
	  endpoint.z = myPos->z+direction.z;
	}
	else { // cant go all the way over the edge...
	  double z_scale = fabs(dist_from_edge/direction.z);
	  // now bridge cable gap
	  vec3 exitBarrel;
	  exitBarrel.x = start_pos.x+(direction.x*z_scale);
	  exitBarrel.y = start_pos.y+(direction.y*z_scale);
	  exitBarrel.z = start_pos.z+(direction.z*z_scale);
	  double gap_scale=fabs((_geomParams->Get_positionEndcapLayer()[1]-fabs(exitBarrel.z))/direction.z);
	  vec3 enterEndcap;
	  enterEndcap.x=exitBarrel.x+(direction.x*gap_scale)*_algoParams->GetXY_gap_transit_factor();
	  enterEndcap.y=exitBarrel.y+(direction.y*gap_scale)*_algoParams->GetXY_gap_transit_factor();
	  if(exitBarrel.z>0)
	    enterEndcap.z=_geomParams->Get_positionEndcapLayer()[0];
	  else
	    enterEndcap.z=-_geomParams->Get_positionEndcapLayer()[0];
	  int enterEndcapInPsLayer=0;
	  double psDist=9999;
	  for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	    if(fabs(enterEndcap.y-_geomParams->Get_positionBarrelLayer()[ilayer])<psDist){
	      psDist=fabs(enterEndcap.y-_geomParams->Get_positionBarrelLayer()[ilayer]);
	      enterEndcapInPsLayer = ilayer;
	    }
	  }
	  if(_algoParams->GetDebug()>2) 
	    cout << "...and would enter Endcap in PSlayer " << enterEndcapInPsLayer << " at " << enterEndcap.x << ", " << enterEndcap.y << ", " << enterEndcap.z << endl;
	  if(enterEndcapInPsLayer>(startPseudoLayer+layersFromStart)) {
	    gapHoleVeto=true;
	    endpoint.x=0;	    
	    endpoint.y=0;	    
	    endpoint.z=0;
	    break;
	  }	    
	  else {
	    if(enterEndcapInPsLayer==(startPseudoLayer+layersFromStart))
	       gapJumpTo=enterEndcapInPsLayer;
	    double y_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[enterEndcapInPsLayer];
	    scale = fabs(y_dist/clusterDirection->y);
	    vec3 localEndpoint;
	    localEndpoint.x = enterEndcap.x+(direction.x*scale);
	    localEndpoint.y = enterEndcap.y+(direction.y*scale);
	    localEndpoint.z = enterEndcap.z+(direction.z*scale);
	    if(_algoParams->GetDebug()>2)
	      cout << "In local frame new point is: " << X_POS << ", " << Y_POS << endl;
	    //no go back to global frame
	    endpoint.x = localEndpoint.x*cos(GAMMA)+localEndpoint.y*sin(GAMMA);
	    endpoint.y = localEndpoint.y*cos(GAMMA)-localEndpoint.x*sin(GAMMA);
	    endpoint.z = localEndpoint.z;
	  }
	}	
      }
    }
    else { // just follow normal endcap procedure
      if(_algoParams->GetDebug()>2)
	cout << "From an endcap hit" << endl;
      for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
	if(fabs(fabs(zhit)-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	  endcapLayerDist = fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer]);
	  bestEndcapLayer = ilayer;
	}
      }
      startPseudoLayer=bestEndcapLayer;
      double z_dist=_geomParams->Get_positionBarrelLayer()[(startPseudoLayer+layersFromStart)]-_geomParams->Get_positionBarrelLayer()[startPseudoLayer];
      scale = fabs(z_dist/clusterDirection->z);
      endpoint.x = myPos->x+(clusterDirection->x*scale);
      endpoint.y = myPos->y+(clusterDirection->y*scale);
      endpoint.z = myPos->z+(clusterDirection->z*scale);
    }
    
    break;
  }
  case CLUS_LOCATION_UNKNOWN:  {
    cout << "ERROR: Cluster location is unknown!!!" << endl;
    break;
  }
    // TODO: treat overlap
  }
  if(_algoParams->GetDebug()>2) 
    cout << "Followed cluster direction to (" << endpoint.x << ", " << endpoint.y << ", " << endpoint.z << ")" << endl;
  return endpoint;
}


void ECALGarlicClusterHelpers::FreeHits(ExtendedCluster *noCluster)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  vector<ExtendedHit* > *clusterHits = &(noCluster->hitVec);
  int NHitsInCluster = clusterHits->size();
  
  if(NHitsInCluster>0) {
    for( int hit_i = 0; hit_i<NHitsInCluster; hit_i++ ) {
      ExtendedHit *a_ext_hit=dynamic_cast<ExtendedHit *>((*clusterHits)[hit_i]);
      if(a_ext_hit) { 
	a_ext_hit->cluster=0;
	a_ext_hit->clusterHitVec=0;
      }
    }
  }
  clusterHits->clear();
  if (_algoParams->GetDebug()>2) cout << NHitsInCluster << " hits freed in cluster, cluster deleted" << endl;
}


void ECALGarlicClusterHelpers::FreeHits(vector<ExtendedHit*> &noCluster)
{
  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();
  
  int NHitsInCluster=noCluster.size();
  for( int hit_i = 0; hit_i<NHitsInCluster; hit_i++ ) {
    ExtendedHit *a_ext_hit=dynamic_cast<ExtendedHit *>(noCluster[hit_i]);
    if(a_ext_hit) { 
      a_ext_hit->cluster=0;
      a_ext_hit->clusterHitVec=0;
    }
  }
  noCluster.clear();
  if(_algoParams->GetDebug()>2) 
    cout << NHitsInCluster << " hits freed in cluster" << endl;
}

bool ECALGarlicClusterHelpers::CheckClusterCriteria(ExtendedCluster *myCluster,vector<ExtendedTrack*> tracks )
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  bool veto=false;
  
  vector<ExtendedHit* > *clusterHits = &(myCluster->hitVec);
  int NClusteredHits = clusterHits->size();

  double Etot=0;
  int has_b_hits = 0;
  int has_e_hits = 0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    Etot+=hit_en;
    CalorimeterHitZone hit_zone = (CalorimeterHitZone)myHit->zone;
    if(hit_zone == CALHITZONE_BARREL)
       has_b_hits = 1;
    else 
      if(hit_zone == CALHITZONE_ENDCAP || hit_zone == CALHITZONE_RING || hit_zone == CALHITZONE_OVERLAP )
	has_e_hits = 1;
  }
  double E_GeV = 0;
  //double E_GeV = Etot;
  //   modified for initial finiding
  //  if(has_e_hits==1 && has_b_hits==0 ) {
  //    //E_GeV = toGeVFctn_EC->Eval(Etot);
  //    E_GeV = ECALGarlicEnergyEstimator::getGeV_EC(Etot);
  //  } else {
  //    //    E_GeV = toGeVFctn->Eval(Etot);
  //    E_GeV = ECALGarlicEnergyEstimator::getGeV(Etot);
  //  }  

  bool isbarrel = ! (has_e_hits==1 && has_b_hits==0);
  E_GeV = _energyEstimator->getGeV(Etot, isbarrel);

  // total energy has to bigger than _minEnergy
  if(E_GeV<_algoParams->GetMinEnergy())
    return 0;

  double shortest_dist_track=9999;

  // distance from track
  int nTracks = tracks.size();
  for(int t_i=0;t_i<nTracks;t_i++) {
    ExtendedTrack *a_track = dynamic_cast<ExtendedTrack*> (tracks[t_i]);
    HelixClass *a_helix = a_track->helix;
    float dist[3];
    for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
      ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
      float hitPos[3];
      hitPos[0]=myHit->hit->getPosition()[0];
      hitPos[1]=myHit->hit->getPosition()[1];
      hitPos[2]=myHit->hit->getPosition()[2];
      float dummy;
      dummy = a_helix->getDistanceToPoint(hitPos,dist);
      if(dist[2]<shortest_dist_track)
	shortest_dist_track=dist[2];
    }
  }
  //if(shortest_dist_track<(3*_geomParams->Get_padSizeEcal() [1])) {
  double allowed_dist = sqrt((0.5*_geomParams->Get_padSizeEcal() [1]*0.5*_geomParams->Get_padSizeEcal() [1])+(1.5*_geomParams->Get_padSizeEcal() [1]*1.5*_geomParams->Get_padSizeEcal() [1]));
  if(shortest_dist_track<(allowed_dist-0.01)) { // changed 11/03/09 
    cout << "Cluster too close to track, E= " << E_GeV << endl;
    return 0;
  }
  if(veto==true)
    return 0;
  if (_algoParams->GetDebug()>1) cout << "Cluster is real, E=" << E_GeV << endl;
  return 1;
}



void ECALGarlicClusterHelpers::GetFreeNeighbours(ExtendedHit *a_hit, vector<ExtendedHit* > &myNeighbours, ExtendedCluster *preCluster, int iteration, double Theta)
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  vector<ExtendedHit* > *hits=&(preCluster->hitVec);
  int search_pl=a_hit->pseudoLayer;
  int search_layer=(*(_geomParams->Get_defaultDecoder()))(a_hit->hit)["K-1"];
  int search_i=(*(_geomParams->Get_defaultDecoder()))(a_hit->hit)["I"];
  int search_j=(*(_geomParams->Get_defaultDecoder()))(a_hit->hit)["J"];
  int search_m=(*(_geomParams->Get_defaultDecoder()))(a_hit->hit)["M"];
  int search_s=(*(_geomParams->Get_defaultDecoder()))(a_hit->hit)["S-1"];
  int search_z=a_hit->zone;
  vec3 hitPos;
  hitPos.x=a_hit->hit->getPosition()[0];
  hitPos.y=a_hit->hit->getPosition()[1];
  hitPos.z=a_hit->hit->getPosition()[2];
  int NHitsInPreCluster=hits->size();
  for( int hit_i = 0; hit_i<NHitsInPreCluster; hit_i++ ) {
    ExtendedHit *a_ext_hit=dynamic_cast<ExtendedHit *>((*hits)[hit_i]);
    if(a_ext_hit) {
      if(a_ext_hit->clusterHitVec!=0)
	continue;
      if(a_ext_hit->pseudoLayer<search_pl+4 && a_ext_hit->pseudoLayer>=search_pl) {
	vec3 nHitPos;
	nHitPos.x=a_ext_hit->hit->getPosition()[0];
	nHitPos.y=a_ext_hit->hit->getPosition()[1];
	nHitPos.z=a_ext_hit->hit->getPosition()[2];
	double dist=9999;

	int cell_i=(*(_geomParams->Get_defaultDecoder()))(a_ext_hit->hit)["I"];
	int cell_j=(*(_geomParams->Get_defaultDecoder()))(a_ext_hit->hit)["J"];
	int cell_m=(*(_geomParams->Get_defaultDecoder()))(a_ext_hit->hit)["M"];
	int cell_s=(*(_geomParams->Get_defaultDecoder()))(a_ext_hit->hit)["S-1"];
	int cell_z=a_ext_hit->zone;
	double dist3d = _geomHelper->Get3dDistance(&hitPos,&nHitPos);
	double dist2dProj= _geomHelper->Get2dProjDistance(&hitPos,&nHitPos);
	if(dist3d>(_geomParams->Get_zOfEndcap()-_geomParams->Get_zOfBarrel())) // eliminates possibility of barrel-endcap comparison
	  continue;
	dist=dist2dProj;

	if(search_m<6 && search_m>0) {//barrel solution to limit search radius
	  if(cell_m==search_m) { //same module...
	    if(cell_s==search_s) { //...and the same stave
	      if( (fabs(float(cell_i-search_i))>10) || (fabs(float(cell_j-search_j))>10) )
		continue;

	      double maxdist = (((_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer]+2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize())+0.01);
	      if(dist<(maxdist/sin(Theta)))
		myNeighbours.push_back(a_ext_hit);
	    }
	    else { // same module but not the same stave
	      double maxdist = (_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer];
	      if(dist<(maxdist/sin(Theta)))
		myNeighbours.push_back(a_ext_hit);
	    }
	  }
	  else { // two different modules...
	    if(fabs(float(cell_m-search_m))>1) // too far
	      continue;
	    if(cell_s==search_s) { // ...but the same stave
	      if(cell_m<search_m) {
		int new_cell_j = (cell_j - 10*_geomParams->Get_nCellsPerWafer());
		if((fabs(float(new_cell_j-search_j))>10) )
		  continue;
	      }
	      else {
		int new_search_j = (search_j - 10*_geomParams->Get_nCellsPerWafer());
		if((fabs(float(cell_j-new_search_j))>10) )
		  continue;
	      }
	      
	      double maxdist = (((_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])+2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule())+0.01;
	      if(dist<(maxdist/sin(Theta))) {
		myNeighbours.push_back(a_ext_hit);
	      }
	    }
	    else { // in different modules and different staves
	      if(cell_m<search_m) {
		int new_cell_j = (cell_j - 10*_geomParams->Get_nCellsPerWafer());
		if((fabs(float(new_cell_j-search_j))>10) )
		  continue;
	      }
	      else {
		int new_search_j = (search_j - 10*_geomParams->Get_nCellsPerWafer());
		if((fabs(float(cell_j-new_search_j))>10) )
		  continue;
	      }
	      double maxdist = (((_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])+2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule())+0.01;
	      if(dist<(maxdist/sin(Theta)))
		myNeighbours.push_back(a_ext_hit);
	    }
	  }
	}
	else { //endcap solution
	  if(cell_m!=search_m)  // avoid +-z comparisons
	    continue;
	  if(search_z!=CALHITZONE_RING && cell_z!=CALHITZONE_RING) { // neglect EcalRing
	    if(cell_s==search_s) { //...the same stave
	      if( (fabs(float(cell_i-search_i))>10) || (fabs(float(cell_j-search_j))>10) )
		continue;

	      double maxdist = (((_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])+2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize())+0.01;
	      if(dist<(maxdist/fabs(cos(Theta))))
		myNeighbours.push_back(a_ext_hit);
	    }
	    else { // not the same stave

	      if(!( ( cell_i<5 && search_j<5) || ( cell_j<5 && search_i>(10*_geomParams->Get_nCellsPerWafer()-6) ) ))
		continue;	     

	      double maxdist=sqrt((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])*
				  (2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer]) +
				  (_geomParams->Get_padSizeEcal() [search_layer]*_geomParams->Get_padSizeEcal() [search_layer]));
	      if(dist<maxdist/fabs(cos(Theta))) {
		myNeighbours.push_back(a_ext_hit);
	      }
	    }
	  }
	  else {  // now EcalRing
	    if(cell_z==CALHITZONE_ENDCAP && search_z==CALHITZONE_RING) {
	      if(cell_i!=0)
		continue;
	      double maxdist=sqrt((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])*
				  (2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer]) +
				  (_geomParams->Get_padSizeEcal() [search_layer]*_geomParams->Get_padSizeEcal() [search_layer]));
	      if(dist<(maxdist/fabs(cos(Theta)))) {
		myNeighbours.push_back(a_ext_hit);
	      }
	    }
	    if(cell_z==CALHITZONE_RING && search_z==CALHITZONE_ENDCAP) {
	      if(search_i!=0)
		continue;
	      double maxdist=sqrt((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer])*
				  (2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+(_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer]) +
				  (_geomParams->Get_padSizeEcal() [search_layer]*_geomParams->Get_padSizeEcal() [search_layer]));
	      if(dist<(maxdist/fabs(cos(Theta)))) {
		myNeighbours.push_back(a_ext_hit);
	      }
	    }
	    if(cell_z==CALHITZONE_RING && search_z==CALHITZONE_RING) {
	      if( (fabs(float(cell_i-search_i))>10) || (fabs(float(cell_j-search_j))>10) )
		continue;

	      double maxdist = ((_algoParams->GetDistanceFactor())*_geomParams->Get_padSizeEcal() [search_layer]+2*_geomParams->Get_guardringSize())+0.01;
	      if(dist<(maxdist/fabs(cos(Theta)))) {
		  myNeighbours.push_back(a_ext_hit);
	      }
	    }
	  }
	}
      }
    }
  }
  if(_algoParams->GetDebug()>2) 
    cout << "Found " << myNeighbours.size() << " neighbours " << endl;
}




void ECALGarlicClusterHelpers::RecAddNeighbourBins(TH2F *histo, int bin_i, vector<int> *seed_c) 
{

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  int nX  = histo->GetXaxis()->GetNbins()+2;
  int bin_x = bin_i%nX;
  int bin_y = bin_i/nX;
  int neigh_bins[4];
  neigh_bins[0] = histo->GetBin(bin_x-1,bin_y);
  neigh_bins[1] = histo->GetBin(bin_x+1,bin_y);
  neigh_bins[2] = histo->GetBin(bin_x,bin_y+1);
  neigh_bins[3] = histo->GetBin(bin_x,bin_y-1);
  for(int n_i=0;n_i<4;n_i++) {
    int n_bin_i = neigh_bins[n_i];
    float bin_cont = histo->GetBinContent(n_bin_i);
    if(bin_cont>0 && bin_cont < (histo->GetBinContent(bin_i))) {
      seed_c->push_back(n_bin_i);
      RecAddNeighbourBins(histo,n_bin_i,seed_c);
    }
  }
}


void ECALGarlicClusterHelpers::IterAddNeighbourBins(TH2F *histo, int bin_i, vector<int> *seed_c , int iteration) 
{  

  if (!_nn_is_setup && _geomParams->Get_symmetry()>0) setupNN();

  if(iteration==5)
    return;
  bool force = false;
  if(iteration<3)
    force=true;

  int nX  = histo->GetXaxis()->GetNbins()+2;
  int bin_x = bin_i%nX;
  int bin_y = bin_i/nX;
  int neigh_bins[4];
  neigh_bins[0] = histo->GetBin(bin_x-1,bin_y);
  neigh_bins[1] = histo->GetBin(bin_x+1,bin_y);
  neigh_bins[2] = histo->GetBin(bin_x,bin_y+1);
  neigh_bins[3] = histo->GetBin(bin_x,bin_y-1);
  for(int n_i=0;n_i<4;n_i++) {
    int n_bin_i = neigh_bins[n_i];
    bool found = false;
    for(unsigned int j=0;j<seed_c->size();j++) {
      if(n_bin_i==(*seed_c)[j]) {
	found=true;
	break;
      }
    }
    if(!found) {      
      float bin_cont = histo->GetBinContent(n_bin_i);
      if(bin_cont==0)
	continue;
      if(!force)
	if(bin_cont > (histo->GetBinContent(bin_i))) 
	  continue;
      seed_c->push_back(n_bin_i);
      IterAddNeighbourBins(histo,n_bin_i,seed_c,iteration+1);
    }
  }
}
