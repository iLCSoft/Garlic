#include "ECALGarlicGeometryHelpers.hh"

#include <assert.h>
#include <cmath>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

double ECALGarlicGeometryHelpers::Get2dProjDistance(vec3 *a,vec3 *b)
{

  assert (_geomParams->Get_zOfBarrel()>-999 && _geomParams->Get_rOfBarrel()>-999 && _geomParams->Get_zOfEndcap()>-999);

  double dist;
  if((a->z<_geomParams->Get_zOfBarrel() && a->z>(-_geomParams->Get_zOfBarrel()) && (b->z<_geomParams->Get_zOfBarrel()) && b->z>(-_geomParams->Get_zOfBarrel()))) {
    double xp_a=a->x/sqrt(a->x*a->x+a->y*a->y)*_geomParams->Get_rOfBarrel();
    double yp_a=a->y/sqrt(a->x*a->x+a->y*a->y)*_geomParams->Get_rOfBarrel();
    double zp_a=a->z/sqrt(a->x*a->x+a->y*a->y)*_geomParams->Get_rOfBarrel();
    double xp_b=b->x/sqrt(b->x*b->x+b->y*b->y)*_geomParams->Get_rOfBarrel();
    double yp_b=b->y/sqrt(b->x*b->x+b->y*b->y)*_geomParams->Get_rOfBarrel();
    double zp_b=b->z/sqrt(b->x*b->x+b->y*b->y)*_geomParams->Get_rOfBarrel();
    dist=sqrt((xp_a-xp_b)*(xp_a-xp_b)+(yp_a-yp_b)*(yp_a-yp_b)+(zp_a-zp_b)*(zp_a-zp_b));
  }
  else {
    double scale_a = fabs(_geomParams->Get_zOfEndcap()/a->z);
    double x_a=(a->x)*scale_a;
    double y_a=(a->y)*scale_a;
    
    double scale_b = fabs(_geomParams->Get_zOfEndcap()/b->z);
    double x_b=(b->x)*scale_b;
    double y_b=(b->y)*scale_b;
    dist=sqrt((x_a-x_b)*(x_a-x_b)+(y_a-y_b)*(y_a-y_b));
  }
  return dist;
}


double ECALGarlicGeometryHelpers::Get3dDistance(vec3 *a,vec3 *b)
{
  double x_a=a->x;
  double y_a=a->y;
  double z_a=a->z;
  double x_b=b->x;
  double y_b=b->y;
  double z_b=b->z;
  double dist=sqrt((x_a-x_b)*(x_a-x_b)+(y_a-y_b)*(y_a-y_b)+(z_a-z_b)*(z_a-z_b));
  return dist;
}

void ECALGarlicGeometryHelpers::GetDistancesBetweenClusters(ExtendedCluster *a,ExtendedCluster *b, double *distances)
{
  ClusterParameters *a_par = a->Parameters;
  vec3 a_cog;
  a_cog.x = a_par->COGx;
  a_cog.y = a_par->COGy;
  a_cog.z = a_par->COGz;
  ClusterParameters *b_par = b->Parameters;
  vec3 b_cog;
  b_cog.x = b_par->COGx;
  b_cog.y = b_par->COGy;
  b_cog.z = b_par->COGz;
  double smallest_dist=9999;
  double firstHitDistance=9999;
  int firstBPSLayer=99;
  int firstAPSLayer=99;
  std::vector<ExtendedHit* > *aHits = &(a->hitVec);
  int NaHits = aHits->size();
  vector<ExtendedHit* > *bHits = &(b->hitVec);
  int NbHits = bHits->size();
  for(int b_hit_i=0;b_hit_i<NbHits;b_hit_i++) {
    ExtendedHit *myBHit = dynamic_cast<ExtendedHit* > ((*bHits)[b_hit_i]);
    vec3 b_pos;
    b_pos.x=(myBHit->hit)->getPosition()[0];
    b_pos.y=(myBHit->hit)->getPosition()[1];
    b_pos.z=(myBHit->hit)->getPosition()[2];
    for(int a_hit_i=0;a_hit_i<NaHits;a_hit_i++) {
      ExtendedHit *myAHit = dynamic_cast<ExtendedHit* > ((*aHits)[a_hit_i]);
      vec3 a_pos;
      a_pos.x=(myAHit->hit)->getPosition()[0];
      a_pos.y=(myAHit->hit)->getPosition()[1];
      a_pos.z=(myAHit->hit)->getPosition()[2];
      double dist=ECALGarlicGeometryHelpers::Get3dDistance(&a_pos,&b_pos);
      if(dist<smallest_dist)
	smallest_dist=dist;
      if(myBHit->pseudoLayer<=firstBPSLayer) {
	firstBPSLayer=myBHit->pseudoLayer;
	if(myAHit->pseudoLayer<=firstAPSLayer) {
	  firstAPSLayer=myAHit->pseudoLayer;
	  double a_firstDistance = ECALGarlicGeometryHelpers::Get2dProjDistance(&a_pos,&b_pos);
	  if(a_firstDistance<firstHitDistance) {
	    firstHitDistance=a_firstDistance;
	  }
	}
      }
    }
  }
  distances[0] = smallest_dist;
  distances[1] = ECALGarlicGeometryHelpers::Get3dDistance(&a_cog,&b_cog);
  distances[2] = firstHitDistance;
}



void ECALGarlicGeometryHelpers::AssignPseudoLayer(ExtendedHit* &a_hit)
{
  //again do the same as Pandora except endcap, here psLayer=physical layer
  const float* pos =(a_hit->hit)->getPosition();
  
  CalorimeterHitZone zone = CALHITZONE_UNKNOWN;
  double prodRadiusMax=0.;
 
  for(int istave = 0; istave < _geomParams->Get_symmetry(); ++istave){
    double prodRadius = pos[0]*_geomParams->Get_barrelStaveDir()[istave].x+pos[1]*_geomParams->Get_barrelStaveDir()[istave].y;
    if(prodRadius>prodRadiusMax)prodRadiusMax=prodRadius;
  }
  float xhit = fabs(pos[0]);
  float yhit = fabs(pos[1]);
  float zhit = fabs(pos[2]);

  if(zhit>_geomParams->Get_zOfBarrel())zone = CALHITZONE_ENDCAP;
  if(zhit<_geomParams->Get_positionEndcapLayer()[0])zone = CALHITZONE_BARREL;
  if(zhit>_geomParams->Get_zOfBarrel() && xhit<_geomParams->Get_rInnerEcalEndcap() && yhit<_geomParams->Get_rInnerEcalEndcap()){zone = CALHITZONE_RING;}
  if(zone==CALHITZONE_UNKNOWN)zone = CALHITZONE_OVERLAP;

  a_hit->zone=zone;

  int   bestLayer=0;
  int   bestBarrelLayer=0;
  double barrelLayerDist = 9999.;
  int   bestEndcapLayer=0;
  double endcapLayerDist = 9999.;

  switch(zone){
  case CALHITZONE_BARREL: 
    // use barrel layer 
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
	barrelLayerDist = fabs(prodRadiusMax-_geomParams->Get_positionBarrelLayer()[ilayer]);
	bestBarrelLayer = ilayer;
      }
    }
    a_hit->pseudoLayer=bestBarrelLayer;
    break;
  case CALHITZONE_ENDCAP: 
    // use end layer 
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    a_hit->pseudoLayer=bestEndcapLayer;
    break;
  case CALHITZONE_RING: // use the same convention as for the endcap 
    // use end layer 
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    a_hit->pseudoLayer=bestEndcapLayer;
    break;
  case CALHITZONE_OVERLAP: // NEW: do not assign barrel ps layers
    for(int ilayer = 0; ilayer <= _geomParams->Get_nPseudoLayers(); ++ilayer){
      if(fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer])<endcapLayerDist){
	endcapLayerDist = fabs(zhit-_geomParams->Get_positionEndcapLayer()[ilayer]);
	bestEndcapLayer = ilayer;
      }
    }
    bestLayer = bestEndcapLayer;
    a_hit->pseudoLayer=bestLayer;
    break;
  default:
    std::cout << " ERROR SHOULD NOT GET HERE " << std::endl;
    break;
  }
}


void ECALGarlicGeometryHelpers::GetNearestHitInPseudoLayer(int ps_layer,vec3 *seed,ExtendedCluster *preCluster, ExtendedHit* &nearestHit)
{
  int NHitsInCluster=(preCluster->hitVec).size();
  double minDist = 9999;
  vector<ExtendedHit* > &hits=preCluster->hitVec;
  //nearestHit = NULL;
  for( int hit_i = 0; hit_i<NHitsInCluster; hit_i++ ) {
    ExtendedHit *a_ext_hit= hits[hit_i];
    if(a_ext_hit) {
      if(a_ext_hit->clusterHitVec!=0 || a_ext_hit->pseudoLayer!=ps_layer) continue;
      vec3 hitPos;
      hitPos.x=a_ext_hit->hit->getPosition()[0];
      hitPos.y=a_ext_hit->hit->getPosition()[1];
      hitPos.z=a_ext_hit->hit->getPosition()[2];
      double dist=ECALGarlicGeometryHelpers::Get2dProjDistance(&hitPos,seed);
      if(dist<minDist && (fabs(hitPos.z-seed->z)<(_geomParams->Get_zOfEndcap()-_geomParams->Get_zOfBarrel())))  {
	minDist=dist;
	nearestHit=a_ext_hit;
      }
    }
  }
  if(nearestHit==0) {
    if(_algoParams->GetDebug()>2)
      cout << "No hit found in layer " << ps_layer << endl;
  }
  else {
    if(_algoParams->GetDebug()>2) {
      vec3 hitPos;
      hitPos.x=nearestHit->hit->getPosition()[0];
      hitPos.y=nearestHit->hit->getPosition()[1];
      hitPos.z=nearestHit->hit->getPosition()[2];
      cout << "Nearest Hit is " << nearestHit->hit->getCellID0() << " at (" << hitPos.x << ", " << hitPos.y << ", " << hitPos.z << ") in pseudo layer " << ps_layer << endl;
    }
  }
}

