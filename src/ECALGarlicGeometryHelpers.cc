#include "ECALGarlicGeometryHelpers.hh"

#include <assert.h>
#include <cmath>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

#include "ECALGarlicGeometryParameters.hh"

void ECALGarlicGeometryHelpers::getFrontFaceProj(const float* a, float* proj) {
  // project a point towards/away from the origin (IP) onto the ECAL front face
  float scale(0);

  float scaleBarrel = ECALGarlicGeometryParameters::Instance().Get_rOfBarrel()/ sqrt(pow(a[0], 2) + pow(a[1], 2));
  float scaleEndcap = fabs(ECALGarlicGeometryParameters::Instance().Get_zOfEndcap()/a[2]);

  if ( fabs ( scaleBarrel*a[2] ) < ECALGarlicGeometryParameters::Instance().Get_zOfEndcap() )  // scaled position in barrel region
    scale = scaleBarrel;
  else
    scale = scaleEndcap;

  for (int i=0; i<3; i++) proj[i]=scale*a[i];

  // // check my new fn
  // float pp[3];
  // ECALGarlicGeometryHelpers::getFrontFaceProj(a, a, pp);
  // 
  // if (  fabs( pp[0]-proj[0] ) > 0.1 ||
  // 	fabs( pp[1]-proj[1] ) > 0.1 ||
  // 	fabs( pp[2]-proj[2] ) > 0.1 ) {
  // 
  //   cout << "check proj: ";
  //   for (int i=0; i<3; i++) 
  //     cout << a[i] << " ";
  //   cout << ": ";
  //   for (int i=0; i<3; i++) 
  //     cout << proj[i] << " ";
  //   cout << " : ";
  //   for (int i=0; i<3; i++) cout << pp[i] << " ";
  //   cout << endl;
  // 
  // }

  return;
}


void ECALGarlicGeometryHelpers::getLineCylinderIntersection( const float* aPoint, const float* aDir, float radius, float* intersec1, float* intersec2 ) {
  // intersection of 3d line with infinite cylinder aligned along z axis, axis at x,y=0,0

  // first get 2-d intersections
  getLineCircleIntersection( aPoint, aDir, radius, intersec1, intersec2 );

  // calculate z position at each x,y point
  if ( fabs(intersec1[0])>0 || fabs(intersec1[1])>0 ) { // found intersection in 2d
    if ( fabs(aDir[0])>0 ) {
      intersec1[2] = aPoint[2] + (aDir[2]/aDir[0])*(intersec1[0]-aPoint[0]);
      intersec2[2] = aPoint[2] + (aDir[2]/aDir[0])*(intersec2[0]-aPoint[0]);
    } else {
      intersec1[2] = aPoint[2] + (aDir[2]/aDir[1])*(intersec1[1]-aPoint[1]);
      intersec2[2] = aPoint[2] + (aDir[2]/aDir[1])*(intersec2[1]-aPoint[1]);
    }
  } else { // no intersection
    intersec1[2]=0;
    intersec2[2]=0;
  }

  //  cout << " zpositions : " << intersec1[2] << " " << intersec2[2] << endl;

  return;
}

void ECALGarlicGeometryHelpers::getLineCircleIntersection( const float* aPoint, const float* aDir, float radius, float* intersec1, float* intersec2 ) {
  // circle centered at 0,0
  // 2-d
  
  // eqn of line: y = a + b*x
  float a = aPoint[1] - aPoint[0]*aDir[1]/aDir[0];
  float b = aDir[1]/aDir[0];

  // then solve for intersection:
  // x**2 + y**2 = r**2
  // y = a + b*x
  // x = ( -2ab +- sqrt ( 4 a**2 b**2 - 4*(1+b**2)(a**2-r**2) ) ) / 2(1+b**2)

  float det = pow( 2.*a*b, 2.) - 4 * ( 1. + pow(b, 2) ) * ( pow(a, 2.) - pow(radius, 2.) );
  if ( det<0 ) { // no intersection
    intersec1[0]=0;
    intersec1[1]=0;

    intersec2[0]=0;
    intersec2[1]=0;
  } else {
    det=sqrt(det);
    intersec1[0] = ( -2.*a*b + det ) / ( 2.*(1.+pow(b, 2.) ) );
    intersec1[1] = a + b*intersec1[0];
    
    intersec2[0] = ( -2.*a*b - det ) / ( 2.*(1.+pow(b, 2.) ) );
    intersec2[1] = a + b*intersec2[0];
  }

  //  cout << "2d intersects: " << intersec1[0] << " " << intersec1[1] << " : " << intersec2[0] << " " << intersec2[1] << endl;

  return;
}


void ECALGarlicGeometryHelpers::getFrontFaceProj(const float* a, const float* dir, float* proj) {

  // cout << "point : ";for (int i=0; i<3; i++) cout << a[i] << " "; 
  // cout << " : dir "; for (int i=0; i<3; i++) cout << dir[i] << " "; 
  // cout << endl;
  // cout << sqrt( a[0]*a[0] + a[1]*a[1] ) << " - " << ECALGarlicGeometryParameters::Instance().Get_rOfBarrel() << endl;

  // project a point along some direction onto the ECAL front face

  // first try the endcaps (easier)
  float pointOnPlane[3]={0, 0, fabs(ECALGarlicGeometryParameters::Instance().Get_zOfEndcap()) };
  if ( a[2]<0 ) pointOnPlane[2]*=-1;
  float planeNorm[3]={0,0,1};
  float endcapProj[3];
  GetGeneralPointPlaneProjection(a, pointOnPlane, planeNorm, dir, endcapProj);
  // check the radius is in the correct range
  float rr = sqrt ( pow( endcapProj[0], 2) + pow( endcapProj[1], 2) );
  bool intersectEndcap = rr <= ECALGarlicGeometryParameters::Instance().Get_rOfBarrel();
   
  // now the barrel
  float barrelIntersection1[3];
  float barrelIntersection2[3];
  getLineCylinderIntersection( a, dir, ECALGarlicGeometryParameters::Instance().Get_rOfBarrel(), barrelIntersection1, barrelIntersection2 );
  // check the barrel intersections are in the correct z range

  bool intersectBarrel1=false;
  bool intersectBarrel2=false;
  for (int i=0; i<3; i++) {
    if ( fabs(barrelIntersection1[i])>0.1 ) { intersectBarrel1=true; break;}
  }
  for (int i=0; i<3; i++) {
    if ( fabs(barrelIntersection2[i])>0.1 ) { intersectBarrel2=true; break;}
  }

  //  cout << " mmm " << intersectBarrel1 << " " << intersectBarrel2 << endl;

  //  bool atLeatOneBarrelIntersection = intersectBarrel1 || intersectBarrel2;
  //  if ( ! atLeatOneBarrelIntersection ) {
  // check that it doesn;t just miss...
  // find closest point approach to axis

  intersectBarrel1 = intersectBarrel1 && fabs ( barrelIntersection1[2] ) <= fabs(ECALGarlicGeometryParameters::Instance().Get_zOfEndcap());
  intersectBarrel2 = intersectBarrel2 && fabs ( barrelIntersection2[2] ) <= fabs(ECALGarlicGeometryParameters::Instance().Get_zOfEndcap());

//  cout << "point " << a[0] << " " << a[1] << " " << a[2] << endl;
//
//  cout << "endcap : " << intersectEndcap << " "; 
//  for (int i=0; i<3; i++) cout << endcapProj[i] << " "; cout << endl;
//  
//  cout << "barrel1: " << intersectBarrel1 << " "; 
//  for (int i=0; i<3; i++) cout << barrelIntersection1[i] << " "; 
//  cout << sqrt( barrelIntersection1[0]*barrelIntersection1[0] + barrelIntersection1[1]*barrelIntersection1[1] ) << endl;
//  
//  cout << "barrel2: " << intersectBarrel2 << " "; 
//  for (int i=0; i<3; i++) cout << barrelIntersection2[i] << " "; 
//  cout << sqrt( barrelIntersection2[0]*barrelIntersection2[0] + barrelIntersection2[1]*barrelIntersection2[1] ) << endl;


  // now choose the remaining intersection which is closest to the original point
  if ( !intersectEndcap && !intersectBarrel1 && !intersectBarrel2 ) {
    // no intersection
    for (int i=0; i<3; i++) proj[i]=0;
  } else {
    float dist[3]={0,0,0};
    if ( !intersectEndcap ) dist[0]=999999;
    if ( !intersectBarrel1 ) dist[1]=999999;
    if ( !intersectBarrel2 ) dist[2]=999999;

    for (int i=0; i<3; i++) {
      if ( intersectEndcap )  dist[0]+=pow( endcapProj[i]-a[i], 2);
      if ( intersectBarrel1 ) dist[1]+=pow( barrelIntersection1[i]-a[i], 2);
      if ( intersectBarrel2 ) dist[2]+=pow( barrelIntersection2[i]-a[i], 2);
    }

    if      ( dist[0]<=dist[1] && dist[0]<=dist[2] ) 
      for (int i=0; i<3; i++) proj[i] = endcapProj[i];
    else if ( dist[1]<=dist[0] && dist[1]<=dist[2] ) 
      for (int i=0; i<3; i++) proj[i] = barrelIntersection1[i];
    else if ( dist[2]<=dist[0] && dist[2]<=dist[1] ) 
      for (int i=0; i<3; i++) proj[i] = barrelIntersection2[i];
    else 
      cout << "ECALGarlicGeometryHelpers::getFrontFaceProj shouldn;t get here? " << dist[0] << " " << dist[1] << " " << dist[2] << endl;
  }
  
  return;
}

double ECALGarlicGeometryHelpers::Get2dProjDistance(const float* a, const float* b) {
  float aProj[3];
  getFrontFaceProj(a, aProj);
  float bProj[3];
  getFrontFaceProj(b, bProj);
  return Get3dDistance(aProj, bProj);
}

void ECALGarlicGeometryHelpers::GetGeneralPointPlaneProjection(const float* point, const float* pointOnPlane, const float* planeNormal, const float* projDirection, float* projection) {

  // project a general point along a general projDirection onto a plane defined by pointOnPlane and planeNormal
  float temp[3];
  diff(pointOnPlane, point, temp);

  float alpha = dot( temp, planeNormal ) / dot( projDirection, planeNormal );

  float temp2[3];
  scale( projDirection, alpha, temp2 );

  sum( point, temp2, projection );
  return;
}


double ECALGarlicGeometryHelpers::Get3dDistance(const float* a, const float* b) {
  return sqrt ( pow(a[0] - b[0],2) + pow(a[1] - b[1],2) + pow(a[2] - b[2],2) );
}

double ECALGarlicGeometryHelpers::GetDistToLine(const float* a, const float* l1, const float* l2) {
  // point a, line between 11, 12
  // dist = | (a-l1) x (a-l2) | / | (l2-l1) |
  float d1[3];
  float d2[3];
  float d3[3];
  float cr[3];
  diff(a, l1, d1);
  diff(a, l2, d2);
  diff(l2, l1, d3);
  cross (d1, d2, cr);
  return mag(cr)/mag (d3);
}

void ECALGarlicGeometryHelpers::cross(const float* a, const float* b, float* vec) {
  vec[0] = a[1]*b[2] - a[2]*b[1];
  vec[1] = a[2]*b[0] - a[0]*b[2];
  vec[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

float ECALGarlicGeometryHelpers::dot(const float* a, const float* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float ECALGarlicGeometryHelpers::mag(const float* a) {
  return sqrt( pow(a[0], 2)+pow(a[1], 2)+pow(a[2], 2) );
}

void ECALGarlicGeometryHelpers::diff(const float* a, const float* b, float* vec) {
  for (int i=0; i<3; i++) vec[i] = a[i]-b[i];
  return;
}

void ECALGarlicGeometryHelpers::sum(const float* a, const float* b, float* vec) {
  for (int i=0; i<3; i++) vec[i] = a[i]+b[i];
  return;
}

float ECALGarlicGeometryHelpers::cosangle(const float* a, const float* b) {
  return dot(a,b)/(mag(a)*mag(b));
}

float ECALGarlicGeometryHelpers::angle(const float* a, const float* b) {
  return acos(cosangle(a, b));
}

float ECALGarlicGeometryHelpers::phi(const float* a) {
  float ph = atan2(a[1],a[0]);
  if (ph<0) ph+=2*acos(-1);
  return ph;
}

float ECALGarlicGeometryHelpers::costheta(const float* a) {
  return a[2]/mag(a);
}

float ECALGarlicGeometryHelpers::theta(const float* a) {
  return acos(costheta(a));
}

float ECALGarlicGeometryHelpers::getPseudoLayerRadius(int psLayer) {
  return ECALGarlicGeometryParameters::Instance().Get_positionBarrelLayer() [psLayer];
}


float ECALGarlicGeometryHelpers::getPseudoLayerDistNormal(int ps1, int ps2) {
  return ECALGarlicGeometryParameters::Instance().Get_positionBarrelLayer()[ps2] - ECALGarlicGeometryParameters::Instance().Get_positionBarrelLayer()[ps2];
}

int ECALGarlicGeometryHelpers::getPseudoStave(const float* pos) {
  return int( (phi(pos) + acos(-1)/ECALGarlicGeometryParameters::Instance().Get_symmetry()) / (2*acos(-1)/ECALGarlicGeometryParameters::Instance().Get_symmetry()) );
}

float ECALGarlicGeometryHelpers::getPseudoStavePhi(int istave) {
  return istave*(2*acos(-1)/ECALGarlicGeometryParameters::Instance().Get_symmetry());
}

void ECALGarlicGeometryHelpers::norm(const float* a, float* norma) {
  scale(a, mag(a), norma);
}

void ECALGarlicGeometryHelpers::scale(const float* a, float factor, float* aprime) {
  for (int i=0; i<3; i++) aprime[i]=a[i]*factor;
}

float ECALGarlicGeometryHelpers::det(const float* abcd) {
  /*

  calculate determinant of (  a   b  ) = ad-cb
                           (  c   d  )

  */
  return abcd[0]*abcd[3] - abcd[2]*abcd[1];
}

void ECALGarlicGeometryHelpers::getLineBarrelPseudoLayerIntersection(const float* point1, const float* point2, int pseudolayer, float* intersection) {
  // get the intersection of 2d line (in xy, linking point1, point2)
  //  with a particular (barrel) pseudolayer

  // start looking in pseudostave containing first point on line
  int pstave = getPseudoStave(point1);

  float pt1[2]={0};
  float pt2[2]={0};
  float intsec[2];
  bool gotit(false);
  for (int dstave=0; dstave<=4; dstave++) { // start looking in same stave, then look further away
    for (int ipm=0; ipm<2; ipm++) { // look on either side of initial stave
      if ( (dstave==0 || dstave==4) && ipm>0) continue;      
      int stave = ipm==0 ? pstave+dstave : pstave-dstave;
      getPointsOnBarrelPseudoLayer(pseudolayer, pstave, pt1, pt2);
      get2dLineIntersection(point1, point2, pt1, pt2, intsec);    
      int intpstave = getPseudoStave(intsec);
      if (intpstave==stave) {
	gotit=true;
	break; // the interaction point is in the considered stave: it's good
      }
    }
    if (gotit) break;
  }

  if (gotit) {
    for (int i=0; i<2; i++)
      intersection[i]=intsec[i];
  } else {
    cout << "ECALGarlicGeometryHelpers::getLineBarrelPseudoLayerIntersection error: " << endl;
    cout << "could not find pslayer-line intersection" << endl;
    cout << "player = " << pseudolayer << " pt1 = " << point1[0] << " " << point1[1] << " pt2 = " << point2[0] << " " << point2[1] << endl;
    intersection[0] = -99999;
    intersection[1] = -99999;
  }
  return;
}



void ECALGarlicGeometryHelpers::getPointsOnBarrelPseudoLayer(int pslayer, int pseudostave, float* point1, float* point2) {
  // get two 2d points on a barrel pseudolayer

  // first point is the central one
  float radius = getPseudoLayerRadius(pslayer);
  float phi = getPseudoStavePhi(pseudostave);

  point1[0] = radius*cos(phi);
  point1[1] = radius*sin(phi);

  float vectorInPLayer[2] = {-sin(phi), cos(phi)};

  for (int i=0; i<2; i++) point2[i]=point1[i]+vectorInPLayer[i];
  return;
}


void  ECALGarlicGeometryHelpers::get2dLineIntersection(const float* p1, const float* p2, const float* p3, const float* p4, float* intersection) {
  // calculate intersection of 2 lines in 3d
  //  line 1 has points p1, p2
  //  line 2 has points p3, p4

  // common denomiator
  float mat[4];
  mat[0] = p1[0]-p2[0];
  mat[1] = p1[1]-p2[1];
  mat[2] = p3[0]-p4[0];
  mat[3] = p3[1]-p4[1];
  float denom = det(mat);

  float temp[4];
  temp[0] = p1[0];
  temp[1] = p1[1];
  temp[2] = p2[0];
  temp[3] = p2[1];
  float r1 = det(temp);
  
  temp[0] = p3[0];
  temp[1] = p3[1];
  temp[2] = p4[0];
  temp[3] = p4[1];
  float r2 = det(temp);

  mat[0] = r1;
  mat[1] = p1[0]-p2[0];
  mat[2] = r2;
  mat[3] = p3[0]-p4[0];
  intersection[0]=det(mat)/denom;

  mat[1] = p1[1]-p2[1];
  mat[3] = p3[1]-p4[1];
  intersection[1]=det(mat)/denom;

  return;
}


void ECALGarlicGeometryHelpers::get3dLinePlaneIntersection(const float* planePoint, const float* planeNormal, 
							   const float* linePoint, const float* lineDirection, 
							   float* intersection) {

  // point in plane p
  // defn of plane:
  //          (p-planePoint).planeNormal = 0
  // point on line q
  // defn of line:
  //          q = linePoint + d*lineDirection
  // at intersection
  //   (linePoint + d*lineDirection - planePoint).planeNormal = 0
  // d = (planePoint-linePoint).planeNormal / lineDirection.planeNormal;

  float pl[3];
  for (int i=0; i<3; i++) pl[i] = planePoint[i] - linePoint[i];

  float pldotnorm(0);
  for (int i=0; i<3; i++) pldotnorm+=pl[i]*planeNormal[i];

  float dirdotnorm(0);
  for (int i=0; i<3; i++) dirdotnorm+=lineDirection[i]*planeNormal[i];

  float d = pldotnorm/dirdotnorm;

  for (int i=0; i<3; i++) intersection[i] = linePoint[i] + d*lineDirection[i];

  return;
}
