#include "ECALGarlicExtendedTrack.hh"
#include "ECALGarlicExtendedHit.hh"
#include <algorithm>

void ExtendedTrack::calculateTrk() {

  float phi0 = _parentTrack->getPhi();
  float d0 = _parentTrack->getD0();
  float z0 = _parentTrack->getZ0();
  float omega = _parentTrack->getOmega();
  float tanLambda = _parentTrack->getTanLambda();
  float bfield = ECALGarlicGeometryParameters::Instance().Get_bField();

  _helix->Initialize_Canonical(phi0, d0, z0, omega, tanLambda, bfield);


  calculate();

  _minDistHitEcalEntry = 99999999;
  TrackerHitVec hits = _parentTrack->getTrackerHits();
  for (size_t i=0; i<hits.size(); i++) {
    float dist(0);
    for (int j=0; j<3; j++) 
      dist+=pow( hits[i]->getPosition()[j] - _ecalEntryPoint[j], 2);
    dist = sqrt(dist);
    if (dist<_minDistHitEcalEntry) _minDistHitEcalEntry=dist;
  }

  return;
}

void ExtendedTrack::calculateMC() {

  float position[3];
  float momentum[3];
  for (int i=0; i<3; i++) {
    momentum[i] = _mcPart->getMomentum()[i];
    position[i] = _mcPart->getVertex()[i];
  }
  float charge = _mcPart->getCharge();
  float bfield = ECALGarlicGeometryParameters::Instance().Get_bField();

  _helix->Initialize_VP(position, momentum, charge, bfield);

  calculate();

  return;
}

void ExtendedTrack::calculate() {

  // initialise helix with track parameters (at closeset approach to IP for now)
  //   should update to get params at the ECAL surface when this is easy to do

  //  float phi0 = _parentTrack->getPhi();
  //  float d0 = _parentTrack->getD0();
  //  float z0 = _parentTrack->getZ0();
  //  float omega = _parentTrack->getOmega();
  //  float tanLambda = _parentTrack->getTanLambda();
  //  float bfield = ECALGarlicGeometryParameters::Instance().Get_bField();
  //
  //  _helix->Initialize_Canonical(phi0, d0, z0, omega, tanLambda, bfield);

  float entryPoint[3]={0,0,0};
  float refPoint[3];
  for (int i=0; i<3; i++)
    refPoint[i] = (_helix->getReferencePoint())[i];

  if( fabs(_helix->getMomentum()[0]) < 0.0000000001 && fabs(_helix->getMomentum()[1]) < 0.0000000001 ) {
    cout << "ExtendedTrack::calculate : Strange Track!" << endl;
  }

  float shortest_time=1.0e+10;

  float intersec[3]={0};
  for (int isec=0; isec<ECALGarlicGeometryParameters::Instance().Get_symmetry()+2; isec++) {
    float time=-999;
    switch (isec) {
    case 0:
      time = _helix->getPointInZ( ECALGarlicGeometryParameters::Instance().Get_zOfEndcap(), refPoint, intersec);
      break;
    case 1:
      time = _helix->getPointInZ(-ECALGarlicGeometryParameters::Instance().Get_zOfEndcap(), refPoint, intersec);
      break;
    default:
      int stave = isec-2;
      float dir[2];
      float point[2];
      float r = ECALGarlicGeometryParameters::Instance().Get_rOfBarrel();
      for (int j=0; j<2; j++) {
	dir[j] = ECALGarlicGeometryParameters::Instance().Get_barrelStaveDir()[stave][j];
	point[j] = dir[j]*r;
      }

      // get helix-plane intersection by hand (the one from HelixClass seems strange/buggy...)
      
      float rr = _helix->getRadius();
      float xc = _helix->getXC();
      float yc = _helix->getYC();

      // change vars to put circle origin @ 0,0
      float point_wrt_centre[2];
      point_wrt_centre[0] = point[0] - xc;
      point_wrt_centre[1] = point[1] - yc;

      float vec_along_line[2];
      vec_along_line[0] =  dir[1];
      vec_along_line[1] = -dir[0];

      // two points on line
      float pts[2][2];
      pts[0][0] = point_wrt_centre[0] + vec_along_line[0];
      pts[0][1] = point_wrt_centre[1] + vec_along_line[1];

      pts[1][0] = point_wrt_centre[0] - vec_along_line[0];
      pts[1][1] = point_wrt_centre[1] - vec_along_line[1];

      float dx, dy;
      dx = pts[1][0] - pts[0][0];
      dy = pts[1][1] - pts[0][1];
      float dr = sqrt( pow(dx, 2) + pow(dy, 2) );
      float D = pts[0][0]*pts[1][1] - pts[1][0]*pts[0][1];

      float delta = pow(rr, 2)*pow(dr, 2) - pow(D, 2);

      if (delta<0) {
	//	cout << "no intersection with this plane" << endl;
      } else {

	float intersections_wrt_centre[2][2]={{0},{0}};

	float sign = dy<0 ? -1. : 1.;

	intersections_wrt_centre[0][0] = ( D*dy + sign*dx*sqrt(delta))/pow(dr, 2);
	intersections_wrt_centre[0][1] = (-D*dx + fabs(dy)*sqrt(delta))/pow(dr, 2);

	intersections_wrt_centre[1][0] = ( D*dy - sign*dx*sqrt(delta))/pow(dr, 2);
	intersections_wrt_centre[1][1] = (-D*dx - fabs(dy)*sqrt(delta))/pow(dr, 2);

	// go back to detector frame (circle centre != 0,0)

	float intersections[2][3]={{0},{0}};

	intersections[0][0]=intersections_wrt_centre[0][0] + xc;
	intersections[0][1]=intersections_wrt_centre[0][1] + yc;

	intersections[1][0]=intersections_wrt_centre[1][0] + xc;
	intersections[1][1]=intersections_wrt_centre[1][1] + yc;

	float phi[2]={0};
	phi[0] = atan2( intersections_wrt_centre[0][1], intersections_wrt_centre[0][0] );
	phi[1] = atan2( intersections_wrt_centre[1][1], intersections_wrt_centre[1][0] );

	float ref_pt_wrt_centre[2]={0};
	ref_pt_wrt_centre[0] = refPoint[0] - xc;
	ref_pt_wrt_centre[1] = refPoint[1] - yc;

	float phi_ref = atan2( ref_pt_wrt_centre[1], ref_pt_wrt_centre[0] );

	// find the corresponding phi values

	float dphi[2];
	dphi[0] = phi[0]-phi_ref;
	dphi[1] = phi[1]-phi_ref;

	float charge = _helix->getCharge();
	for (int i=0; i<2; i++) {
	  if (dphi[i]<0 && charge<0) {
	    dphi[i]+=2.*acos(-1);
	  } else if (dphi[i]>0 && charge>0) {
	    dphi[i]-=2.*acos(-1);
	  }
	}
	float times[2]={0};
	for (int i=0; i<2; i++) {
	  times[i] = -charge*dphi[i]*rr/_helix->getPXY();
	  intersections[i][2] = refPoint[2] + times[i]*_helix->getMomentum()[2];
	}
	int itim = times[0] < times[1] ? 0 : 1;
	for (int k=0; k<3; k++) intersec[k] = intersections[itim][k];
	time = times[itim];
      }


      break;
    }

    if (time>0 && time<shortest_time) {
      for (int i=0; i<3; i++) entryPoint[i]=intersec[i];
      shortest_time=time;
    }

  }


  for (int k=0; k<3; k++) _ecalEntryPoint[k] = entryPoint[k];

  _helix->getExtrapolatedMomentum(entryPoint, _ecalEntryDirection);


  return;
}


float ExtendedTrack::getDistanceToPoint(const float* point) {

  PointInCalo p(point);
  int zone = p.getZone();

  float refPt[3]={0};
  float intersection[6]={0};

  //  cout << "getDistanceToPoint " << _helix << endl;
  //  cout << "- " << _helix->getReferencePoint() << endl;

  for (int i=0; i<3; i++) refPt[i] = _helix->getReferencePoint()[i];

  float time(0);
  float radius(0);
  
  switch (zone) {
  case CALHITZONE_BARREL: // find intersection with cylinder at same R as point
    for (int i=0; i<2; i++) radius+=pow(point[i], 2);
    radius=sqrt(radius);
    time = _helix->getPointOnCircle(radius, refPt, intersection);
    break;
  default: // assume endcap/ring - z plane
    time = _helix->getPointInZ(point[2], refPt, intersection);
    break;
  }

  float dist(0);
  for (int i=0; i<3; i++) 
    dist+=pow( point[i] - intersection[i], 2 );
  dist = sqrt(dist);

  //  cout << "point " << point[0] << " " << point[1] << " " << point[2] ;
  //  cout << " intersec " << intersection[0] << " " << intersection[1] << " " << intersection[2] << " : dist " << dist << endl;

  return dist;
}
