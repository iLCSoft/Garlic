#include "ECALGarlicExtendedTrack.hh"
#include "ECALGarlicExtendedHit.hh"
#include <algorithm>

void ExtendedTrack::calculateTrk() {
  const TrackState* tst = 0;
  float bfield = ECALGarlicGeometryParameters::Instance().Get_bField();

  //---------------------------
  // track at IP
  //---------------------------
  tst = _parentTrack->getTrackState( TrackState::AtIP );
  float phi0      = tst->getPhi();
  float d0        = tst->getD0();
  float z0        = tst->getZ0();
  float omega     = tst->getOmega();
  float tanLambda = tst->getTanLambda();
  _helixIP->Initialize_Canonical(phi0, d0, z0, omega, tanLambda, bfield);

  for (int i=0; i<3; i++) 
    _helixIPRefPoint[i] = tst->getReferencePoint()[i];

  //---------------------------
  // track @ CALO entrance
  //---------------------------
  tst = _parentTrack->getTrackState( TrackState::AtCalorimeter );
  if ( !tst ) {
    cout << "warning, could not get track state @ ECAL entrance, trying last hit..." << endl;
    tst = _parentTrack->getTrackState( TrackState::AtLastHit );
  }
  if ( !tst ) {
    cout << "warning, could not get track state @ last hit, trying IP..." << endl;
    tst = _parentTrack->getTrackState( TrackState::AtIP );
  }
  if ( !tst ) {
    cout << "ERROR, could not get any track state! giving up..." << endl;
    return;
  }

  // float phi0      = _parentTrack->getPhi();
  // float d0        = _parentTrack->getD0();
  // float z0        = _parentTrack->getZ0();
  // float omega     = _parentTrack->getOmega();
  // float tanLambda = _parentTrack->getTanLambda();

  phi0      = tst->getPhi();
  d0        = tst->getD0();
  z0        = tst->getZ0();
  omega     = tst->getOmega();
  tanLambda = tst->getTanLambda();


  _helixECAL->Initialize_Canonical(phi0, d0, z0, omega, tanLambda, bfield);

  for (int i=0; i<3; i++) 
    _helixECALRefPoint[i] = tst->getReferencePoint()[i];

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

  _helixECAL->Initialize_VP(position, momentum, charge, bfield);

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

  if (_ttverbose) {
    cout << "hello from ExtendedTrack::calculate" << endl;

    cout << "reference point: " << _helixECALRefPoint[0] << " " << _helixECALRefPoint[1] << " " << _helixECALRefPoint[2] << endl;
    cout << "rel ref pt: " << _helixECAL->getReferencePoint()[0] << " " << 
      _helixECAL->getReferencePoint()[1] << " " << _helixECAL->getReferencePoint()[2] << endl;
  }

  float entryPoint[3]={0,0,0};
  float refPoint[3];
  for (int i=0; i<3; i++) {
    refPoint[i] = (_helixECAL->getReferencePoint())[i];  // this is in Helix's frame (0,0) = center of helix
  }

  if( fabs(_helixECAL->getMomentum()[0]) < 0.0000000001 && fabs(_helixECAL->getMomentum()[1]) < 0.0000000001 ) {
    cout << "ExtendedTrack::calculate : Strange Track!" << endl;
  }

  float shortest_time=1.0e+10;

  float intersec[3]={0};
  for (int isec=0; isec<ECALGarlicGeometryParameters::Instance().Get_symmetry()+2; isec++) {
    float time=-999;
    float xyz[3];

    if (_ttverbose) {
      cout << "testing intersection : " << isec << endl;
    }

    switch (isec) {
    case 0:
      //      time = _helix->getPointInZ( ECALGarlicGeometryParameters::Instance().Get_zOfEndcap(), refPoint, intersec);
      xyz[2]=ECALGarlicGeometryParameters::Instance().Get_zOfEndcap() - _helixECALRefPoint[2]; // in Helix frame
      time = _helixECAL->getPointInZ( xyz[2], refPoint, intersec);
      // for (int i=0; i<3; i++) intersec[i]+=_helixECALRefPoint[i];
      break;
    case 1:
      xyz[2]=-ECALGarlicGeometryParameters::Instance().Get_zOfEndcap() - _helixECALRefPoint[2];
      time = _helixECAL->getPointInZ(xyz[2], refPoint, intersec); // in helix frame
      break;
    default:

      int stave = isec-2;
      float dir[2];
      float point[2];
      float r = ECALGarlicGeometryParameters::Instance().Get_rOfBarrel();
      for (int j=0; j<2; j++) {
	dir[j] = ECALGarlicGeometryParameters::Instance().Get_barrelStaveDir()[stave][j];
	//	point[j] = dir[j]*r;
	point[j] = dir[j]*r - _helixECALRefPoint[j]; // translate to helix frame
      }

      if (_ttverbose) {
	cout << "DIR: " << dir[0] << " " << dir[1] << endl;
	cout << "INT: " << dir[0]*r << " " << dir[1]*r << endl;
	cout << "PNT: " << point[0] << " " << point[1] << endl;
      }

      // get helix-plane intersection by hand (the one from HelixClass seems strange/buggy...)
      
      float rr = _helixECAL->getRadius();
      float xc = _helixECAL->getXC();
      float yc = _helixECAL->getYC();

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


	if (_ttverbose) {
	  cout << "no intersection with this plane" << endl;
	}

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

	if (_ttverbose) {
	  cout << "initial dPhi: " << dphi[0] << " " << dphi[1] << endl;
	}

	float charge = _helixECAL->getCharge();

	// this assumes we go only forward: OK if track referenced from IP,
	// but not generally the best thing to do, 
	//    particularly if ref point is beyond the ECAL surface
	// probably want to go to closest one

	// for (int i=0; i<2; i++) {
	//   if (dphi[i]<0 && charge<0) {
	//     dphi[i]+=2.*acos(-1);
	//   } else if (dphi[i]>0 && charge>0) {
	//     dphi[i]-=2.*acos(-1);
	//   }
	// }

	// if (_ttverbose) {
	//   cout << "corrected... dPhi: " << dphi[0] << " " << dphi[1] << endl;
	// }

	float times[2]={0};
	for (int i=0; i<2; i++) {
	  times[i] = -charge*dphi[i]*rr/_helixECAL->getPXY();
	  intersections[i][2] = refPoint[2] + times[i]*_helixECAL->getMomentum()[2];
	}


	if (_ttverbose) {
	  cout << "the 2 times: " << times[0] << " " << times[1] << endl;
	  cout << "the 2 intersections: " << 
	    intersections[0][0]<< " " << intersections[0][1] << " " << intersections[0][2] << " , " <<
	    intersections[1][0]<< " " << intersections[1][1] << " " << intersections[1][2] << endl;
	}								   

	//	int itim = times[0] < times[1] ? 0 : 1;
	int itim = fabs(times[0]) < fabs(times[1]) ? 0 : 1;
	for (int k=0; k<3; k++) {
	  intersec[k] = intersections[itim][k];
	}
	time = times[itim];
      }

      break;
    }

    if (_ttverbose) {
      float rr(0);
      for (int i=0; i<3; i++)
	rr+=pow(intersec[i]+_helixECALRefPoint[i],2);
      rr=sqrt(rr);
    
      cout << intersec[0] << " " << intersec[1] << " " << intersec[2] << " : " << time << endl;
      cout << "  " << intersec[0]+_helixECALRefPoint[0] << " " <<
	intersec[1]+_helixECALRefPoint[1] << " " <<
	intersec[2]+_helixECALRefPoint[2] << " : " << rr << endl;
    }

    //    if (time>0 && time<shortest_time) {
    if (fabs(time)<fabs(shortest_time)) {
      for (int i=0; i<3; i++) entryPoint[i]=intersec[i];
      shortest_time=time;
    }

  }


  for (int k=0; k<3; k++) _ecalEntryPoint[k] = entryPoint[k] + _helixECALRefPoint[k]; // move back to detector frame (wrt IP)

  _helixECAL->getExtrapolatedMomentum(entryPoint, _ecalEntryDirection);

  if (_ttverbose) {
    cout << "entry point: " << _ecalEntryPoint[0] << " " << _ecalEntryPoint[1] << " " << _ecalEntryPoint[2] << endl;
    cout << "ref   point: " << _helixECALRefPoint[0] << " " << _helixECALRefPoint[1] << " " << _helixECALRefPoint[2] << endl;
    cout << "entry dirn : " << _ecalEntryDirection[0] << " " << _ecalEntryDirection[1] << " " << _ecalEntryDirection[2] << endl;
  }

  return;
}


float ExtendedTrack::getDistanceToPoint(const float* point, bool atIP) {

  HelixClass* helix = atIP ? _helixIP : _helixECAL ;
  float * helixRefPoint = atIP ? _helixIPRefPoint : _helixECALRefPoint;

  float point_helixRF[3];
  for (int i=0; i<3; i++)
    point_helixRF[i]=point[i]-helixRefPoint[i];

  PointInCalo p(point);
  int zone = p.getZone();

  // check if the zone is "unknown"
  // usually means not inside the calo volume
  // just use z to decide
  if ( zone == CALHITZONE_UNKNOWN ) { 
    if ( fabs ( p.getPosition()[2] ) < ECALGarlicGeometryParameters::Instance().Get_zOfEndcap() ) // treat as barrel
      zone=CALHITZONE_BARREL;
    else 
      zone=CALHITZONE_ENDCAP; // otherwise endcap
  }


  float refPt[3]={0};
  float intersection[6]={0};

  for (int i=0; i<3; i++) refPt[i] = helix->getReferencePoint()[i];

  float time(0);
  float radius(0);
  
  switch (zone) {
  case CALHITZONE_BARREL: // find intersection with cylinder at same R as point
    //    for (int i=0; i<2; i++) radius+=pow(point[i], 2);
    for (int i=0; i<2; i++) radius+=pow(point_helixRF[i], 2);
    radius=sqrt(radius);
    //    cout << " in barrel...radius = " << radius << endl;
    time = helix->getPointOnCircle(radius, refPt, intersection);
    break;
  default: // assume endcap/ring - z plane
    //    cout << " in endcap! " << endl;
    time = helix->getPointInZ(point_helixRF[2], refPt, intersection);
    break;
  }

  float dist(0);
  for (int i=0; i<3; i++) 
    dist+=pow( point_helixRF[i] - intersection[i], 2 );
  dist = sqrt(dist);

  //  cout << "point " << point[0] << " " << point[1] << " " << point[2] ;
  //  cout << " intersec " << intersection[0] << " " << intersection[1] << " " << intersection[2] << " : dist " << dist << endl;

  return dist;
}
