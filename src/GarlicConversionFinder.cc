#include "GarlicConversionFinder.hh"
#include <EVENT/TrackerHit.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TDirectory.h>

using std::cout;
using std::endl;

GarlicConversionFinder::GarlicConversionFinder(float b_field, bool makePlots) {

  _bField=b_field;

  _max_dz=5.;

  _max_dxy=2.5;

  _max_angle=0.05;

  _min_conv_radius=15;

  _min_hit_deltaR=100;
  _suspicious_hit_deltaR=100;
  _suspicious_hit_distance=200;
  _max_hits_inside=2;
  _veto_hits_inside=10;

  _makePlots = makePlots;
  if ( _makePlots ) {
    _origdir = gDirectory;
    _fout = new TFile("conversions.root","recreate");
    _h_dz_dxy  = new TH2F("dz_dxy",   "dz_dxy",   100,0.,25.,100,0.,25.);
    _h_dz_ang  = new TH2F("dz_angle" ,"dz_angle" ,100,0.,25.,100,0.,0.5);
    _h_dxy_ang = new TH2F("dxy_angle","dxy_angle",100,0.,25.,100,0.,0.5);
    _h_z_rad   = new TH2F("z_rad",    "z_rad",    100,0,2000,100,0,2500);

    _h_sel_dz_dxy  = new TH2F("sel_dz_dxy",   "sel_dz_dxy",   100,0.,25.,100,0.,25.);
    _h_sel_dz_ang  = new TH2F("sel_dz_angle" ,"sel_dz_angle" ,100,0.,25.,100,0.,0.5);
    _h_sel_dxy_ang = new TH2F("sel_dxy_angle","sel_dxy_angle",100,0.,25.,100,0.,0.5);
    _h_sel_z_rad   = new TH2F("sel_z_rad",    "sel_z_rad",    100,0,2000,100,0,2500);
    _origdir->cd();
  }

}

GarlicConversionFinder::~GarlicConversionFinder() {
  if ( _makePlots ) {
    cout << "hello from GarlicConversionFinder destructor" << endl;
    _fout->cd();
    _fout->Write();
    _fout->Close();
  }
}


void GarlicConversionFinder::calculateConversions() {

  if ( _inputTracks.size()<2 ) return;

  ExtendedTrack* itrk[2]={0,0};
  for (size_t i1=0; i1<_inputTracks.size()-1; i1++) {
    itrk[0] = _inputTracks[i1];
    for (size_t i2=i1+1; i2<_inputTracks.size(); i2++) {
      itrk[1] = _inputTracks[i2];

      //      cout << "track pair " << i1 << " electron? " << _inputTracks[i1]->getElectronSel() << " , " << i2 << " electron? " << _inputTracks[i2]->getElectronSel() << endl;
      //      cout << "    " << itrk[0] << " " << itrk[1] << endl;

      int numberOfElectrons(0);
      if ( _inputTracks[i1]->getElectronSel()>0 ) numberOfElectrons++;
      if ( _inputTracks[i2]->getElectronSel()>0 ) numberOfElectrons++;

      GarlicConversionInfo convInfo = isConversion(itrk[0], itrk[1]);
      convInfo.electronID[0]=_inputTracks[i1]->getElectronSel();
      convInfo.electronID[1]=_inputTracks[i2]->getElectronSel();

      if (_makePlots) {
        _fout->cd();
        _h_dz_dxy->Fill( convInfo.dz, convInfo.dxy );
        _h_dz_ang ->Fill( convInfo.dz, convInfo.angle );
        _h_dxy_ang->Fill( convInfo.dxy, convInfo.angle );
        _h_z_rad->Fill( fabs(convInfo.position.Z()), convInfo.radius );

        if ( convInfo.selLoose ) {
          _h_sel_dz_dxy->Fill( convInfo.dz, convInfo.dxy );
          _h_sel_dz_ang ->Fill( convInfo.dz, convInfo.angle );
          _h_sel_dxy_ang->Fill( convInfo.dxy, convInfo.angle );
          _h_sel_z_rad->Fill( fabs(convInfo.position.Z()), convInfo.radius );
        }
      }

      if ( (convInfo.selLoose || convInfo.selTight) && numberOfElectrons>0 ) {
	//        cout << "selected conversion with " << numberOfElectrons << " electrons, tight ? " << convInfo.selTight << endl;
        _conversions.push_back( convInfo );
      }
    }
  }

  // check for any doubly-used tracks, resolve ambiguities if necessary

//  cout << "now check for double use of tracks in conversions: " << endl;
//  for (size_t k=0; k<_conversions.size(); k++) {
//    cout << k << " = " << _conversions[k].etrks[0] << " " << _conversions[k].etrks[1] << endl;
//  }

  while (1) {

    std::vector < ExtendedTrack* > multiplyUsedTracks;
    if ( _conversions.size()>1 ) {
      for (size_t k=0; k<_conversions.size()-1; k++) {
        for (size_t k2=k+1; k2<_conversions.size(); k2++) {
          if ( _conversions[k2].etrks[0]==_conversions[k].etrks[0] ||
               _conversions[k2].etrks[0]==_conversions[k].etrks[1] )
            multiplyUsedTracks.push_back( _conversions[k2].etrks[0] );
          if ( _conversions[k2].etrks[1]==_conversions[k].etrks[0] ||
               _conversions[k2].etrks[1]==_conversions[k].etrks[1] )
            multiplyUsedTracks.push_back( _conversions[k2].etrks[1] );
        }
      }
    }


    //cout << "number of multiply used tracks: " << multiplyUsedTracks.size() << endl;

    if ( multiplyUsedTracks.size()==0 ) break;

    // now resolve any ambiguities
    //    cout << multiplyUsedTracks.size() << " multiply used tracks" << endl;

    std::map <float, ExtendedTrack*> orderedTracks;
    // first order them by momentum
    for (size_t i=0; i<multiplyUsedTracks.size(); i++) {
      Track* trk = multiplyUsedTracks[i]->getTrack();
      float pp = fabs( sqrt( 1. + pow( trk->getTanLambda(), 2 ) ) / trk->getOmega() ); // proportional to total momentum
      orderedTracks[pp] = multiplyUsedTracks[i];
    }

    for ( std::map <float, ExtendedTrack*>::reverse_iterator itt=orderedTracks.rbegin(); itt!=orderedTracks.rend(); itt++) {

      std::map <float, int> affectedConvs; affectedConvs.clear();

      for (size_t ll=0; ll<_conversions.size(); ll++) {
        if ( _conversions[ll].etrks[0]==itt->second || _conversions[ll].etrks[1]==itt->second ) {
	  // make up a score to rank the possibilities
	  float score(0);
	  if ( _conversions[ll].selTight )        score+=10000;          // tight gets the highest score
	  if ( _conversions[ll].electronID[0]>0 ) score+=1000;    // as many ID'd electrons as possible
	  if ( _conversions[ll].electronID[1]>0 ) score+=1000;
	  score -= 10*_conversions[ll].nhitsInside[0]; // hits before the conversion position get minus points
	  score -= 10*_conversions[ll].nhitsInside[1];
	  score += _conversions[ll].momentum.Mag()/1000.; // use momentum to break any remaining symmetries
	  affectedConvs[score]=ll;
	}
      }
      int ireject = affectedConvs.begin()->second;
      _conversions.erase( _conversions.begin() + ireject );
    }
  }
}

GarlicConversionInfo GarlicConversionFinder::isConversion(ExtendedTrack* trk1, ExtendedTrack* trk2) {

  GarlicConversionInfo aa = isConversion( trk1->getTrack(), trk2->getTrack() );
  aa.etrks[0]=trk1;
  aa.etrks[1]=trk2;

  return aa;
}

GarlicConversionInfo GarlicConversionFinder::isConversion(Track* trk1, Track* trk2) {

  int rejected=0;

  GarlicConversionInfo info;

  info.trks[0] = trk1;
  info.trks[1] = trk2;


  const TrackState* tks[2]={0,0};

  for (int i=0; i<2; i++) {
    tks[i] = info.trks[i]->getTrackState( TrackState::AtFirstHit );    
    if ( ! tks[i] ) {
      cout << "WARNING, could not get track state @ first hit, getting from IP" << endl;
      tks[i] = info.trks[i]->getTrackState( TrackState::AtIP );    
    }
  }

  info.selLoose=false;
  info.selTight=false;
  info.mass_e_e=0;
  info.mass_pi_pi=0;
  info.mass_p_pi=0;
  info.dxy=-1;
  info.dz=-1;
  info.radius=-1;
  info.angle=-1;
  info.nhitsInside[0]=0;
  info.nhitsInside[1]=0;


  //  cout << "hello from GarlicConversionFinder::isConversion" << endl;

  // first check for opposite charge
  if ( tks[0]->getOmega()*tks[1]->getOmega()>0 ) {
    //    cout << "rejecting conversion: same charge pair" << endl;
    rejected=1;
    return info;
  };


  // cout << "momenta: ";
  // for (int itt=0; itt<2; itt++) 
  //   cout << 2.98e-4 * fabs(_bField / tks[itt]->getOmega() ) * sqrt ( 1. + pow( tks[itt]->getTanLambda(), 2) ) << " ";
  // cout << endl;

  TVector2 circleCenter[2];

  float distBetweenCenters(0);

  if ( rejected==0 ) {
    // look at XY information
    // find centre of X-Y circles
    for (int itt=0; itt<2; itt++) 
      circleCenter[itt]=getCircleCentre(tks[itt]);

    // distance between centres
    distBetweenCenters = (circleCenter[0]-circleCenter[1]).Mod();;

    // look at how closely the circles touch
    // comapre distance between centres with sum of radii
    // should be 0 for perfectly measured conversion
    float touchDist = fabs( distBetweenCenters -
                            ( fabs(1./tks[0]->getOmega()) + fabs(1./tks[1]->getOmega()) ) );

    //    cout << "dists: " << distBetweenCenters << " " << touchDist << endl;

    // if too far, reject
    if ( touchDist > 15. ) {
      //      cout << "rejecting conversion: too far apart " << touchDist << endl;
      rejected=2;
    }
  }

  TVector3 positionAtConv[2];
  TVector3 momentumAtConv[2];

  if ( rejected==0 ) {
    // get the estimated conversion position and for each track, and its direction at that position
    for (int itt=0; itt<2; itt++) {

      // from LCIO & L3 notes

      // track parameters
      float C = tks[itt]->getOmega();
      float phi0 = tks[itt]->getPhi();
      float d0 = tks[itt]->getD0();
      float z0 = tks[itt]->getZ0();
      float tanL = tks[itt]->getTanLambda();

      // track reference point
      float xr = tks[itt]->getReferencePoint()[0];
      float yr = tks[itt]->getReferencePoint()[1];
      float zr = tks[itt]->getReferencePoint()[2];

      // point of closest approach
      float x0 = xr - d0*sin(phi0);
      float y0 = yr + d0*sin(phi0);

      // point on track nearest parallel intersection
      float x = circleCenter[itt].X() + ( circleCenter[1-itt].X() - circleCenter[itt].X() )*fabs(1./C)/distBetweenCenters;
      float y = circleCenter[itt].Y() + ( circleCenter[1-itt].Y() - circleCenter[itt].Y() )*fabs(1./C)/distBetweenCenters;

      float alpha =   - C*(x-x0)*cos(phi0) - C*(y-y0)*sin(phi0);
      float beta  = 1 - C*(x-x0)*sin(phi0) + C*(y-y0)*cos(phi0);
      float s = atan2(-alpha, beta)/C;
      float z = z0 + zr + s*tanL;

      positionAtConv[itt].SetXYZ(x,y,z);
      float totalmomentum = 2.98e-4 * fabs(_bField / C ) * sqrt ( 1. + pow( tanL, 2) );
      float phiAtConvPos = atan2( sin(phi0) - C*(x-x0), cos(phi0) + C*(y-y0) );
      momentumAtConv[itt].SetXYZ( cos(phiAtConvPos), sin(phiAtConvPos), tks[itt]->getTanLambda() );
      momentumAtConv[itt]*=totalmomentum/momentumAtConv[itt].Mag();
    }

    info.dxy = (positionAtConv[0]-positionAtConv[1]).Perp();
  }

  //  TVector3 aveConvPos;
  if ( rejected==0 ) {
    // take simple average of 2 estimates
    //    ( would be more correct to do a kinematic fit. )
    info.position = 0.5*(positionAtConv[0] + positionAtConv[1]);
    info.radius = info.position.Perp();

    // assume from IP to conversion point.
    info.momentum = info.position;
    info.momentum *= 1./info.momentum.Mag();
    info.momentum *= ( momentumAtConv[0].Mag() + momentumAtConv[1].Mag() );

    if ( info.radius < _min_conv_radius ) {
      //      cout << "rejecting conversion: inside beampipe... radius = " << info.radius << endl;
      rejected=3;
    }
  }


  if ( rejected==0 ) {
    // check that z positions at this point are consistent
    info.dz = fabs( positionAtConv[0].Z() - positionAtConv[1].Z() );
    if ( info.dz > 15. ) {
      //      cout << "rejecting conversion:: conversion candidate tracks no not match in z! " << info.dz << " mm " << endl;
      rejected=4;
    }
  }

  bool looksBad = false;

  if ( rejected==0 ) {
    // check angle between 3-vectors
    info.angle = momentumAtConv[0].Angle( momentumAtConv[1] );

    if ( info.angle > 0.1 ) {
      //      cout << "rejecting conversion: angle " << info.angle << endl;
      rejected=5;
    }

  }


  if ( rejected==0 ) {
    // check that first hit of tracks is not too far from this point
    // check that it's radius is not much inside the conv position

    bool suspiciousHits=false;

    for (int itt=0; itt<2; itt++) {
      TrackerHitVec allHits = info.trks[itt]->getTrackerHits();
      for (size_t ih=0; ih<allHits.size(); ih++) {
        float hitradius(0);
        for (int i=0; i<2; i++) {
          hitradius+=pow(allHits[ih]->getPosition()[i], 2);
        }
        hitradius=sqrt(hitradius);
        if (hitradius<info.radius) info.nhitsInside[itt]++;
        else break; // in case the track loops...
      }
      TrackerHit* firstHit = allHits[0];// assume element0 is the "first" hit
      float hitradius(0);
      for (int i=0; i<2; i++) {
        hitradius+=pow(firstHit->getPosition()[i], 2);
      }
      hitradius=sqrt(hitradius);
      float distfromconv(0);
      for (int i=0; i<3; i++) {
        distfromconv+=pow( firstHit->getPosition()[i]-info.position[i], 2);
      }
      distfromconv=sqrt(distfromconv);
      //      cout << "first hit: radius " << hitradius << "conv point radius = " << info.radius << " ; distance from conv pos = " << distfromconv << endl;
      //      cout << "nhits inside conversion = " << info.nhitsInside[itt] << endl;


      if ( hitradius < info.radius - _min_hit_deltaR  ) looksBad=true;

      if ( hitradius < info.radius - _suspicious_hit_deltaR ||
           distfromconv > _suspicious_hit_distance ) suspiciousHits=true;

    }

    if ( suspiciousHits ) {
      //      cout << "WARNING, suspicious hits on conversion tracks..." << endl;
    }
  }

  if ( rejected==0 ) {
    if (looksBad) {
      //      cout << "hits much too early on track" << endl;
      rejected=6;
    }
  }

  if ( rejected==0 ) {
    if ( (info.nhitsInside[0]>_max_hits_inside && info.nhitsInside[1]>_max_hits_inside) ||
         info.nhitsInside[0]>_veto_hits_inside || info.nhitsInside[1]>_veto_hits_inside) {
      //      cout << "rejecting conversion: too many hits before the conversion position!" <<
      //        info.nhitsInside[0] << " " << info.nhitsInside[1] << endl;
      rejected=7;
    }
  }

  if ( rejected==0 ) {
    // calculate invariant mass in a few hypotheses
    const float m_e = 0.000511;
    const float m_pi = 0.140;
    const float m_p = 0.938;

    const float m_rho = 0.775;
    const float w_rho = 0.146;

    const float m_k0 = 0.498;

    const float m_lambda = 1.115;


    TLorentzVector v[2];
    for (int i=0; i<2; i++) v[i].SetVectM( momentumAtConv[i], m_e );
    info.mass_e_e = (v[0]+v[1]).M();

    for (int i=0; i<2; i++) v[i].SetVectM( momentumAtConv[i], m_pi );
    info.mass_pi_pi = (v[0]+v[1]).M();

    int hiMom = momentumAtConv[0].Mag()>momentumAtConv[1].Mag() ? 0 : 1;
    v[0].SetVectM ( momentumAtConv[hiMom], m_p );
    v[1].SetVectM ( momentumAtConv[1-hiMom], m_pi );
    info.mass_p_pi = (v[0]+v[1]).M();

    //    cout << "SELECTED A CONVERSION!! position " << info.position[0] << " " << info.position[1] << " " << info.position[2] <<
    //      " ; mismatch in XY " << info.dxy << " angle = " << info.angle << " hit inside: " << info.nhitsInside[0] << " , " << info.nhitsInside[1] << " : ";
    // cout << "masses: e-e, pi-pi, p-pi " << info.mass_e_e << " " << info.mass_pi_pi << " " << info.mass_p_pi << endl;
  }

  //  info.sel = (rejected==0);

  info.selLoose =
    rejected==0 &&
    info.dxy   < 2.*_max_dxy &&
    info.dz    < 2.*_max_dz &&
    info.angle < 2.*_max_angle;

  info.selTight =
    rejected==0 &&
    info.dxy   < _max_dxy &&
    info.dz    < _max_dz &&
    info.angle < _max_angle;


  return info;

}

TVector2 GarlicConversionFinder::getCircleCentre(const TrackState* trk) {
  // get centre of track helix in XY

  float p_ref_x = trk->getReferencePoint()[0];
  float p_ref_y = trk->getReferencePoint()[1];

  float Omega = trk->getOmega();
  float d0 = trk->getD0();
  float phi0 = trk->getPhi();

  float p_c_x = p_ref_x + ( 1./Omega - d0 ) * sin ( phi0 );
  float p_c_y = p_ref_y - ( 1./Omega - d0 ) * cos ( phi0 );

  return TVector2( p_c_x, p_c_y );
}
