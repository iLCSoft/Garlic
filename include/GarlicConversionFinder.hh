#ifndef _GARLIC_CONV_FIND_
#define _GARLIC_CONV_FIND_

#include "ECALGarlicExtendedTrack.hh"
#include <EVENT/Track.h>
#include <vector>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

struct GarlicConversionInfo {
  ExtendedTrack* etrks[2];
  Track* trks[2];
  bool selTight;
  bool selLoose;
  float mass_e_e;
  float mass_pi_pi;
  float mass_p_pi;
  float dxy;
  float dz;
  float radius;
  float angle;
  int nhitsInside[2];
  bool electronID[2];
  TVector3 position;
  TVector3 momentum;
};


class GarlicConversionFinder {

public:
  GarlicConversionFinder(float b_field, bool makePlots=false);
  ~GarlicConversionFinder();

  void setTracks( std::vector<ExtendedTrack* > & tracks ) {
    _conversions.clear();
    _inputTracks.clear();
    _inputTracks = tracks;
  }

  //  std::vector < std::pair < ExtendedTrack* , ExtendedTrack* > > getConversions() {
  std::vector < GarlicConversionInfo > getConversions() {
    calculateConversions();
    return _conversions;
  }

private:

  std::vector < ExtendedTrack* > _inputTracks;

  //  std::vector < std::pair < ExtendedTrack* , ExtendedTrack* > > _conversions;

  std::vector < GarlicConversionInfo > _conversions;

  TVector2 getCircleCentre(const TrackState* trk);

  GarlicConversionInfo isConversion(ExtendedTrack* trk1, ExtendedTrack* trk2);
  GarlicConversionInfo isConversion(Track* trk1, Track* trk2);
  void calculateConversions();

  float _bField;
  float _max_dz;
  float _max_dxy;
  float _max_angle;
  float _min_conv_radius;

  float _min_hit_deltaR;
  float _suspicious_hit_deltaR;
  float _suspicious_hit_distance;

  int _max_hits_inside;
  int _veto_hits_inside;

  bool _makePlots;

  TDirectory* _origdir;
  TFile* _fout;
  TH2F* _h_dz_dxy;
  TH2F* _h_dxy_ang;
  TH2F* _h_dz_ang;
  TH2F* _h_z_rad;

  TH2F* _h_sel_dz_dxy;
  TH2F* _h_sel_dxy_ang;
  TH2F* _h_sel_dz_ang;
  TH2F* _h_sel_z_rad;


};

#endif
