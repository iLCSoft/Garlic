#ifndef ECALGARLICEXTENDEDTRACK_HH_
#define ECALGARLICEXTENDEDTRACK_HH_

#include <lcio.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <HelixClass.h>
#include <iostream>
#include "ECALGarlicGeometryParameters.hh"

using std::cout;
using std::endl;

class ExtendedTrack {

public:

  ExtendedTrack(lcio::Track* track) {
    _parentTrack=track;
    _mcPart=NULL;
    _helix = new HelixClass();
    calculateTrk();
  }

  ExtendedTrack(lcio::MCParticle* mcp) {
    _mcPart = mcp;
    _parentTrack=NULL;
    _helix = new HelixClass();
    calculateMC();
  }


  ~ExtendedTrack(){
    if (_helix) {
      delete _helix;
      _helix=NULL;
    }
  }

  const float* getMomentum() {return _helix->getMomentum();}
  float getTotalMomentum() {
    float mom(0);
    for (int i=0; i<3; i++) mom+=pow(_helix->getMomentum()[i], 2);
    mom=sqrt(mom);
    return mom;
  }
  HelixClass* getHelix() {return _helix;}
  lcio::Track* getTrack() {return _parentTrack;}
  const float* getEcalEntryPos() {return _ecalEntryPoint;}
  const float* getEcalEntryDir() {return _ecalEntryDirection;}

  float getDistanceToPoint(const float* point);
  float getMinDist_HitEcalEntry() {return _minDistHitEcalEntry;}
  
private:    

  lcio::Track* _parentTrack;
  lcio::MCParticle* _mcPart;
  
  HelixClass* _helix;

  float _ecalEntryPoint[3];
  float _ecalEntryDirection[3];
  float _minDistHitEcalEntry;

  void calculateTrk();
  void calculateMC();
  void calculate();

};

#endif
