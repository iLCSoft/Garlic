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

  ExtendedTrack(lcio::Track* track, bool verb=false) {
    _ttverbose=verb;
    _parentTrack=track;
    _mcPart=NULL;
    _electronSelection=0;
    _helixECAL = new HelixClass();
    _helixIP = new HelixClass();
    calculateTrk();
  }

  ExtendedTrack(lcio::MCParticle* mcp, bool verb=false) {
    _ttverbose=verb;
    _mcPart = mcp;
    _parentTrack=NULL;
    _helixECAL = new HelixClass();
    _helixIP = new HelixClass();
    calculateMC();
  }


  ~ExtendedTrack(){
    if (_helixECAL) {
      delete _helixECAL;
      _helixECAL=NULL;
    }
    if (_helixIP) {
      delete _helixIP;
      _helixIP=NULL;
    }
  }

  const float* getMomentum(bool atIP=true) {
    return atIP ? _helixIP->getMomentum() : _helixECAL->getMomentum();
  }

  float getTotalMomentum(bool atIP=true) {
    float mom(0);
    for (int i=0; i<3; i++) mom+= atIP ? pow(_helixIP->getMomentum()[i], 2) : pow(_helixECAL->getMomentum()[i], 2);
    mom=sqrt(mom);
    return mom;
  }

  HelixClass* getHelix() {return _helixECAL;}
  HelixClass* getHelixIP() {return _helixIP;}
  lcio::Track* getTrack() {return _parentTrack;}
  const float* getEcalEntryPos() {return _ecalEntryPoint;}
  const float* getEcalEntryDir() {return _ecalEntryDirection;}

  float getDistanceToPoint(const float* point, bool atIP=false);
  float getMinDist_HitEcalEntry() {return _minDistHitEcalEntry;}

  void setElectronSel(int esel) {_electronSelection=esel;}
  int getElectronSel() {return _electronSelection;}

private:    

  bool _ttverbose;

  lcio::Track* _parentTrack;
  lcio::MCParticle* _mcPart;
  
  HelixClass* _helixECAL;
  HelixClass* _helixIP;

  float _ecalEntryPoint[3];
  float _ecalEntryDirection[3];
  float _minDistHitEcalEntry;

  float _helixECALRefPoint[3];
  float _helixIPRefPoint[3];

  void calculateTrk();
  void calculateMC();
  void calculate();

  int _electronSelection;

};

#endif
