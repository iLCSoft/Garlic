#ifndef __MCPhoton_HH__
#define __MCPhoton_HH__

#include "TObject.h"
#include "TBuffer.h"

class MCPhoton : public TObject {

public:
  
  MCPhoton();
  
  ~MCPhoton() {}
  
 
  float cosTheta;
  float phi;
  float Etot;
  float E_GeV;
  float RecEnRatio;
  int nHits;
  int zone;
  bool interaction;
  bool rec;
  float smallestDistToTrack;
  float distToTrack;
  float smallestDistToNextPhoton;
 
  ClassDef(MCPhoton,1);
};

#endif
  
