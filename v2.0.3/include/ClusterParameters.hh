#ifndef __ClusterParameters_HH__
#define __ClusterParameters_HH__

#include "TObject.h"
#include "TBuffer.h"

class ClusterParameters : public TObject {

public:
  
  ClusterParameters();
  
  ~ClusterParameters() {}
  
 
  int ID;
  double theta;
  double phi;
  double Etot;
  double Etot_g;
  double Etot_g_noLC;
  double Etot_g_phi;
  double Etot_g_theta;
  double Etot_g_pos;
  double E_GeV;
  double E_GeV_pre;
  double E_GeV_pre_noPhi;
  double E_GeV_noLC;
  double E_GeV_noTheta;
  double E_GeV_noPhi;
  double E_GeV_en;
  double E_GeV_hits;
  double E_GeV_mix;
  double E_GeV_opt;
  double Es1;
  double Es2;
  double Es3;
  double EnPs;
  double En1odd;
  double En1even;
  double En2odd;
  double En2even;
  double NPs;
  double N1odd;
  double N1even;
  double N2odd;
  double N2even;
  double start;
  double end;
  double depth;
  double COGx;
  double COGy;
  double COGz;
  double POSx;
  double POSy;
  double POSz;
  double distToBiggest;
  double smallestDistToBiggest;
  double distFirstCell;
  double distToTrack;
  double smallestDistToTrack;
  double E9C;
  double E4C;
  double E1C;
  double Chi2_long;
  int nHits;
  int nGhostHits;
  int zone;
  int psHits;
  int pl0Hits;
  double hitDensity;
  double enDensity;
  double dirErr;
  double seedDirErr;
  double Eccentricity;
  double Width;
  double photonProb;
  double photonProbMult;
  double bckgrndProbMult;
  double signalProbMult;
  int surroundingLayers;
  int particleID;
  bool isReal;
  double MLP;
  double finalCut;
  double Volume;
  int AssTo;

  ClassDef(ClusterParameters,1);
};

#endif
  
