#ifndef GarlicClusterPropertiesProcessor_h
#define GarlicClusterPropertiesProcessor_h
#include "marlin/Processor.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"

using namespace marlin ;


class GarlicClusterPropertiesProcessor : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new GarlicClusterPropertiesProcessor ; }
  
  GarlicClusterPropertiesProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
    
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

 protected:

  virtual void initVars();

  float _maxCosth;
  float _maxDistPull;
  float _hitFraction;

  TFile* _fout;
  TTree* _tree;

  int _evtN;

  float _mcTheta;
  float _mcPhi;

  float _mcE;
  int   _mcPDG;

  int _convertedPhoton;
  int _bremsElectron;

  float _ecalEn;

  float _clusEn;
  int _clusHits;

  int _totEcalHits;

  float _pointAng;
  float _eccen;
  float _width;
  float _vol;
  float _start;
  float _end;
  float _depth;
  float _reldepth;

  float _tubeE[5];
  float _tubeN[5];
  float _longE[5];
  float _longN[5];
  float _relLongE[5];
  float _relLongN[5];
  float _relrelLongE[5];

  float _hitEnMean;
  float _hitEnRMS;
  float _hitEnQ1;
  float _hitEnQ2;
  float _hitEnQ3;

  float _mol90a;
  float _mol90b;

  float _earlymol90a;
  float _earlymol90b;

  float _fracDim[3];
  
  float _transRMSa;
  float _transRMSb;

  int _nTrack;
  float _trackMom;

  int _nLay;
  int _nPLay;

  float _fracPLay;
  float _pLayHole;

  float _clMass;

  std::string _outfilename;

  void setupGeom();
  bool _geomSetup;

  std::map < std::pair < int, float >, std::vector < TH2F*> > _varHistos;
  std::vector < TH2F*> makeHistos( std::pair < int, float > pdgEn );
  void fillHistos( std::vector < TH2F*> hh );


  std::vector < std::string > _vns;

};


#endif



