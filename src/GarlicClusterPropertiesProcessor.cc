#include "GarlicClusterPropertiesProcessor.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <string>
#include <marlin/Global.h>
#include "lcio.h"
#include "EVENT/MCParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCCollection.h"

#include "GarlicGeometryParameters.hh"
#include "GarlicAlgorithmParameters.hh"
#include "GarlicExtendedHit.hh"
#include "GarlicExtendedTrack.hh"
#include "GarlicExtendedCluster.hh"

#include <UTIL/CellIDDecoder.h>

#include <gear/GEAR.h>
#include <gear/GearParameters.h>
#include <gear/BField.h>
#include <gear/CalorimeterParameters.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/LayerLayout.h>

#include "TMath.h"
#include "TF1.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace lcio ;
using namespace marlin ;
using std::cout;
using std::endl;

GarlicClusterPropertiesProcessor aGarlicClusterPropertiesProcessor ;

GarlicClusterPropertiesProcessor::GarlicClusterPropertiesProcessor() : Processor("garlicCutCalculator") {

  // modify processor description
  _description = "GarlicClusterPropertiesProcessor makes histograms of EM shower properties";
  // designed to be used on samples of monoenergetic photons and electrons/positrons

  registerProcessorParameter("outputFile",
                             "name of file in which to save clustering histograms",
                             _outfilename,
                             std::string("photonClusters.root"));

  registerProcessorParameter("maxCosth",
                             "max costh to consider",
                             _maxCosth,
                             float(0.98));

  registerProcessorParameter("maxDistPull",
                             "maximum distance to pull cluster by one hit",
                             _maxDistPull, float(10.));

  registerProcessorParameter("hitFraction",
                             "cluster parameters for this fraction of hits",
                             _hitFraction,
                             float(0.9));

}


void GarlicClusterPropertiesProcessor::init() {
  cout << "hello from GarlicClusterPropertiesProcessor::init" << endl;
  _geomSetup=false;

  _fout = new TFile(_outfilename.c_str(),"recreate");
  _tree = new TTree("clus","clus");

  _tree->Branch("eventNumber", &_evtN);

  _tree->Branch("mcEnergy", &_mcE);
  _tree->Branch("mcPDG",    &_mcPDG);

  _tree->Branch("mcTheta", &_mcTheta);
  _tree->Branch("mcPhi", &_mcPhi);

  _tree->Branch("bremsElectron", &_bremsElectron);
  _tree->Branch("convertedPhoton", &_convertedPhoton);
  _tree->Branch("nTrack", &_nTrack);
  _tree->Branch("trackMom", &_trackMom);

  _tree->Branch("ecalEn",&_ecalEn);

  _tree->Branch("clusEn",&_clusEn);
  _tree->Branch("clusHits",&_clusHits);

  _tree->Branch("totEcalHits",&_totEcalHits);

  _tree->Branch( "pointAng", &_pointAng);
  _tree->Branch( "eccen", &_eccen);
  _tree->Branch( "width", &_width);
  _tree->Branch( "vol", &_vol);
  _tree->Branch( "start", &_start);
  _tree->Branch( "end", &_end);
  _tree->Branch( "depth", &_depth);
  _tree->Branch( "reldepth", &_reldepth);
  _tree->Branch( "tubeE", _tubeE, "tubeE[5]/F");
  _tree->Branch( "tubeN", _tubeN, "tubeN[5]/F");
  _tree->Branch( "longE", _longE, "longE[5]/F");
  _tree->Branch( "longN", _longN, "longN[5]/F");
  _tree->Branch( "relLongE", _relLongE, "relLongE[5]/F");
  _tree->Branch( "relLongN", _relLongN, "relLongN[5]/F");
  _tree->Branch( "relrelLongE", _relrelLongE, "relrelLongE[5]/F");
  _tree->Branch( "hitEnMean", &_hitEnMean);
  _tree->Branch( "hitEnRMS", &_hitEnRMS);
  _tree->Branch( "hitEnQ1", &_hitEnQ1);
  _tree->Branch( "hitEnQ2", &_hitEnQ2);
  _tree->Branch( "hitEnQ3", &_hitEnQ3);
  _tree->Branch( "mol90a", &_mol90a);
  _tree->Branch( "mol90b", &_mol90b);
  _tree->Branch( "earlymol90a", &_earlymol90a);
  _tree->Branch( "earlymol90b", &_earlymol90b);
  _tree->Branch( "fracDim", _fracDim, "fracDim[3]/F");
  _tree->Branch( "transRMSa", &_transRMSa);
  _tree->Branch( "transRMSb", &_transRMSb);
  _tree->Branch( "fracPLay", &_fracPLay);
  _tree->Branch( "pLayHole", &_pLayHole);
  _tree->Branch( "nLay", &_nLay);
  _tree->Branch( "nPLay", &_nPLay);
  _tree->Branch( "clMass", &_clMass);

  return;
}

void GarlicClusterPropertiesProcessor::processRunHeader( LCRunHeader* run) {
  cout << "hello from GarlicClusterPropertiesProcessor::processRunHeader" << endl;
}

std::vector < TH2F*> GarlicClusterPropertiesProcessor::makeHistos( std::pair < int, float > pdgEn ) {
  std::vector < TH2F*> hh;
  TString lab = "pdg"; lab+=pdgEn.first; // lab+="_e"; lab+=pdgEn.second; 
  lab+="_";
  TH2F* h;
  std::string vn;

  const float maxLay=35;
  const float maxR  =65;
  const int nbins = 400;

  bool firstTime = _vns.size()==0;

  vn="tube0_E"     ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins, -0.1,1.1); hh.push_back(h);

  cout << lab+vn << " " << lab+vn << " " << nbins << " " << 0.5*pdgEn.second << " " << 2*pdgEn.second << " " << nbins << " " <<  -0.1 << " " << 1.1 << endl;

  vn="tube0_N"     ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins, -0.1,1.1); hh.push_back(h);
  vn="eOnP"        ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins, -0.1,3.); hh.push_back(h);
  vn="logPointAng" ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-5,1); hh.push_back(h);
  vn="start"       ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-1,maxLay); hh.push_back(h);
  vn="reldepth"    ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-1,maxLay); hh.push_back(h);
  vn="relrelLongE0"; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="relrelLongE1"; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="relrelLongE2"; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="hitEnDistr"  ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,5.); hh.push_back(h);
  vn="fracDim"     ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="transRmsA"   ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,maxR); hh.push_back(h);
  vn="transRmsB"   ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,maxR); hh.push_back(h);
  vn="eccen"       ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="molA"        ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,maxR); hh.push_back(h);
  vn="molB"        ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,maxR); hh.push_back(h);
  vn="eccenMol"    ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="earlyMolB"   ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,maxR); hh.push_back(h);
  vn="fracPLay"    ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,1.1); hh.push_back(h);
  vn="pLayHoles"   ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-1,maxLay); hh.push_back(h);
  vn="nLay"        ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-1,maxLay); hh.push_back(h);
  vn="nPLay"       ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-1,maxLay); hh.push_back(h);
  vn="clMass"      ; if ( firstTime ) _vns.push_back(vn); h = new TH2F(lab+vn,lab+vn,nbins,0.5*pdgEn.second,2*pdgEn.second, nbins,-0.1,2.5); hh.push_back(h);

  return hh;
}

void GarlicClusterPropertiesProcessor::fillHistos( std::vector < TH2F*> hh ) {

  //  cout << _mol90a << endl;


  hh[0 ]->Fill( _clusEn, _tubeE[0] );
  hh[1 ]->Fill( _clusEn, _tubeN[0] );
  if ( _trackMom>0 ) hh[2]->Fill( _clusEn, _clusEn/_trackMom );
  hh[3 ]->Fill( _clusEn, log10( _pointAng ) );
  hh[4 ]->Fill( _clusEn, _start );
  hh[5 ]->Fill( _clusEn, _reldepth );
  hh[6 ]->Fill( _clusEn, _relrelLongE[0] );
  hh[7 ]->Fill( _clusEn, _relrelLongE[1] );
  hh[8 ]->Fill( _clusEn, _relrelLongE[2] );
  hh[9 ]->Fill( _clusEn, (_hitEnQ3-_hitEnQ1)/_hitEnQ2 );
  hh[10]->Fill( _clusEn, _fracDim[0] );
  hh[11]->Fill( _clusEn, _transRMSa );
  hh[12]->Fill( _clusEn, _transRMSb );
  hh[13]->Fill( _clusEn, (_transRMSb-_transRMSa)/(_transRMSb+_transRMSa) );
  hh[14]->Fill( _clusEn, _mol90a );
  hh[15]->Fill( _clusEn, _mol90b );
  hh[16]->Fill( _clusEn, (_mol90b-_mol90a)/(_mol90a+_mol90b) );
  hh[17]->Fill( _clusEn, _earlymol90b );
  hh[18]->Fill( _clusEn, _fracPLay );
  hh[19]->Fill( _clusEn, _pLayHole );
  hh[20]->Fill( _clusEn, _nLay );
  hh[21]->Fill( _clusEn, _nPLay );
  hh[22]->Fill( _clusEn, _clMass );

  return;
}


void GarlicClusterPropertiesProcessor::processEvent( LCEvent * evt ) {

  initVars();

  if (evt->getEventNumber()%100==0)
    cout << "hello from GarlicClusterPropertiesProcessor::processEvent " << evt->getEventNumber() << endl;

  _evtN = evt->getEventNumber();

  if (!_geomSetup) setupGeom();

  _mcTheta=-999;
  _mcPhi=-999;

  _convertedPhoton=0;
  _bremsElectron=0;

  //  cout << "hello1" << endl;

  try {
    LCCollection* mcColl = evt->getCollection( "MCParticle" );
    int nmc = mcColl->getNumberOfElements();
    if (nmc>0) {
      MCParticle* mc0 = dynamic_cast<MCParticle*>(mcColl->getElementAt(0));

      _mcTheta=acos(mc0->getMomentum()[2]/mc0->getEnergy());
      _mcPhi = atan2( mc0->getMomentum()[1], mc0->getMomentum()[0] );

      _mcE = mc0->getEnergy();
      _mcE = int(_mcE*100.)/100.;

      _mcPDG = mc0->getPDG();

      if ( mc0->getPDG()==22 && mc0->isDecayedInTracker() ) _convertedPhoton=1;

      // look for on-prompt photon (e.g. brems)
      int nmc = mcColl->getNumberOfElements();
      for (int ih=0; ih<nmc; ih++) {
        MCParticle* mcp = dynamic_cast<MCParticle*>(mcColl->getElementAt(ih));
        if ( mcp->getPDG()==22 ) {
          float posR(0);
          float posZ = fabs( mcp->getVertex()[2] );
          for (int i=0; i<2; i++) {
            posR+=pow(mcp->getVertex()[i], 2);
          }
          posR=sqrt(posR);

          if ( posZ>0.1 && posR>0.1 && posZ<2225 && posR<1808 ) { // not at origin, before TPC endplate/cage
            _bremsElectron=1;
            break;
          }
        }
      }
    }
  }
  catch(DataNotAvailableException err) {};


  if ( abs(_mcPDG)==22 && _convertedPhoton>=1 ) return;
  if ( abs(_mcPDG)==11 && _bremsElectron>=1 ) return;
  if ( fabs(cos(_mcTheta)) >= _maxCosth ) return;



  bool crazyEvent=false;

  float singleTrackDirectionAtEcal[3]={0};
  _nTrack=0;
  _trackMom=-999;
  GarlicExtendedTrack *a_ext_track(NULL);

  // look for any tracks
  TVector3 entryPos(0,0,0);
  try {
    LCCollection* trackColl = evt->getCollection( "MarlinTrkTracks" );
    _nTrack=trackColl->getNumberOfElements();
    if ( trackColl->getNumberOfElements() == 1 ) {
      Track *a_track = dynamic_cast<Track*>(trackColl->getElementAt( 0 ));
      a_ext_track = new GarlicExtendedTrack(a_track);
      for (int i=0; i<3; i++) singleTrackDirectionAtEcal[i] = a_ext_track->getEcalEntryDir()[i];
      _trackMom = a_ext_track->getTotalMomentum();
      for (int l=0; l<3; l++) entryPos[l] = a_ext_track->getEcalEntryPos()[l];
    }
  } catch (DataNotAvailableException err) {};

  if ( abs(_mcPDG)==22 && _nTrack>0 ) return;
  if ( abs(_mcPDG)==11 && _nTrack!=1 ) return;


  // get impact position on ecal
  if ( _nTrack==0 ) {
    entryPos.SetMagThetaPhi( 2000, _mcTheta, _mcPhi );
  }

  // total energy
  std::vector <std::string> hitcolls;
  hitcolls.push_back("ECALBarrel");
  hitcolls.push_back("ECALEndcap");
  _totEcalHits=0;
  _ecalEn=0;
  for (size_t i=0; i<hitcolls.size(); i++) {
    try {
      LCCollection* hitColl = evt->getCollection( hitcolls[i] );
      int nhits = hitColl->getNumberOfElements();
      _totEcalHits+=hitColl->getNumberOfElements();
      for (int ih=0; ih<nhits; ih++) {
        CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(hitColl->getElementAt(ih));
        _ecalEn+=hit->getEnergy();
      }
    }
    catch(DataNotAvailableException err) {};
  }

  // first calculate mean position of all ECAL hits
  // within angular window from entry point
  float meanpos[3]={0};
  float toten(0);
  int tothits(0);

  float closestDist=999;

  std::vector < CalorimeterHit* > vetoHits;
  std::map < float, CalorimeterHit* > hitDists;
  while ( closestDist>_maxDistPull) {
    for (int i=0; i<3; i++) meanpos[i]=0;
    toten=0;
    hitDists.clear();
    for (size_t i=0; i<hitcolls.size(); i++) {
      try {
        LCCollection* hitColl = evt->getCollection( hitcolls[i] );
        int nhits = hitColl->getNumberOfElements();

        // get the hit ID decoder
        if (nhits>0 && !GarlicGeometryParameters::Instance().Get_defaultDecoder() ) {
          CellIDDecoder<CalorimeterHit>* dec = new CellIDDecoder<CalorimeterHit> (hitColl);
          GarlicGeometryParameters::Instance().Set_defaultDecoder(dec);
        }

        // loop over the hits
        for (int ih=0; ih<nhits; ih++) {
          CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(hitColl->getElementAt(ih));
          if ( find( vetoHits.begin(), vetoHits.end(), hit) != vetoHits.end() ) continue;

          TVector3 pp(0,0,0);
          for (int l=0; l<3; l++) pp[l] = hit->getPosition()[l];

          if ( pp.Angle( entryPos ) > 5.*3.14159/180. ) {
            vetoHits.push_back(hit);
            continue;
          }

          float ehit = hit->getEnergy();
          toten+=ehit;
          for (int j=0; j<3; j++) {
            meanpos[j]+=ehit*hit->getPosition()[j];
          }
          tothits++;
        }
      }
      catch(DataNotAvailableException err) {};
    }
    for (int j=0; j<3; j++) {
      meanpos[j]/=toten;
    }


    // now order hits according to distance to line between IP and cog
    float modc(0);
    for (int j=0; j<3; j++) {
      modc+=pow(meanpos[j], 2);
    }
    modc=sqrt(modc);
    for (size_t i=0; i<hitcolls.size(); i++) {
      try {
        LCCollection* hitColl = evt->getCollection( hitcolls[i] );
        int nhits = hitColl->getNumberOfElements();
        for (int ih=0; ih<nhits; ih++) {

          CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(hitColl->getElementAt(ih));
          if ( find( vetoHits.begin(), vetoHits.end(), hit) != vetoHits.end() ) continue;

          float dot(0);
          float modh(0);
          for (int j=0; j<3; j++) {
            dot += hit->getPosition()[j]*meanpos[j];
            modh += pow( hit->getPosition()[j], 2 );
          }
          modh=sqrt(modh);
          float temp = dot / (modc*modh);
          if (temp>1) temp=1;
          else if (temp<-1) temp=-1;
          float angle = acos ( temp );
          float dist = modh*sin(angle);
          while ( hitDists.find(dist)!=hitDists.end() ) {
            dist+=0.01;
          }
          hitDists[dist] = hit;
        }
      }
      catch(DataNotAvailableException err) {};
    }

    closestDist = hitDists.begin() -> first;
    if ( closestDist > _maxDistPull ) { // mean pos has been pulled away from main cluster by some outlying hits
      // remove furthest hit
      vetoHits.push_back(hitDists.rbegin()->second);
      if (vetoHits.size()>5) {
        cout << "warning, vetoing " << vetoHits.size() << " hits far from cluster " << vetoHits.size() << ", event # " << _evtN << endl;
        cout << " entry point: " << entryPos.X() << " " << entryPos.Y() << " " << entryPos.Z() << endl;
        cout << " mcdirection: " << _mcTheta << " " << _mcPhi << " ; " << cos(_mcPhi)*sin(_mcTheta) << " " << sin(_mcPhi)*sin(_mcTheta) << " " << cos(_mcTheta) << endl;
        cout << " mean cluster pos " << meanpos[0] << " " << meanpos[1] << " " << meanpos[2] << endl;

        crazyEvent=true;
        break;
      }
    }
  }


  if ( ! crazyEvent ) {

    // now add 90% of closest hits to cluster
    GarlicExtendedCluster* exClus90  = new GarlicExtendedCluster();
    GarlicExtendedCluster* exClus2RM = new GarlicExtendedCluster();

    if ( _nTrack==1 ) {
      exClus90->setReferenceDirection( singleTrackDirectionAtEcal );
      exClus90->setAssociatedTrack( a_ext_track );

      exClus2RM->setReferenceDirection( singleTrackDirectionAtEcal );
      exClus2RM->setAssociatedTrack( a_ext_track );
    }

    int nhitstoadd = int( _hitFraction * hitDists.size() );

    float totalenergy(0);
    for ( std::map < float, CalorimeterHit* >::iterator itt = hitDists.begin(); itt!=hitDists.end(); itt++ ) {
      totalenergy+=itt->second->getEnergy();
    }



    int nadd(0);
    float eadd(0);
    for ( std::map < float, CalorimeterHit* >::iterator itt = hitDists.begin(); itt!=hitDists.end(); itt++ ) {

      bool addTo90 = nadd++ <= nhitstoadd;
      bool addTo2RM = itt->first < 2*GarlicAlgorithmParameters::Instance().GetMoliereRadius();

      if ( !addTo90 && !addTo2RM ) continue; // completely ignore this hit


      if ( itt->first > 200 ) {
        cout << "WARNING, trying to add very far-away hit to central cluster... dist = " << itt->first << " : " << evt->getEventNumber() << " " << cos(_mcTheta) << " " << hitDists.size() <<
          " eventN = " << _evtN << endl;
        cout << "    " << eadd/totalenergy << endl;
        cout << " entry point: " << entryPos.X() << " " << entryPos.Y() << " " << entryPos.Z() << endl;
        cout << " mcdirection: " << _mcTheta << " " << _mcPhi << " ; " << cos(_mcPhi)*sin(_mcTheta) << " " << sin(_mcPhi)*sin(_mcTheta) << " " << cos(_mcTheta) << endl;
        cout << "vetoed: " << vetoHits.size() << " selected " << hitDists.size() << endl;
        cout << " mean cluster pos " << meanpos[0] << " " << meanpos[1] << " " << meanpos[2] << endl;
        cout << "add90, add2RM " << addTo90 << " " << addTo2RM << endl;
        crazyEvent=true;
        break;
      }



      GarlicExtendedHit* exhit = new GarlicExtendedHit( itt->second );

      if ( addTo90 ) exClus90->addHit(exhit);
      if ( addTo2RM ) exClus2RM->addHit(exhit);

    }

    if ( crazyEvent ) {
      for ( std::map < float, CalorimeterHit* >::iterator itt = hitDists.begin(); itt!=hitDists.end(); itt++ ) {
        cout << itt->first << " ";
      }
      cout << endl;
    }


    if (_mcE>0 && nadd>5 && !crazyEvent) {

      std::pair < int, float > pdgEn( abs(_mcPDG), _mcE );

      if ( _varHistos.find(pdgEn)==_varHistos.end() ) {
        _varHistos[pdgEn]=makeHistos(pdgEn);
      }


      for (int ijk=0; ijk<2; ijk++) {

        GarlicExtendedCluster* exClus = ijk==0 ? exClus90 : exClus2RM;

        _clusEn   = exClus->getEnergy();
        _clusHits = exClus->getNhits();

        _pointAng = exClus->getClusterPointAngle();
        _eccen    = exClus->getEccentricity();
        _width    = exClus->getWidth();
        _vol      = exClus->getVolume();
        _start    = exClus->getStart();
        _end      = exClus->getEnd();
        _depth    = exClus->getMeanDepth();
        _reldepth = exClus->getRelMeanDepth();

        _clMass   = exClus->getClusterMass();

        for (int i=0; i<5; i++) {
          _tubeE[i] = exClus->getTubeEn()[i];
          _tubeN[i] = exClus->getTubeN()[i];
          _longE[i] = exClus->getLongEn()[i];
          _longN[i] = exClus->getLongN()[i];
          _relLongE[i] = exClus->getRelLongEn()[i];
          _relLongN[i] = exClus->getRelLongN()[i];
          _relrelLongE[i] = exClus->getRelRelLongEn()[i];
        }

        _hitEnMean = exClus->getHitMeanEn();
        _hitEnRMS  = exClus->getHitRMSEn();
        _hitEnQ1   = exClus->getHitQ1En();
        _hitEnQ2   = exClus->getHitQ2En();
        _hitEnQ3   = exClus->getHitQ3En();

        _mol90a = exClus->get1dMol90().first;
        _mol90b = exClus->get1dMol90().second;

        _earlymol90a = exClus->getEarly1dMol90().first;
        _earlymol90b = exClus->getEarly1dMol90().second;

        _fracPLay = exClus->getFracPseudoLayers();
        _pLayHole = exClus->getBiggestPlayHole();

        _nLay = exClus->getNLayers();
        _nPLay = exClus->getNPseudoLayers();

        for (int i=0; i<3; i++)
          _fracDim[i] = exClus->getFractalDimension()[i];

        _transRMSa = exClus->getTransverseRMS().first;
        _transRMSb = exClus->getTransverseRMS().second;

        fillHistos(_varHistos[pdgEn]);

      }

    }

    if (exClus90) {delete exClus90; exClus90=0;}
    if (exClus2RM) {delete exClus2RM; exClus2RM=0;}

  }

  return;
}

void GarlicClusterPropertiesProcessor::check( LCEvent * evt ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void GarlicClusterPropertiesProcessor::end(){


  // delete histos with too few entries
  for ( std::map < std::pair < int, float >, std::vector < TH2F*> >::iterator itt=_varHistos.begin(); itt!=_varHistos.end(); itt++ ) {
    if ( itt->second[0]->GetEntries()<100 ) {
      for (size_t j=0; j<itt->second.size(); j++) {
	delete itt->second[j];
      }
      _varHistos.erase(itt);
    }
  }

  if (0) {

    const int MAXPDG=2;
    const int MAXV=50;
    const int MAXR=4;
    const int MAXE=50;

    if ( MAXV< int(_vns.size()) ) cout << "ERROR increase MAXV!" << endl;

    float evtFracs[MAXR] = {0.01, 0.05, 0.95, 0.99};

    float range[MAXPDG][MAXV][MAXR][MAXE];
    float energies[MAXPDG][MAXE];
    int npt[MAXPDG];
    TGraph* gr[MAXPDG][MAXV][MAXR];

    for (int i1=0; i1<MAXPDG; i1++) {
      npt[i1]=0;
      for (int i4=0; i4<MAXE; i4++) {
        energies[i1][i4]=0;
        for (int i2=0; i2<MAXV; i2++) {
          for (int i3=0; i3<MAXR; i3++) {
            range[i1][i2][i3][i4]=-999;
            gr[i1][i2][i3]=NULL;
          }
        }
      }
    }



    TString plname="clusterParms.ps";
    TCanvas* cc = new TCanvas();
    cc->Print(plname+"[");

    std::ofstream myfile;
    myfile.open ("garlicCuts.txt");

    // with remaining ones, get ranges
    for (int ipdg=0; ipdg<2; ipdg++) {

      int ien(0);

      for ( std::map < std::pair < int, float >, std::vector < TH2F*> >::iterator itt=_varHistos.begin(); itt!=_varHistos.end(); itt++ ) {

        if      ( ipdg==0 && abs( itt->first.first ) !=22 ) continue;
        else if ( ipdg==1 && abs( itt->first.first ) !=11 ) continue;

        energies[ipdg][ien]=log10( itt->first.second );

        // get the energy cut
        TH2F* h2 = itt->second[0];
        TH1D* hd = h2->ProjectionX();
        hd->Fit("gaus");
        float mu =  hd->GetFunction("gaus")->GetParameter(1);
        float sig = hd->GetFunction("gaus")->GetParameter(2);
        int emin=hd->GetXaxis()->FindBin(mu-3.*sig);
        int emax=hd->GetXaxis()->FindBin(mu+3.*sig);

        cout << "for energy " << itt->first.second << " use these energy cuts: " << mu << " " << sig << " , " << mu-3.*sig << " " << mu+3.*sig << " , " << emin << " " << emax << endl;

        for (size_t j=0; j<itt->second.size(); j++) {

          TH2F* hh = itt->second[j];

          TH1D* h = hh->ProjectionY("_py", emin, emax);

          float tot = h->Integral();

          //float onePercentLow(999);
          //float twoPercentLow(999);
          //float onePercentHi(999);
          //float twoPercentHi(999);

          float evtFracPoints[MAXR]; for (int i=0; i<MAXR; i++) evtFracPoints[i]=999;

          float lower(0);
          for (int i=1; i<=h->GetNbinsX(); i++) {
            lower+=h->GetBinContent(i)/tot;
            for (int jj=0; jj<MAXR; jj++) {
              if ( evtFracPoints[jj]>998 && lower > evtFracs[jj] ) {
                if ( evtFracs[jj]<0.5 ) {
                  evtFracPoints[jj]=h->GetBinLowEdge(i) - h->GetBinWidth(i);
                } else {
                  evtFracPoints[jj]=h->GetBinLowEdge(i) + 2*h->GetBinWidth(i);
                }
              }
            }
          }

          cout << itt->first.first << " " << itt->first.second << " : " << j << " ";
          for (int jj=0; jj<MAXR; jj++) {
            cout << evtFracPoints[jj] << " ";
          }
          cout << endl;

          for (int jj=0; jj<MAXR; jj++) {
            range[ipdg][j][jj][ien] = evtFracPoints[jj];
          }

        }
        ien++;

        npt[ipdg]=ien;

      }


      for (size_t i=0; i<_vns.size(); i++) {
        if ( npt[ipdg]<=0 ) continue;

        float ymin(9999999999);
        float ymax(-9999999999);

        for (int k=0; k<4; k++) {
          gr[ipdg][i][k] = new TGraph( npt[ipdg], energies[ipdg], range[ipdg][i][k] );
          TString grn = "graph_pdg"; grn+=ipdg; grn+="_"; grn+=_vns[i]; grn+="_range"; grn+=k;
          gr[ipdg][i][k]->SetNameTitle(grn, grn);

          gr[ipdg][i][k]->SetLineColor(k+1);
          gr[ipdg][i][k]->SetMarkerColor(k+1);

          gr[ipdg][i][k]->Fit("pol3","q");
          gr[ipdg][i][k]->GetFunction("pol3")->SetLineColor(k+1);
          gr[ipdg][i][k]->GetFunction("pol3")->SetLineStyle(2);




          if ( ipdg==0 ) myfile << "22 ";
          else           myfile << "11 ";
          myfile << _vns[i] << " ";
          myfile << evtFracs[k] << " ";
          myfile << gr[ipdg][i][k]->GetFunction("pol3")->GetParameter(0) << " ";
          myfile << gr[ipdg][i][k]->GetFunction("pol3")->GetParameter(1) << " ";
          myfile << gr[ipdg][i][k]->GetFunction("pol3")->GetParameter(2) << " ";
          myfile << gr[ipdg][i][k]->GetFunction("pol3")->GetParameter(3) << " ";
          myfile << endl;

          for (int ll=0; ll<npt[ipdg]; ll++) {
            ymin = std::min( ymin, range[ipdg][i][k][ll] );
            ymax = std::max( ymax, range[ipdg][i][k][ll] );
          }

        }

        cc->Clear();
        gr[ipdg][i][0]->Draw("apl");

        gr[ipdg][i][0]->GetHistogram()->GetYaxis()->SetRangeUser( ymin - (ymax-ymin)/5. , ymax + (ymax-ymin)/5. );

        for (int k=1; k<4; k++) {
          gr[ipdg][i][k]->Draw("samepl");
        }
        cc->Print(plname);

      }

    }

    cc->Print(plname+"]");

    myfile.close();

    _fout->cd();

    for (int i1=0; i1<MAXPDG; i1++) {
      for (int i2=0; i2<MAXV; i2++) {
        for (int i3=0; i3<MAXR; i3++) {
          if ( gr[i1][i2][i3] ) gr[i1][i2][i3]->Write();
        }
      }
    }

  }

  _fout->Write(0);
  _fout->Close();

  std::cout << "GarlicClusterPropertiesProcessor::end()  " << name()
            << std::endl ;
}

void GarlicClusterPropertiesProcessor::setupGeom() {

  GarlicAlgorithmParameters::Instance().SetMoliereRadius(20.);
  GarlicAlgorithmParameters::Instance().SetEnergyMIPconversion(140.);

  //  float _x_activeThickness = 0.5; // si thickness

  // size of dead zones, ncells, si thickness, (not in GEAR file, specified in steering file...)
  //  GarlicGeometryParameters::Instance().Set_activeThickness (_x_activeThickness);
  GarlicGeometryParameters::Instance().Set_defaultDecoder  (NULL);

  // b field
  GarlicGeometryParameters::Instance().Set_bField( Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z() );

  // Calorimeter geometry from GEAR
  const gear::CalorimeterParameters& pEcalBarrel = Global::GEAR->getEcalBarrelParameters();
  const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
  const gear::LayerLayout& ecalBarrelLayout = pEcalBarrel.getLayerLayout();
  const gear::LayerLayout& ecalEndcapLayout = pEcalEndcap.getLayerLayout();

  GarlicGeometryParameters::Instance().Set_symmetry ( pEcalBarrel.getSymmetryOrder() );




  std::vector< vector < float > > barrelStaveDir;
  if(pEcalBarrel.getSymmetryOrder() > 1){
    float phi0 = pEcalBarrel.getPhi0();
    for(int i=0; i<pEcalBarrel.getSymmetryOrder(); i++){
      float phi  = phi0 + i*2*acos(-1.)/pEcalBarrel.getSymmetryOrder();
      vector < float > staveDir;
      staveDir.push_back( -sin(phi) );
      staveDir.push_back( cos(phi) );
      staveDir.push_back( 0. );
      barrelStaveDir.push_back(staveDir);
    }
  }
  GarlicGeometryParameters::Instance().Set_barrelStaveDir(barrelStaveDir);
  // radius of barrel front face
  GarlicGeometryParameters::Instance().Set_rOfBarrel ( pEcalBarrel.getExtent()[0] );
  GarlicGeometryParameters::Instance().Set_rMaxOfBarrel ( pEcalBarrel.getExtent()[1] + 30.); // take account of extra spacein corners...

  // z of barrel end
  GarlicGeometryParameters::Instance().Set_zOfBarrel ( pEcalBarrel.getExtent()[3] );

  // inner r of endcap
  GarlicGeometryParameters::Instance().Set_rInnerEcalEndcap ( pEcalEndcap.getExtent()[0] );

  // outer r of endcap
  GarlicGeometryParameters::Instance().Set_rOfEndcap ( pEcalEndcap.getExtent()[1] );

  // z of endcap front face
  GarlicGeometryParameters::Instance().Set_zOfEndcap ( pEcalEndcap.getExtent()[2] );

  // number of layers (include the preshower layer)
  GarlicGeometryParameters::Instance().Set_nBarrelEcalLayers( ecalBarrelLayout.getNLayers() );
  GarlicGeometryParameters::Instance().Set_nEndcapEcalLayers( ecalEndcapLayout.getNLayers() );
  GarlicGeometryParameters::Instance().Set_nPseudoLayers( std::max(ecalBarrelLayout.getNLayers(), ecalEndcapLayout.getNLayers()) );

  GarlicGeometryParameters::Instance().Set_absorberRadiationLength(3.5);

  // layer0 is "preshower" layer
  // layer1 is "first" ECAL layer

  // absorber thickness arrays in barrel and endcap

  float absThick[MAX_NUMBER_OF_LAYERS]={0};
  float cellSize[MAX_NUMBER_OF_LAYERS]={0};
  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    absThick[i] = ecalBarrelLayout.getAbsorberThickness(i);
    cellSize[i] = ecalBarrelLayout.getCellSize0(i);
  }
  GarlicGeometryParameters::Instance().Set_absThicknessBarrelLayer(absThick);
  GarlicGeometryParameters::Instance().Set_padSizeEcal(cellSize);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) absThick[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    absThick[i] = ecalEndcapLayout.getAbsorberThickness(i);
  }
  GarlicGeometryParameters::Instance().Set_absThicknessEndcapLayer(absThick);

  // layer positions: this should be approx position of centre of silicon layer

  float positions[MAX_NUMBER_OF_LAYERS]={0};
  for (int i=0; i<ecalBarrelLayout.getNLayers(); i++) {
    positions[i] = ecalBarrelLayout.getDistance(i); // + ecalBarrelLayout.getThickness(i)/2;
  }
  GarlicGeometryParameters::Instance().Set_positionBarrelLayer(positions);

  for (int i=0; i<MAX_NUMBER_OF_LAYERS; i++) positions[i]=0;
  for (int i=0; i<ecalEndcapLayout.getNLayers(); i++) {
    positions[i] = ecalEndcapLayout.getDistance(i); // + ecalEndcapLayout.getThickness(i)/2;
  }

  GarlicGeometryParameters::Instance().Set_positionEndcapLayer(positions);

  _geomSetup=true;

  return;
}


void GarlicClusterPropertiesProcessor::initVars() {

  _evtN = 0;

  _mcTheta = -999;
  _mcPhi   = -999;
  _mcE=0;
  _mcPDG=0;

  _convertedPhoton = 0;
  _bremsElectron = 0;

  _ecalEn = 0;

  _clusEn = 0;
  _clusHits = 0;

  _totEcalHits = 0;

  _clMass=0;

  _pointAng = 0;
  _eccen = 0;
  _width = 0;
  _vol = 0;
  _start = 0;
  _end = 0;
  _depth = 0;
  _reldepth = 0;

  for (int i=0; i<5; i++) {
    _tubeE[i] = 0;
    _tubeN[i] = 0;
    _longE[i] = 0;
    _longN[i] = 0;
    _relLongE[i] = 0;
    _relLongN[i] = 0;
    _relrelLongE[i] = 0;
  }

  _hitEnMean = 0;
  _hitEnRMS = 0;
  _hitEnQ1 = 0;
  _hitEnQ2 = 0;
  _hitEnQ3 = 0;

  _mol90a = 0;
  _mol90b = 0;

  _earlymol90a = 0;
  _earlymol90b = 0;

  for (int i=0; i<3; i++)
    _fracDim[i] = 0;

  _transRMSa = 0;
  _transRMSb = 0;

  _nTrack = 0;
  _trackMom = 0;

  return;
}
