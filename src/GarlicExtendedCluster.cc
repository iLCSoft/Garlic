#include "GarlicExtendedCluster.hh"
#include "GarlicExtendedHit.hh"
#include "GarlicClusterVarsGenericObject.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <assert.h>

#include <UTIL/LCRelationNavigator.h>
#include "EVENT/SimCalorimeterHit.h"

#include "GarlicAlgorithmParameters.hh"
#include "GarlicGeometryParameters.hh"

//#include "ECALGarlicNeuralNetwork.hh"
//#include "pdf_photon_v5.hh"
//#include "pdf_electron_v5.hh"
//#include "pdf_photon_v4.hh"
//#include "pdf_electron_v4.hh"

using std::vector;
using std::min;
using std::max;

#include "TVector3.h"
#include "TF1.h"
#include "TFile.h"

#include "CLHEP/Vector/ThreeVector.h"

using std::cout;
using std::endl;

int GarlicExtendedCluster::counter=0;


void GarlicExtendedCluster::calculateTrackLikeVars() {

  if (_validTrackLike) return;

  std::map < int, float > playerMIPEnergies;

  // see if a cluster looks like a MIP
  // typically a "track cluster" made around track extrapolation
  int nply = getNPseudoLayers();
  int nhits = getNhits();
  float mipenergy(0);
  if ( getHits()->size()>0 ) {
    for (size_t i=0; i<getHits()->size(); i++) {
      float emip = getHits()->at(i)->getEnergyMIP();
      mipenergy+=emip;

      int iplay=int(getHits()->at(i)->getPseudoLayer());
      if ( playerMIPEnergies.find(iplay)!=playerMIPEnergies.end() ) 
	playerMIPEnergies[iplay]+=emip;
      else
	playerMIPEnergies[iplay] =emip;

    }
  }
  float aveEn = mipenergy/nply;
  // cout << "looks like MIP ? nPseudoLay = " << nply << " nhits = " << nhits << " aveHitEn (MIP) = " << aveEn << endl;
  bool miplike(false);
  if ( nply>5 && float(nhits)/float(nply)<1.5 && nhits < 45 && aveEn < 3.) miplike=true;


  // look for 2 consecutive layers with > 5 MIP energy
  float prevEn(-1);
  float prevLay(-1);
  int iteractionPlay(-1);
  for ( std::map < int, float >::iterator itt=playerMIPEnergies.begin(); itt!=playerMIPEnergies.end(); itt++ ) {
    //cout << itt->first << " " << itt->second << " ; ";
  }
  //  cout << endl;

  for ( std::map < int, float >::iterator itt=playerMIPEnergies.begin(); itt!=playerMIPEnergies.end(); itt++ ) {
    float thisEn = itt->second;
    int thisLay =  itt->first;
    if ( thisLay-prevLay==1 && thisEn>5 && prevEn>5 ) {
      iteractionPlay=itt->first;
      break;
    }
    prevLay=thisLay;
    prevEn=thisEn;
  }

  //  cout << "     interaction p-layer = " << iteractionPlay << endl;

  _isMipLike = miplike;
  _interactionPlayer = iteractionPlay;
  _validTrackLike = true;

  for (int ii=0; ii<3; ii++)
    _interactionPosition[ii]=0;

  if ( iteractionPlay>=0 ) {
    float tote(0);
    for (size_t i=0; i<getHits()->size(); i++) {
      int iplay=int(getHits()->at(i)->getPseudoLayer());
      if ( iplay==iteractionPlay ) {
	float emip = getHits()->at(i)->getEnergyMIP();
	tote+=emip;
	for (int ii=0; ii<3; ii++) {
	  _interactionPosition[ii]+=emip*(getHits()->at(i)->getCaloHit()->getPosition()[ii]);
	}
      }
    }
    for (int ii=0; ii<3; ii++)
      _interactionPosition[ii]/=tote;
  }

  //   cout << "   interaction position = ";
  //   for (int i=0; i<3; i++)
  //     cout << _interactionPosition[i] << " ";
  //   cout << endl;

  return;
}


void GarlicExtendedCluster::calculateClusterMass() {

  if ( _validClusterMass ) return;

  _clusterMass=0;

  float forMom[4]={0};

  for (size_t ih = 0; ih<_hitVec.size(); ih++) {
    CalorimeterHit* hit = _hitVec[ih]->getCaloHit();
    float energy = hit->getEnergy();

    float pos[3];
    float absPos(0);
    for (int j=0; j<3; j++) {
      pos[j] = hit->getPosition()[j];
      absPos+=pow(pos[j],2);
    }
    absPos=sqrt(absPos);

    for (int j=0; j<3; j++) {
      forMom[j]+=energy*pos[j]/absPos;
    }
    forMom[3]+=energy;
  }

  _clusterMass = pow(forMom[3],2);
  for (int j=0; j<3; j++) {
    _clusterMass-=pow(forMom[j],2);
  }
  _clusterMass=sqrt(_clusterMass);

  _validClusterMass=true;
  return;
}

void GarlicExtendedCluster::calculateClusterStartPosition() {
  // position of earliest hit near to shower axis

  if (_validClusterStart) return;

  for (int i=0; i<3; i++) _clusterStartPos[i]=-999;

  // assumption==0: line from cog to IP
  //             1: measured cluster axis
  //             2: reference direction
  const int assump = 0;

  float maxDist = 2.; // *GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  int minPl(999999999);
  std::vector < int > bestHits;

  while ( bestHits.size()==0 && maxDist<10 ) {
    //cout << "trying maxdist = " << maxDist << endl;
    for (size_t ih = 0; ih<_hitVec.size(); ih++) {
      //cout << " ..... hit " << ih << endl;
      CalorimeterHit* hit = _hitVec[ih]->getCaloHit();
      float dist = this->getDistToClusterAxis ( hit->getPosition(), assump );
      if ( dist < maxDist*GarlicGeometryParameters::Instance().Get_padSizeEcal()[0] ) {
	int play = _hitVec[ih]->getPseudoLayer();
	if ( play < minPl ) {
	  bestHits.clear();
	  bestHits.push_back(ih);
	  minPl=play;
	} else if ( play==minPl ) {
	  bestHits.push_back(ih);
	}
      }
    }
    maxDist+=1.;
  }

  if ( bestHits.size()>0 ) {
    for (int i=0; i<3; i++) {
      _clusterStartPos[i] = 0;
    }
    _clusterStartPlayer = minPl;
    float toten(0);
    for (size_t j=0; j<bestHits.size(); j++) {
      CalorimeterHit* hh = _hitVec[bestHits[j]]->getCaloHit();
      toten += hh->getEnergy();
      for (int i=0; i<3; i++) {
	_clusterStartPos[i] += hh->getEnergy()*hh->getPosition()[i];
      }
    }
    for (int i=0; i<3; i++) {
      _clusterStartPos[i] /= toten;
    }
  } else {
    cout << "WARNING from GarlicExtendedCluster::calculateClusterStartPosition : could not find earliest hit!!" << endl;
  }

  _validClusterStart=true;
  return;
}



void GarlicExtendedCluster::calculatePlayerHoles(float min_energy_fraction) {

  if ( _validPlayerHoles ) return;

  std::map <int, GarlicExtendedHit*> players;
  for (size_t i=0; i<getHits()->size(); i++) {
    int play = getHits()->at(i)->getPseudoLayer();
    players[play]=getHits()->at(i);
  }
  // find largest hole
  int lastplay(-1);
  int biggestHole(0);
  int beforeHole(-1);
  for ( std::map <int, GarlicExtendedHit*>::iterator itt=players.begin(); itt!=players.end(); itt++) {
    //if (verbose)
    //  cout << "CHECKIT! " << itt->first << " " << lastplay << " " << itt->first-lastplay << " " << biggestHole << endl;

    int thishole = itt->first - lastplay - 1;

    //    if ( lastplay>=0 && thishole>0 ) cout << "got a hole..." << lastplay << " " << itt->first << endl;

    if ( lastplay>=0 && thishole>0 && thishole > biggestHole ) {

      // get energy before and after hole
      float ebefore(0);
      float eafter(0);
      for ( std::map <int, GarlicExtendedHit*>::iterator jtt=players.begin(); jtt!=players.end(); jtt++) {
	//    if (verbose)
	//      cout << "hit in play " << jtt->first << " " << jtt->second->getCaloHit()->getEnergy() << endl;
	if ( jtt->first <= lastplay ) {
	  ebefore+=jtt->second->getCaloHit()->getEnergy();
	} else {
	  eafter+=jtt->second->getCaloHit()->getEnergy();
	}      
      }
      ebefore/=(ebefore+eafter);
      eafter/=(ebefore+eafter);

      //      cout << "before, after: " << ebefore << " " << eafter << endl;

      if ( ebefore < min_energy_fraction || eafter < min_energy_fraction ) { // hole at end/begin of shower
	// do nothing
      } else { // a significant hole
	biggestHole=thishole;
	beforeHole=lastplay;
      }

    }
    lastplay=itt->first;
  }


  // get energy before and after hole
  float ebefore(0);
  float eafter(0);
  for ( std::map <int, GarlicExtendedHit*>::iterator itt=players.begin(); itt!=players.end(); itt++) {
    //    if (verbose)
    //      cout << "hit in play " << itt->first << " " << itt->second->getCaloHit()->getEnergy() << endl;
    if ( itt->first <= beforeHole ) {
      ebefore+=itt->second->getCaloHit()->getEnergy();
    } else {
      eafter+=itt->second->getCaloHit()->getEnergy();
    }      
  }
  ebefore/=(ebefore+eafter);
  eafter/=(ebefore+eafter);

  _biggestPlayHole = biggestHole;
  _biggestPlayHole_eFracBefore=ebefore;
  _biggestPlayHole_eFracAfter=ebefore;

  _validPlayerHoles=true;

  //  cout << "biggest hole..." << biggestHole << " from lay " << beforeHole << " e before/after = " << ebefore << " " << eafter << endl;

  return;
}

void GarlicExtendedCluster::makeLongitudinalProjectionGraph() {

  if (_validLongProj) return;
  else _validLongProj=true;

  if (_hitVec.size()>0) {
    // group hits into pseudolayers
    std::map <int, std::vector <GarlicExtendedHit*> > playerHits;
    for (size_t ih=0; ih<_hitVec.size(); ih++) {
      GarlicExtendedHit* eh=_hitVec[ih];
      int play = eh->getPseudoLayer();

      if (playerHits.find(play)==playerHits.end()) {
        std::vector <GarlicExtendedHit*> temp;
        temp.push_back(eh);
        playerHits[play]=temp;
      } else {
        playerHits[play].push_back(eh);
      }
    }

    float depths[100]={0};
    float energies[100]={0};
    float fluc[100]={0};
    //    float zero[100]={0};
    int npoints(0);

    // calculate energy and average pos/depth in each player
    float earliestx0(999);
    float latestx0(0);

    for (std::map <int, std::vector <GarlicExtendedHit*> >::iterator tt=playerHits.begin(); tt!=playerHits.end(); tt++) {

      float etot(0);
      float meanPos[3]={0,0,0};
      for (size_t i=0; i<tt->second.size(); i++) {
        GarlicExtendedHit* ehit = (tt->second)[i];
        float en = ehit->getCaloHit()->getEnergy();
        etot+=en;
        for (int j=0; j<3; j++) {
          meanPos[j]+=en*ehit->getCaloHit()->getPosition()[j];
        }
      }
      for (int j=0; j<3; j++)
        meanPos[j]/=etot;

      PointInCalo pt(meanPos);
      float x0depth = pt.getX0depth();
      if (x0depth<earliestx0) earliestx0=x0depth;
      if (x0depth>latestx0) latestx0=x0depth;

      // for fluctuations
      float stoch = 0.17;
      float a = etot/pow(stoch, 2.);
      float b = 1./pow(stoch, 2.);
      float sig = sqrt(a)/b;

      depths[npoints]=x0depth;
      energies[npoints]=etot;
      fluc[npoints]=sig;
      npoints++;

    }

    if (earliestx0>998) {
      cout << "WARNING, could not correctly get first point in long profile: x0 " << earliestx0 << endl;
      earliestx0=0;
    }

    /*
      if (_longProf) {delete _longProf; _longProf=NULL;}
      _longProf = new TGraphErrors(npoints, depths, energies, zero, fluc);
      TF1 f1("emshow","x-[3]<0 ? 0 : [0]*(([1]-1)/[2])*pow(([1]-1)/[2]*(x-[3]), [1]-1)*exp(-(([1]-1)/[2])*(x-[3]))/TMath::Gamma([1])", earliestx0, latestx0);
      f1.SetParameter(0, this->getEnergy());
      f1.SetParError(0, this->getEnergy()/2.);
      f1.SetParameter(1, 5.0);
      f1.SetParError(1, 3.0);
      f1.SetParameter(2, 5.0);
      f1.SetParError(2, 3.0);
      f1.FixParameter(3, earliestx0);
      _longProf->Fit("emshow","q");
      TF1* ff = _longProf->GetFunction("emshow");
      _longFitPars[0] = ff->GetParameter(0);
      _longFitPars[1] = ff->GetParameter(1);
      _longFitPars[2] = ff->GetParameter(2);
      _longFitPars[3] = ff->GetParameter(3);
      _longFitPars[4] = ff->GetProb();
      if ( (2*_longFitPars[0] - this->getEnergy())/this->getEnergy() > 5.) ff->SetLineStyle(2);
    */

  }
  return;
}


void GarlicExtendedCluster::calculateTransverseProjectionShapes() {

  //  cout << "hello from calculateTransverseProjectionShapes! " << _validTransClusterShapes << " : " << this->getHits()->size() << endl;


  if (_validTransClusterShapes) return;
  else _validTransClusterShapes=true;

  calculateClusterShapes();

  _transProjAxisLengths[0]=-999;
  _transProjAxisLengths[1]=-999;

  if (this->getHits()->size()>0) {

    // only early part of ECAL:
    // from start of shower, to 5 X0 layer
    float endX0 = getStart() + 5; // this is in x0
    //    std::pair < TH2F*, TH2I* > earlyProjs = this->GetProjectionHistos(0, endX0);
    MakeProjectionHistos(0, endX0);
    std::pair < TH2F*, TH2I* > earlyProjs (_hEnergy, _hHits );

    if ( earlyProjs.first ) {
      std::vector < std::pair < float, float > > radii = getTransRadii(earlyProjs.first);
      _earlyMinMol95=radii[2].first; _earlyMaxMol95=radii[2].second;
      _earlyMinMol90=radii[3].first; _earlyMaxMol90=radii[3].second;
    }

    //    cout << "early : " << _earlyMinMol95 << " " << _earlyMinMol90 << endl;

    cleanUpHistograms();
    earlyProjs.first=NULL;
    earlyProjs.second=NULL;

    // now in all ECAL depth

    // std::pair < TH2F*, TH2I* > projs = this->GetProjectionHistos();
    MakeProjectionHistos(0, 1000); // use hits in all layers
    std::pair < TH2F*, TH2I* > projs (_hEnergy, _hHits );

    TH2F* hh = projs.first;

    if (!hh) {
      cout << "ERROR could not get proj histos! " << this->getHits()->size() << endl;
    } else {
      std::vector < std::pair < float, float > > radii = getTransRadii(hh);
      if ( radii.size()==6 ) {
        _transProjAxisLengths[0] = radii[0].first; _transProjAxisLengths[1] = radii[0].second;
        _transRmsMin =             radii[1].first; _transRmsMax             = radii[1].second;
        _minMol95 =                radii[2].first; _maxMol95                = radii[2].second;
        _minMol90 =                radii[3].first; _maxMol90                = radii[3].second;
        _minMol80 =                radii[4].first; _maxMol80                = radii[4].second;
        _minMol60 =                radii[5].first; _maxMol60                = radii[5].second;

      } else {
        cout << "ERROR couldnt get tranverse radii...? " << radii.size() << endl;
      }
    }
  }

  //  cout << "bye from calculateTransverseProjectionShapes" << endl;

  return;
}


std::vector < std::pair < float, float > > GarlicExtendedCluster::getTransRadii(TH2F* hh) {
  std::vector < std::pair < float, float > > shapes;
  if (!hh) {
    cout << "ERROR could not get proj histos! " << this->getHits()->size() << endl;
  } else {
    float rmsx = hh->GetRMS(1);
    float rmsy = hh->GetRMS(2);
    float cor = hh->GetCorrelationFactor();

    float den = pow(rmsx,2) - pow(rmsy,2);
    float tan2a = den>0 ? 2.*cor*rmsx*rmsy/den : 0;
    float a  = atan(tan2a)/2.;

    float p1 = pow(rmsx*rmsy, 2)*(1. - pow(cor,2) );
    float p2(p1);

    p1/=pow(rmsy*cos(a), 2) - 2*cor*rmsx*rmsy*sin(a)*cos(a) + pow(rmsx*sin(a), 2);
    p2/=pow(rmsy*sin(a), 2) + 2*cor*rmsx*rmsy*sin(a)*cos(a) + pow(rmsx*cos(a), 2);

    p1=sqrt(p1);
    p2=sqrt(p2);

    shapes.push_back( std::pair < float, float > ( min(p1, p2), max(p1, p2) ) );

    //_transProjAxisLengths[0] = min(p1, p2);
    //   _transProjAxisLengths[1] = max(p1, p2);

    TString name = hh->GetName();

    float xx = max( fabs(hh->GetXaxis()->GetBinCenter(1)), fabs(hh->GetXaxis()->GetBinCenter(hh->GetNbinsX())) );
    float yy = max( fabs(hh->GetYaxis()->GetBinCenter(1)), fabs(hh->GetYaxis()->GetBinCenter(hh->GetNbinsY())) );
    float mm = sqrt( pow(xx, 2) + pow(yy, 2) );

    TH1F h_xp ("projx_"+name, "projx_"+name,2*int(mm),-mm,mm);
    TH1F h_yp ("projy_"+name, "projy_"+name,2*int(mm),-mm,mm);
    TH2F h_xyp("projxy_"+name, "projxy_"+name,2*int(mm),-mm,mm,2*int(mm),-mm,mm);

    for (int ibx=1; ibx<=hh->GetNbinsX(); ibx++) {
      float x = hh->GetXaxis()->GetBinCenter(ibx);
      for (int iby=1; iby<=hh->GetNbinsY(); iby++) {
        float en = hh->GetBinContent(ibx, iby);

	// cout << ibx << " " << iby << " " << en << " : ";

        if (en>0) {
          float y = hh->GetYaxis()->GetBinCenter(iby);
          float x_prime = x*cos(a) + y*sin(a);
          float y_prime = -x*sin(a) + y*cos(a);

          h_xp.Fill(x_prime, en);
          h_yp.Fill(y_prime, en);
          h_xyp.Fill(x_prime, y_prime, en);
        }

      }
    }

    //_transRmsMin = std::min(h_xp.GetRMS(), h_yp.GetRMS());
    //_transRmsMax = std::max(h_xp.GetRMS(), h_yp.GetRMS());
    shapes.push_back( std::pair < float, float > ( std::min(h_xp.GetRMS(), h_yp.GetRMS()), std::max(h_xp.GetRMS(), h_yp.GetRMS()) ) );

    std::vector < float > molFracs;
    molFracs.push_back( 0.95 );
    molFracs.push_back( 0.90 );
    molFracs.push_back( 0.80 );
    molFracs.push_back( 0.60 );
    for (size_t i=0; i<molFracs.size(); i++) {
      float mol1 = get1dmoliere(&h_xp, molFracs[i]);
      float mol2 = get1dmoliere(&h_yp, molFracs[i]);
      shapes.push_back( std::pair < float, float > ( std::min( mol1, mol2 ), std::max( mol1, mol2 ) ) );
    }

    if (p1<0 || p2<0) {
      cout << "WARNING strange transverse profile: ";
      cout << this->getHits()->size() << " rmsx,y = " << rmsx << " " << rmsy <<
        " correl = " << cor << " , axis lengths = " << p1 << " " << p2 << endl;
    }

  }

  return shapes;

}



float GarlicExtendedCluster::get1dmoliere(TH1F* h, float fraction) {
  float total = h->Integral();
  float target = fraction*total;

  float minwidth = 999999;
  for (int i=1; i<h->GetNbinsX(); i++) {
    float integ(0);
    for (int j=i; j<=h->GetNbinsX(); j++) {
      if (h->GetBinContent(j)>0) {
        integ+=h->GetBinContent(j);
        if (integ>target) {
          float width = h->GetXaxis()->GetBinCenter(j) - h->GetXaxis()->GetBinCenter(i);
          if (width<minwidth) minwidth=width;
          break;
        }
      }
    }
  }

  return minwidth;
}


void GarlicExtendedCluster::calculateFractalDimension() {

  if (_validFractal) return;
  _validFractal=true;

  std::vector < std::vector < int > > allCells;
  std::vector <int> MSKIJ;

  int ng[nfrac] = {2, 4, 8};

  for (int j=0; j<nfrac; j++) {

    int n = ng[j];

    allCells.clear();

    for (size_t i=0; i<_hitVec.size(); i++) {

      CalorimeterHit* hit = dynamic_cast <CalorimeterHit*> (_hitVec[i]->getCaloHit());

      int M = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["M"];
      int S = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["S-1"];
      int K = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["K-1"];

      int I = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["I"];
      int J = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["J"];

      I/=n;
      J/=n;

      MSKIJ.clear();
      MSKIJ.push_back(M);
      MSKIJ.push_back(S);
      MSKIJ.push_back(K);
      MSKIJ.push_back(I);
      MSKIJ.push_back(J);

      if (find(allCells.begin(), allCells.end(), MSKIJ)==allCells.end())
        allCells.push_back(MSKIJ);

    }

    float hitRatio = float(_hitVec.size())/allCells.size();
    float fracDim = log(hitRatio)/log(n);

    _fractalDimension[j] = fracDim;
  }

  return;
}

void GarlicExtendedCluster::calculateEMFit() {

  makeLongitudinalProjectionGraph();

  return;
}


void GarlicExtendedCluster::calculateClusterShapes() {

  if (_validClusterShapes) {
    return;
  }

  if (_hitVec.size()>MAXHITS) cout << "WARNING, ignoring hits: nhits =" << _hitVec.size() << " MAXHITS=" << MAXHITS << endl;

  if (_hitVec.size()==0) {
    cout << "WARNING, cluster with only " << _hitVec.size() << " hits!" << endl;
    return;
  }

  float hitE[MAXHITS]={0};
  float hitX[MAXHITS]={0};
  float hitY[MAXHITS]={0};
  float hitZ[MAXHITS]={0};

  int nhit(0);
  _energy=0;
  _nhits=_hitVec.size();
  std::vector <int> players;
  std::vector <int> layers;

  for (size_t i=0; i<_hitVec.size(); i++) {
    int player = _hitVec[i]->getPseudoLayer();
    if ( find( players.begin(), players.end(), player)==players.end() ) players.push_back(player);
    int layer = _hitVec[i]->getLayer();
    if ( find( layers.begin(), layers.end(), layer)==layers.end() ) layers.push_back(layer);


    CalorimeterHit* chit = dynamic_cast <CalorimeterHit*> (_hitVec[i]->getCaloHit());
    _energy+=chit->getEnergy();
    if (nhit<MAXHITS) {
      hitE[nhit]=chit->getEnergy();
      hitX[nhit]=chit->getPosition()[0];
      hitY[nhit]=chit->getPosition()[1];
      hitZ[nhit]=chit->getPosition()[2];
      nhit++;
    }
  }

  _nPLayers=players.size();
  _nLayers=layers.size();

  int pmax(-9999);
  int pmin(9999);
  for (int i=0; i<_nPLayers; i++) {
    if ( players[i]>pmax ) pmax=players[i];
    if ( players[i]<pmin ) pmin=players[i];
  }
  _fracPLayers=float(_nPLayers)/float(pmax-pmin+1);

  //  cout << _nPLayers << " " << pmin << " " << pmax << " " << _fracPLayers << endl;
  //  if ( _fracPLayers>1 ) {
  //    for (int i=0; i<_nPLayers; i++) {
  //      cout << players[i] << " ";
  //    }
  //    cout << endl;
  //  }

  if (_clusterShapes) {delete _clusterShapes; _clusterShapes=NULL;}
  _clusterShapes = new ClusterShapes(nhit, hitE, hitX, hitY, hitZ);

  _fitted_width         = _clusterShapes->getWidth();
  _fitted_eccentricity  = _clusterShapes->getElipsoid_eccentricity();
  _fitted_volume        = _clusterShapes->getElipsoid_vol();

  float dist(0);
  for (int i=0; i<3; i++) {
    float temp = _clusterShapes->getCentreOfGravity()[i];
    _centreOfGravity[i] = temp;
    dist+=pow(temp, 2);
  }
  dist = sqrt(dist);
  _theta = acos( _centreOfGravity[2]/dist);
  _phi   = atan2( _centreOfGravity[1], _centreOfGravity[0] );

  //  getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _projectedPosition);
  // projection onto front ECAL face
  // have to be careful about the direction:
  // if it is not externally specified, it is calculated by this function
  // ...so we can't use the getReferenceDirection function to avoid some strange circular logic
  if ( !_validDirection ) { // assume pointing
    getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _clusterShapes->getCentreOfGravity(), _projectedPosition);
    _clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(), _clusterShapes->getCentreOfGravity());
  } else {
    getFrontFaceProj(_clusterShapes->getCentreOfGravity(), getReferenceDirection(), _projectedPosition);

    if ( fabs(_projectedPosition[0])+fabs(_projectedPosition[1])+fabs(_projectedPosition[2]) < 0.01 ) { // no intersection found when using cluster direction, project to IP
      getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _clusterShapes->getCentreOfGravity(), _projectedPosition);
      _clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(), _clusterShapes->getCentreOfGravity());
    } else { // we found a good intersection using cluster direction, use that
      _clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(),  getReferenceDirection() );

      float dist(0);
      for (int i=0; i<3; i++) dist+=pow( _projectedPosition[i]-_clusterShapes->getCentreOfGravity()[i], 2);
      dist=sqrt(dist);
      if ( dist>250 ) { // seems too far away, use proj to ip
	getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _clusterShapes->getCentreOfGravity(), _projectedPosition);
	_clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(), _clusterShapes->getCentreOfGravity());
      }
    }

    // check it's reasonable....
    // calculate distance beween projected and COG
    float dist(0);
    for (int i=0; i<3; i++) dist+=pow( _projectedPosition[i]-_clusterShapes->getCentreOfGravity()[i], 2);
    dist=sqrt(dist);

    if ( dist>500 ) { // seems too far away! maybe initial direction was somewhat crazy? assume pointing...

      cout << "WARNING, some problems projecting cluster COG onto ECAL front face along requested direction" << endl;
      cout << "    too large projection - cog distance: " << dist << endl;
      cout << "--- COG   " << _clusterShapes->getCentreOfGravity()[0] << " " <<
        _clusterShapes->getCentreOfGravity()[1] << " " << _clusterShapes->getCentreOfGravity()[2] << endl;
      cout << "--refDir  " << getReferenceDirection()[0] << " " << getReferenceDirection()[1] << " " << getReferenceDirection()[2] << endl;
      cout << "----proj  " << _projectedPosition[0] << " " << _projectedPosition[1] << " " << _projectedPosition[2] << endl;
      cout << "    assuming direction pointing to IP" << endl;

      float zero[3]={0,0,0};
      setReferenceDirection(zero);
      getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _clusterShapes->getCentreOfGravity(), _projectedPosition);
      _clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(), _clusterShapes->getCentreOfGravity());
    }

  }

  _validClusterShapes=true;

  return;
}

void GarlicExtendedCluster::calculateTrackDistances() {

  if (_validTrackDistances) return;
  _validTrackDistances=true;

  _trackDist_cog=999;
  _trackDist_proj=999;
  _trackDist_min=999;
  _trackDist_first=999;

  

  _closestTrack=0;

  if (_tracks && _tracks->size()>0) {
    for(size_t t_i=0; t_i<_tracks->size(); t_i++) {

      GarlicExtendedTrack *a_track = dynamic_cast<GarlicExtendedTrack*> (_tracks->at(t_i) );

      float dist1 = a_track->getDistanceToPoint( _centreOfGravity );
      _trackDist_cog=min(_trackDist_cog, dist1);

      //      dist = a_track->getDistanceToPoint( _projectedPosition );
      float dist2 = a_track->getDistanceToPoint( getProjectedPosition() );


      float tt(0);
      for (int i=0; i<3; i++) tt+=pow( _centreOfGravity[i]-getProjectedPosition()[i], 2);
      tt=sqrt(tt);

      if (dist2<_trackDist_proj) {
        _trackDist_proj=dist2;
        _closestTrack=a_track;
      }
    }

    for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
      GarlicExtendedHit*   myHit = _hitVec[ihit];
      CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
      for(size_t t_i=0; t_i<_tracks->size(); t_i++) {
        GarlicExtendedTrack *a_track = dynamic_cast<GarlicExtendedTrack*> (_tracks->at(t_i) );
        float pos[3];
        for (int i=0; i<3; i++) pos[i]=a_hit->getPosition()[i];
        float dist = a_track->getDistanceToPoint(pos);
        _trackDist_min=min(_trackDist_min, dist);
      }
    }
  }

  if (_closestTrack) {
    _trackAng_proj = GarlicGeometryHelpers::angle(_closestTrack->getEcalEntryDir(), this->getClusterShape()->getEigenVecInertia() );
  } else
    _trackAng_proj=-999;

  return;
}



void GarlicExtendedCluster::calculateTubeFracs() {

  if (_validTubeFracs) return;
  _validTubeFracs=true;


  calculateClusterShapes();

  assert (_radii.size() == ntubes);

  //  const float origin[3]={0,0,0}; // tube passing through origin and COG

  for (size_t i=0; i<ntubes; i++) {
    _tube_eFracs[i]=0;
    _tube_nFracs[i]=0;
  }

  float etot(0);
  int ntot(0);

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {

    GarlicExtendedHit*   myHit = _hitVec[ihit];
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());

    float hit_en(a_hit->getEnergy());

    etot+=hit_en;
    ntot++;

    //    float distToLine = GetDistToLine( a_hit->getPosition(), origin, _clusterShapes->getCentreOfGravity() );
    float point2[3]={0};
    if ( _validDirection ) { // if ref direction has been set...otherwise, assume pointing
      for (int i=0; i<3; i++) {
        point2[i]=_clusterShapes->getCentreOfGravity()[i] + getReferenceDirection()[i];
      }
    }

    float distToLine = GetDistToLine( a_hit->getPosition(), point2, _clusterShapes->getCentreOfGravity() );

    int iring(-1);

    // check that Moliere radius looks sensible
    float rmol = GarlicAlgorithmParameters::Instance().GetMoliereRadius();
    if ( rmol<0 || rmol>100 ) cout << "WARNING, looks like crazy Moliere radius..." << rmol << " mm " << endl;


    for (size_t ir=0; ir<_radii.size(); ir++) {
      if (distToLine< _radii[ir]*rmol ) {
        iring=ir;
        break;
      }
    }

    if (iring>=0) {
      _tube_eFracs[iring]+=hit_en;
      _tube_nFracs[iring]++;
    }
  } // hit loop

  if (etot>0 && ntot>0) {
    for (size_t ir=0; ir<_radii.size(); ir++) {
      _tube_eFracs[ir]/=etot;
      _tube_nFracs[ir]/=ntot;
    }
  } else {
    cout << "etot or ntot zero ! " << etot << " " << ntot << endl;
  }

  return;
}

void GarlicExtendedCluster::calculateLongShape() {
  // longitudinal shower shape

  if (_validLongShape) return;
  _validLongShape=true;

  assert (_X0s.size() == nlong);

  for (size_t i=0; i<_X0s.size(); i++) {
    _long_eFracs[i]=0;
    _long_nFracs[i]=0;
    _long_rel_eFracs[i]=0;
    _long_rel_rel_eFracs[i]=0;
    _long_rel_nFracs[i]=0;
  }

  float etot(0);
  int ntot(0);

  _start=999;
  _meandepth=0;
  _relmeandepth=0;
  _end=0;

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    GarlicExtendedHit*   myHit = _hitVec[ihit];
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
    float hit_en(a_hit->getEnergy());
    float thisHitx0 = myHit->getX0depth();
    etot+=hit_en;
    ntot++;
    _start = min( _start, thisHitx0);
    _end   = max( _end, thisHitx0);
    _meandepth += thisHitx0*hit_en;
  }
  if (etot>0 && ntot>0) {
    _meandepth/=etot;
  }

  if ( _meandepth<0 ) {
    cout << "GarlicExtendedCluster::calculateLongShape: strange depth... " << _meandepth << endl;
    for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
      GarlicExtendedHit*   myHit = _hitVec[ihit];
      CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
      float hit_en(a_hit->getEnergy());
      float thisHitx0 = myHit->getX0depth();

      cout << "GarlicExtendedCluster::calculateLongShape: en, x0 = " << hit_en << " " << thisHitx0 <<
        " , position : " << a_hit->getPosition()[0] << " " << a_hit->getPosition()[1] << " " << a_hit->getPosition()[2] << " : ";
      cout << " zone   " << myHit->getZone(true) << endl;

    }
  }



  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    GarlicExtendedHit*   myHit = _hitVec[ihit];
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
    float hit_en(a_hit->getEnergy());
    float thisHitx0 = myHit->getX0depth();
    _relmeandepth += (thisHitx0-_start)*hit_en;
  }
  if (etot>0 && ntot>0) {
    _relmeandepth/=etot;
  }

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    GarlicExtendedHit*   myHit = _hitVec[ihit];
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
    float hit_en(a_hit->getEnergy());
    float thisHitx0 = myHit->getX0depth();

    float relativex0 = thisHitx0 - _start;

    int idepth(-1);
    int idepthrel(-1);
    for (size_t ir=0; ir<_X0s.size(); ir++) {
      if (thisHitx0<_X0s[ir]) {
        idepth=ir;
        break;
      }
    }
    for (size_t ir=0; ir<_X0s.size(); ir++) {
      if (relativex0<_X0s[ir]) {
        idepthrel=ir;
        break;
      }
    }
    if (idepth>=0) {
      _long_eFracs[idepth]+=hit_en;
      _long_nFracs[idepth]++;
    }
    if (idepthrel>=0) {
      _long_rel_eFracs[idepthrel]+=hit_en;
      _long_rel_nFracs[idepthrel]++;
    }


    // relative to mean depth from shower start
    idepthrel=-1;
    if (_relmeandepth==0) {
      idepthrel=1;
    } else if ( relativex0 < 0.5*_relmeandepth ) {
      idepthrel=0;
    } else if ( relativex0 < 1.0*_relmeandepth ) {
      idepthrel=1;
    } else if ( relativex0 < 1.5*_relmeandepth ) {
      idepthrel=2;
    } else if ( relativex0 < 2.*_relmeandepth ) {
      idepthrel=3;
    } else {
      idepthrel=4;
    }

    if (idepthrel>=0) {
      _long_rel_rel_eFracs[idepthrel]+=hit_en;
    }


  } // hit loop


  if (etot>0 && ntot>0) {

    for (size_t ir=0; ir<_X0s.size(); ir++) {
      _long_eFracs[ir]/=etot;
      _long_nFracs[ir]/=ntot;

      _long_rel_rel_eFracs[ir]/=etot;

      _long_rel_eFracs[ir]/=etot;
      _long_rel_nFracs[ir]/=ntot;
    }

    //cout << _relmeandepth << " : ";

    //cout << _long_rel_rel_eFracs[0]+_long_rel_rel_eFracs[1] << " , ";
    //cout << _long_rel_rel_eFracs[2]+_long_rel_rel_eFracs[3]+_long_rel_rel_eFracs[4] << endl;

    //cout << endl;

  } else {
    cout << "GarlicExtendedCluster::calculateLongShape: ERROR etot or ntot zero ! " << etot << " " << ntot << endl;
  }

  return;
}


void GarlicExtendedCluster::calculateZone() {

  if (_validZone) return;
  _validZone=true;

  bool hasBarrelHits = false;
  bool hasEndcapHits = false;

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (_hitVec[ihit]->getCaloHit());
    int module = (*GarlicGeometryParameters::Instance().Get_defaultDecoder()) (a_hit)["M"];
    if(0 < module && module < 6) hasBarrelHits=true;
    else                         hasEndcapHits=true;
  }

  if      ( hasBarrelHits && !hasEndcapHits) _zone=CLUS2_LOCATION_BARREL;
  else if (!hasBarrelHits &&  hasEndcapHits) _zone=CLUS2_LOCATION_ENDCAP;
  else if ( hasBarrelHits &&  hasEndcapHits) _zone=CLUS2_LOCATION_OVERLAP;
  else _zone=CLUS2_LOCATION_UNKNOWN;

  return;
}


void GarlicExtendedCluster::calculateHitEnMeasures() {

  if (_validHitEnMeasures) return;
  _validHitEnMeasures=true;

  //  cout << "recalculting hit en measures..." << this->getHits()->size() << endl;


  std::vector <float> hitenergies;
  float etot(0);
  float ensqsum(0);

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    if (! _hitVec[ihit]) cout << "ERROR GarlicExtendedCluster::calculateHitEnMeasures could not get hit " << ihit << " " << _hitVec[ihit] << endl;
    else {
      CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (_hitVec[ihit]->getCaloHit());
      if (!a_hit) cout << "ERROR GarlicExtendedCluster::calculateHitEnMeasures could not find calohit " << ihit << " " << _hitVec[ihit] << " " << _hitVec[ihit]->getCaloHit() << endl;
      else {
        float en = a_hit->getEnergy();
        etot+=en;
        hitenergies.push_back(en);
        ensqsum+=pow(en, 2.);
      }
    }
  }

  sort(hitenergies.begin(), hitenergies.end());

  if (_hitVec.size()>0) {
    _hitMeanEn = etot/_hitVec.size();
    _hitRMSEn = sqrt( ensqsum/_hitVec.size() - pow(_hitMeanEn, 2) );

    int qq = int( _hitVec.size()/4. );
    int iq=qq;
    if (iq<0) iq=0;
    if (iq>=(int) _hitVec.size()) iq=_hitVec.size()-1;
    _hitQ1En = hitenergies[iq];

    iq=int(_hitVec.size()/2);
    if (iq<0) iq=0;
    if (iq>= (int) _hitVec.size()) iq=_hitVec.size()-1;
    _hitQ2En = hitenergies[iq];

    iq=_hitVec.size()-qq;
    if (iq<0) iq=0;
    if (iq>= (int) _hitVec.size()) iq=_hitVec.size()-1;
    _hitQ3En = hitenergies[iq];
  }

  return;
}

void GarlicExtendedCluster::calculateNNResults() {

  //  if (_validNNResults) return;
  //  _validNNResults=true;

//   std::pair <float, bool> nnres = ECALGarlicNeuralNetwork::Instance().getNNoutput(this);
// 
//   _NNval=nnres.first;
//   _NNsel=nnres.second;
// 
//   //  cout << "calculateNNResults: " << _NNval << " " << _NNsel << endl;

  return;
}

void GarlicExtendedCluster::calculatePhotonCutSelResults() {

  if (_validPhotonCutSelResults) return;
  _validPhotonCutSelResults=true;

  //  bool verbosesel=false;
  GarlicClusterSelector* psel = GarlicAlgorithmParameters::Instance().GetClusterSelector(); 
  // (GarlicAlgorithmParameters::Instance().GetPhotonCutFile(), verbosesel);

  psel->setPdg(22); // for photon cuts

  std::map < std::string, float > phoSel = psel->cluster_select(this);

  //  _photonCutSel= psel->cluster_select(this);
  _photonCutSel  = phoSel["TOTAL"];

  _photonCutSelT = phoSel["TRANS"];
  _photonCutSelL = phoSel["LONGT"];
  _photonCutSelP = phoSel["POINT"];
  _photonCutSelM = phoSel["MISCL"];

  //  cout << "GarlicExtendedCluster::calculatePhotonCutSelResults: photon selection results: " <<
  //    "TOTAL:" << _photonCutSel << " T:" << _photonCutSelT << " L:" << _photonCutSelL << " P:" << _photonCutSelP << " M:" << _photonCutSelM << endl;


  return;
}

void GarlicExtendedCluster::calculateElectronCutSelResults() {

  // for the moment, same as photon selection.
  // needs to be updated

  if (_validElectronCutSelResults) return;
  _validElectronCutSelResults=true;

  //  cout << "recalculating calculateElectronCutSelResults..." << this->getHits()->size() << endl;

  GarlicClusterSelector* psel = GarlicAlgorithmParameters::Instance().GetClusterSelector(); // (GarlicAlgorithmParameters::Instance().GetPhotonCutFile(), verbosesel);

  psel->setPdg(11); // for electron cuts

  std::map < std::string, float > eleSel = psel->cluster_select(this);
  //  _electronCutSel = psel->cluster_select(this); // this does not include E/P

  _electronCutSel = eleSel["TOTAL"];

  _electronCutSelT  = eleSel["TRANS"];
  _electronCutSelL  = eleSel["LONGT"];
  _electronCutSelP  = eleSel["POINT"];
  _electronCutSelM  = eleSel["MISCL"];
  _electronCutSelEP = eleSel["E_ONP"];

  //cout << "GarlicExtendedCluster::calculateElectronCutSelResults: tot " << _electronCutSel << " T " << _electronCutSelT << " L " << _electronCutSelL << 
  //    " P " << _electronCutSelP << " M " << _electronCutSelM << " EP " << _electronCutSelEP << endl;

  return;
}


void GarlicExtendedCluster::freeHits() {
  for (size_t i=0; i<_hitVec.size(); i++)
    _hitVec[i]->setCluster(NULL);
  _hitVec.clear();
  changed();
  return;
}



void GarlicExtendedCluster::GetDistancesToCluster(GarlicExtendedCluster *b, double *distances)
{
  double smallest_dist=9999;
  double firstHitDistance=9999;
  int firstBPSLayer=99;
  int firstAPSLayer=99;

  for(size_t b_hit_i=0; b_hit_i<b->getHits()->size(); b_hit_i++) {
    GarlicExtendedHit *myBHit = b->getHits()->at(b_hit_i);

    for(size_t a_hit_i=0; a_hit_i<_hitVec.size(); a_hit_i++) {
      GarlicExtendedHit *myAHit = _hitVec[a_hit_i];

      double dist=Get3dDistance(myAHit->getCaloHit()->getPosition(), myBHit->getCaloHit()->getPosition());

      if(dist<smallest_dist) smallest_dist=dist;

      if(myBHit->getPseudoLayer()<=firstBPSLayer) {
        firstBPSLayer=myBHit->getPseudoLayer();
        if(myAHit->getPseudoLayer()<=firstAPSLayer) {
          firstAPSLayer=myAHit->getPseudoLayer();
          double a_firstDistance = Get2dProjDistance(myAHit->getCaloHit()->getPosition(), myBHit->getCaloHit()->getPosition());
          if(a_firstDistance<firstHitDistance) {
            firstHitDistance=a_firstDistance;
          }
        }
      }
    }
  }
  distances[0] = smallest_dist;
  distances[1] = Get3dDistance(this->getCentreOfGravity(),b->getCentreOfGravity());
  distances[2] = firstHitDistance;
}

// float GarlicExtendedCluster::getDistToClusterLine(GarlicExtendedCluster* cl) {
//   // line from cog to IP
//   float origin[3]={0};
//   return GetDistToLine(this->getCentreOfGravity(), origin, cl->getCentreOfGravity());
// }

//float GarlicExtendedCluster::getDistToClusterAxis ( const float* pos, bool assumeIpPointing ) {
float GarlicExtendedCluster::getDistToClusterAxis ( const float* pos, int assumption ) {
  // assumption==0: line from cog to IP
  //             1: measured cluster axis
  //             2: reference direction
  if ( assumption<0 || assumption>2 ) cout << "ERROR unknown direction assumption asked for  GarlicExtendedCluster::getDistToClusterAxis " << assumption << endl;

  float origin[3]={0};

  if ( assumption==1 ) {
    for (int i=0; i<3; i++) {
      origin[i] = this->getCentreOfGravity()[i] - this->getClusterShape()->getEigenVecInertia()[i];
    }
  } else if ( assumption==2 ) {
    if ( _validDirection ) {
      for (int i=0; i<3; i++) {
        origin[i] = this->getCentreOfGravity()[i] - getReferenceDirection()[i];
      }
    }
  }
  return GetDistToLine( pos, origin, this->getCentreOfGravity() );
}


//void GarlicExtendedCluster::MakeProjectionHistos(float encut, int maxPseudolayerCut) {
void GarlicExtendedCluster::MakeProjectionHistos(float encut, float maxX0) {
  // algo changed a little by daniel
  // take plane normal to cluster COG at front face

  calculateClusterShapes(); // to get COG

  // clean up pre-existing histograms
  cleanUpHistograms();

  const float zero[3]={0,0,0};
  const float zaxis[3]={0,0,1};

  //cout << "hello from MakeProjectionHistos " << encut << " " << maxX0 << endl;
  //cout << " ref direction: " << getReferenceDirection()[0] << " " << getReferenceDirection()[1] << " " << getReferenceDirection()[2] << endl;

  // define the axes in the plane normal to cluster cog projection to front face
  // _xprimeaxis[0] = _origin[1]*zaxis[2] - _origin[2]*zaxis[1];
  // _xprimeaxis[1] = _origin[2]*zaxis[0] - _origin[0]*zaxis[2];
  // _xprimeaxis[2] = _origin[0]*zaxis[1] - _origin[1]*zaxis[0];
  _xprimeaxis[0] = getReferenceDirection()[1]*zaxis[2] - getReferenceDirection()[2]*zaxis[1];
  _xprimeaxis[1] = getReferenceDirection()[2]*zaxis[0] - getReferenceDirection()[0]*zaxis[2];
  _xprimeaxis[2] = getReferenceDirection()[0]*zaxis[1] - getReferenceDirection()[1]*zaxis[0];
  float mag(0);
  for (int i=0; i<3; i++) mag+=pow(_xprimeaxis[i], 2);
  mag=sqrt(mag);
  for (int i=0; i<3; i++) _xprimeaxis[i]/=mag;

  //  _yprimeaxis[0] = _origin[1]*_xprimeaxis[2] - _origin[2]*_xprimeaxis[1];
  //  _yprimeaxis[1] = _origin[2]*_xprimeaxis[0] - _origin[0]*_xprimeaxis[2];
  //  _yprimeaxis[2] = _origin[0]*_xprimeaxis[1] - _origin[1]*_xprimeaxis[0];
  _yprimeaxis[0] = getReferenceDirection()[1]*_xprimeaxis[2] - getReferenceDirection()[2]*_xprimeaxis[1];
  _yprimeaxis[1] = getReferenceDirection()[2]*_xprimeaxis[0] - getReferenceDirection()[0]*_xprimeaxis[2];
  _yprimeaxis[2] = getReferenceDirection()[0]*_xprimeaxis[1] - getReferenceDirection()[1]*_xprimeaxis[0];
  mag=0;
  for (int i=0; i<3; i++) mag+=pow(_yprimeaxis[i], 2);
  mag=sqrt(mag);
  for (int i=0; i<3; i++) _yprimeaxis[i]/=mag;

  // get the seeding hits
  vector < vector < float > > hitLocalposEn;

  float min_xprime( 9999999);
  float max_xprime(-9999999);
  float min_yprime( 9999999);
  float max_yprime(-9999999);

  int nforSeed(0);

  for(size_t hit_i=0;hit_i<_hitVec.size();hit_i++) {
    GarlicExtendedHit *a_ext_hit=_hitVec[hit_i];

    float enmip = a_ext_hit->getEnergyMIP();

    if ( enmip<encut ) continue;
    //    if ( a_ext_hit->getPseudoLayer()<=maxPseudolayerCut ) {
    if ( a_ext_hit->getX0depth()<=maxX0 ) {

      nforSeed++;
      // project onto plane
      float projection[3] = {0};
      //      GarlicGeometryHelpers::get3dLinePlaneIntersection(_origin, _origin, zero, a_ext_hit->getCaloHit()->getPosition(), projection);
      GarlicGeometryHelpers::get3dLinePlaneIntersection(getReferencePoint(), getReferenceDirection(), zero, a_ext_hit->getCaloHit()->getPosition(), projection);

      // now translate to local coordinates
      // y' axis is in "theta" direction
      // x' perp to plane normal and y'
      float xprime(0);
      float yprime(0);
      for (int i=0; i<3; i++) {
        xprime += projection[i]*_xprimeaxis[i];
        yprime += projection[i]*_yprimeaxis[i];
      }

      if (xprime>max_xprime) max_xprime=xprime;
      if (xprime<min_xprime) min_xprime=xprime;

      if (yprime>max_yprime) max_yprime=yprime;
      if (yprime<min_yprime) min_yprime=yprime;

      vector <float> hitinfo;
      hitinfo.push_back(xprime);
      hitinfo.push_back(yprime);
      hitinfo.push_back(enmip);

      hitLocalposEn.push_back(hitinfo);
    }
  }

  if (nforSeed>0) {

    //  now define the histograms
    float cellsize = GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

    // make bins about 1mm size
    //  int nbinsx = int ( (max_xprime-min_xprime)/cellsize ); // this was for one bin per cell
    //  int nbinsy = int ( (max_yprime-min_yprime)/cellsize );
    //int nbinsx = int ( (max_xprime-min_xprime) ); // this is for one bin per mm
    //int nbinsy = int ( (max_yprime-min_yprime) );
    int nbinsx = int ( (max_xprime-min_xprime)/(cellsize/2) ); // bin size half cell size -> Daniel changed, to avoid too many bins 2016/04/01
    int nbinsy = int ( (max_yprime-min_yprime)/(cellsize/2) );

    // make sure this is divisible by int(cellsize)
    // for possible later rebinning to "cellsize" per bin (for seeding)
    int residx = nbinsx%int(cellsize);
    int residy = nbinsy%int(cellsize);
    if (residx>0) residx=int(cellsize)-residx;
    if (residy>0) residy=int(cellsize)-residy;

    nbinsx+=residx;
    max_xprime+=residx;

    nbinsy+=residy;
    max_yprime+=residy;

    nbinsx+=2*int(cellsize);
    nbinsy+=2*int(cellsize);

    cleanUpHistograms();

    if ( nbinsx>1500 || nbinsy>1500 ) {
      cout << "WARNING TOO MANY BINS! cluster too big " << nbinsx << " " << nbinsy << endl;
      assert(0);
    } else if ( nbinsx<1 || nbinsy<1 ) {
      cout << "WARNING TOO FEW BINS! " << nbinsx << " " << nbinsy << endl;
      assert(0);
    } else {


      TString hname = "enHist"; hname+=_thisCounter;
      _hEnergy = new TH2F(hname, hname,
                          nbinsx, min_xprime-cellsize, max_xprime+cellsize,
                          nbinsy, min_yprime-cellsize, max_yprime+cellsize);

      hname = "nHist"; hname+=_thisCounter;
      _hHits = new TH2I(hname, hname,
                        nbinsx, min_xprime-cellsize, max_xprime+cellsize,
                        nbinsy, min_yprime-cellsize, max_yprime+cellsize);


      // and fill this histo
      for (size_t i=0; i<hitLocalposEn.size(); i++) {
        _hEnergy->Fill( hitLocalposEn[i][0], hitLocalposEn[i][1], hitLocalposEn[i][2]);
        _hHits->Fill( hitLocalposEn[i][0], hitLocalposEn[i][1]);
      }
    }
  }
  return;
}

void GarlicExtendedCluster::GetGlobalPositionFromLocal(float localx, float localy, float* pos) {
  //  for (int i=0; i<3; i++) pos[i] = _origin[i] + localx*_xprimeaxis[i] + localy*_yprimeaxis[i];
  for (int i=0; i<3; i++) pos[i] = getReferencePoint()[i] + localx*_xprimeaxis[i] + localy*_yprimeaxis[i];
  return;
}

IMPL::ClusterImpl GarlicExtendedCluster::getClusterImpl() {

  ClusterImpl cluster;
  for (size_t ih=0; ih<_hitVec.size(); ih++) {
    CalorimeterHit* calohit = _hitVec[ih]->getCaloHit();
    cluster.addHit(calohit, 1.);
  }

  cluster.setEnergy( getEnergy() );
  cluster.setPosition( getCentreOfGravity() );

  return cluster;
}


void GarlicExtendedCluster::calculateMCnature() {

  if (_validMC) return;
  else _validMC=true;

  if (_hitVec.size()==0) return;
  if (!_caloHitSimHitRelation) return;
  LCRelationNavigator nav(_caloHitSimHitRelation);

  std::map < int, float > pdgContribsGen;
  std::map < int, float > pdgContribsCal;

  for(size_t hit_i=0;hit_i<_hitVec.size();hit_i++) {
    CalorimeterHit* calohit=_hitVec[hit_i]->getCaloHit();
    EVENT::LCObjectVec simhits = nav.getRelatedToObjects(calohit);
    if (simhits.size()>0) {
      for (size_t iv=0; iv<simhits.size(); iv++) {
        SimCalorimeterHit* simhit = dynamic_cast <SimCalorimeterHit*> (simhits[iv]);
        for (int im=0; im<simhit->getNMCContributions(); im++) {
          MCParticle* mcp = simhit->getParticleCont(im);
          float en = simhit->getEnergyCont(im);

          MCParticle* gen = lastGeneratedMCAncestor(mcp);
          if (gen) {
            int pdg = abs(gen->getPDG());
            if ( pdgContribsGen.find(pdg)!=pdgContribsGen.end() ) pdgContribsGen[pdg]+=en;
            else pdgContribsGen[pdg]=en;
          }

          MCParticle* cal = firstDecInCaloMCAncestor(mcp);
          if (cal) {
            int pdg = abs(cal->getPDG());
            if ( pdgContribsCal.find(pdg)!=pdgContribsCal.end() ) pdgContribsCal[pdg]+=en;
            else pdgContribsCal[pdg]=en;
          }

        }
      }
    }
  }

  float totgenen(0);
  for ( std::map < int, float >::iterator ij=pdgContribsGen.begin(); ij!=pdgContribsGen.end(); ij++) totgenen+=ij->second;
  for ( std::map < int, float >::iterator ij=pdgContribsGen.begin(); ij!=pdgContribsGen.end(); ij++) {
    float frac = ij->second/totgenen;
    if (frac>_MCGenFrac[0]) {
      _MCGenFrac[1] = _MCGenFrac[0];
      _MCGenPdg[1]  = _MCGenPdg[0];
      _MCGenFrac[0]=frac;
      _MCGenPdg[0]=ij->first;
    } else if (frac>_MCGenFrac[1]) {
      _MCGenFrac[1]=frac;
      _MCGenPdg[1]=ij->first;
    }
  }

  float totcalen(0);
  for ( std::map < int, float >::iterator ij=pdgContribsCal.begin(); ij!=pdgContribsCal.end(); ij++) totcalen+=ij->second;
  for ( std::map < int, float >::iterator ij=pdgContribsCal.begin(); ij!=pdgContribsCal.end(); ij++) {
    float frac = ij->second/totcalen;
    if (frac>_MCCalFrac[0]) {
      _MCCalFrac[1]=_MCCalFrac[0];
      _MCCalPdg[1]=_MCCalPdg[0];
      _MCCalFrac[0]=frac;
      _MCCalPdg[0]=ij->first;
    } else if (frac>_MCCalFrac[1]) {
      _MCCalFrac[1]=frac;
      _MCCalPdg[1]=ij->first;
    }
  }

  return;
}


MCParticle* GarlicExtendedCluster::lastGeneratedMCAncestor( MCParticle* mcp ) {
  MCParticle* genAnc(0);
  MCParticleVec ancs = allMCAncestors(mcp);
  for (size_t i=0; i<ancs.size(); i++) {
    MCParticle* mc = ancs[i];
    if (!mc->isCreatedInSimulation()) {
      genAnc=mc;
      break;
    }
  }
  return genAnc;
}


MCParticle* GarlicExtendedCluster::firstDecInCaloMCAncestor( MCParticle* mcp ) {
  MCParticle* lastDecInCalo(0);
  MCParticleVec ancs = allMCAncestors(mcp);
  for (size_t i=0; i<ancs.size(); i++) {
    MCParticle* mc = ancs[i];
    if (mc->isDecayedInCalorimeter()) {
      lastDecInCalo=mc;
    } else {
      break;
    }
  }
  return lastDecInCalo;
}

MCParticleVec GarlicExtendedCluster::generatedDecInCaloMCAncestors( MCParticle* mcp ) {
  MCParticleVec ancs = allMCAncestors(mcp);
  MCParticleVec genDecInCalo;
  for (size_t i=0; i<ancs.size(); i++) {
    MCParticle* mc = ancs[i];
    if (!mc->isCreatedInSimulation() && mc->isDecayedInCalorimeter()) genDecInCalo.push_back(mc);
  }
  return genDecInCalo;
}


MCParticleVec GarlicExtendedCluster::allMCAncestors( MCParticle* mcp ) {
  MCParticleVec ancs;
  MCParticleVec temp;
  temp.push_back(mcp);
  while (1) {
    MCParticleVec temp2;
    for (size_t i=0; i<temp.size(); i++) {
      MCParticle* mc = temp[i];
      int pdg = abs(mc->getPDG());
      if ( ! (pdg<=6 || pdg==92 ) ) {
        ancs.push_back(mc);
        MCParticleVec pars = mc->getParents();
        for (size_t j=0; j<pars.size(); j++) {
	  int genstatus = pars[j]->getGeneratorStatus();
	  if ( genstatus>0 && genstatus<3 ) // added by DJ, dec 2015. otherwise can produce long (inf?) loops in some cases
	    temp2.push_back(pars[j]);
        }
      }
    }
    temp=temp2;
    if (temp.size()==0) break;
  }
  return ancs;
}



IMPL::LCGenericObjectImpl GarlicExtendedCluster::getGenericObject() {

  GarlicClusterVarsGenericObject obj;

  obj.setID    (getID());
  obj.setNHITS (getNhits());
  obj.setZONE  (getZone());
  obj.setNNSEL (getNNsel());

  obj.setENERGY              (getEnergy());
  obj.setTHETA               (getTheta());
  obj.setPHI                 (getPhi());
  obj.setTRKDIST_COG         (getTrackDist_cog());
  obj.setTRKDIST_PROJ        (getTrackDist_proj());
  obj.setTRKDIST_MIN         (getTrackDist_min());
  obj.setTRKDIST_FIRST       (getTrackDist_first());

  obj.setPOINTING_ANGLE      (getClusterPointAngle());
  obj.setECCENTRICITY        (getEccentricity());
  obj.setWIDTH               (getWidth());
  obj.setVOLUME              (getVolume());
  obj.setSTART               (getStart());
  obj.setEND                 (getEnd());
  obj.setMEAN_DEPTH          (getMeanDepth());
  obj.setREL_MEAN_DEPTH      (getRelMeanDepth());
  obj.setHITEN_MEAN          (getHitMeanEn());
  obj.setHITEN_RMS           (getHitRMSEn());
  obj.setHITEN_Q1            (getHitQ1En());
  obj.setHITEN_Q2            (getHitQ2En());
  obj.setHITEN_Q3            (getHitQ3En());
  obj.setNNOUT               (getNNval());
  obj.setTRACKANGLE_PROJ     (getTrackAng_proj());


  float* f = getFractalDimension();
  obj.setFRACDIM_2           (f[0]);
  obj.setFRACDIM_4           (f[1]);
  obj.setFRACDIM_8           (f[2]);

  std::pair <float, float> ff = getTransverseAxisLengths();
  obj.setTRANSAXISLENGTH_MIN (ff.first);
  obj.setTRANSAXISLENGTH_MAX (ff.second);

  ff = getTransverseRMS();
  obj.setTRANSVERSERMS_MIN   (ff.first);
  obj.setTRANSVERSERMS_MAX   (ff.second);

  ff = get1dMol60();
  obj.setMOL60_MIN           (ff.first);
  obj.setMOL60_MAX           (ff.second);
  ff = get1dMol80();
  obj.setMOL80_MIN           (ff.first);
  obj.setMOL80_MAX           (ff.second);
  ff = get1dMol90();
  obj.setMOL90_MIN           (ff.first);
  obj.setMOL90_MAX           (ff.second);
  ff = get1dMol95();
  obj.setMOL95_MIN           (ff.first);
  obj.setMOL95_MAX           (ff.second);

  f = getTubeEn();
  obj.setTUBE_EN(f);

  f = getTubeN();
  obj.setTUBE_N(f);

  f = getLongEn();
  obj.setLONG_EN(f);

  f = getLongN();
  obj.setLONG_N(f);

  f = getRelLongEn();
  obj.setRELLONG_EN(f);

  f = getRelRelLongEn();
  obj.setRELRELLONG_EN(f);

  f = getRelLongN();
  obj.setRELLONG_N(f);

  // added extra vars: DJ 2016/03/29

  float eonp = getAssociatedTrack() ? getEnergy()/getAssociatedTrack()->getTotalMomentum() : 0 ;
  obj.setELECTRON_EONP     (  eonp  );

  obj.setREL_MEAN_DEPTH    (getRelMeanDepth()  );
  obj.setEARLY1DMOL90_MIN  (getEarly1dMol90().first  );
  obj.setEARLY1DMOL90_MAX  (getEarly1dMol90().second  );
  obj.setFRACPLAYERS       (getFracPseudoLayers()  );
  obj.setMAXPLAYERHOLE     (getBiggestPlayHole()  );
  obj.setCLMASS            (getClusterMass()  );
  obj.setNLAYERS           (getNLayers()  );
  obj.setNPLAYERS          (getNPseudoLayers() );



  return obj.getGenericObject();

}

void GarlicExtendedCluster::Print() {

  cout << " _nhits " << _nhits;
  cout << " _energy " << _energy;
  cout << " _theta " << _theta;
  cout << " _phi " << _phi;
  cout << " _zone " << _zone;
  cout << " _trackDist_cog " << _trackDist_cog;
  cout << " _trackDist_proj " << _trackDist_proj;
  cout << " _trackDist_min " << _trackDist_min;
  cout << " _trackDist_first " << _trackDist_first;
  cout << " _trackAng_proj " << _trackAng_proj;
  cout << " _clusterPointingAngle " << _clusterPointingAngle;
  cout << " _fitted_eccentricity " << _fitted_eccentricity;
  cout << " _fitted_width " << _fitted_width;
  cout << " _fitted_volume " << _fitted_volume;
  cout << " _start " << _start;
  cout << " _end " << _end;
  cout << " _meandepth " << _meandepth;

  cout << " _tube_eFracs " ;
  for (int j=0; j<ntubes; j++) cout << _tube_eFracs[j] << " ";
  cout << endl;

  cout << " _tube_nFracs " ;
  for (int j=0; j<ntubes; j++) cout << _tube_nFracs[j] << " ";
  cout << endl;

  cout << " _long_eFracs ";
  for (int j=0; j<nlong; j++) cout << _long_eFracs[j] << " ";
  cout << endl;

  cout << " _hitMeanEn " <<  _hitMeanEn << endl;
  cout << " _hitRMSEn  " <<  _hitRMSEn  << endl;
  cout << " _hitQ1En   " <<  _hitQ1En   << endl;
  cout << " _hitQ3En   " <<  _hitQ3En   << endl;
}

std::pair < float, float > GarlicExtendedCluster::getCombinedInvariantMass( GarlicExtendedCluster* cl2 ) {
  // calculate invariant mass of this cluster with another one

  float dot(0);
  float mod1(0);
  float mod2(0);
  for (int i=0; i<3; i++) {
    mod1+=pow(this->getCentreOfGravity()[i],2);
    mod2+=pow(cl2->getCentreOfGravity()[i],2);
    dot+=this->getCentreOfGravity()[i]*cl2->getCentreOfGravity()[i];
  }
  mod1=sqrt(mod1);
  mod2=sqrt(mod2);

  float cosTheta = dot/(mod1*mod2);  // angle between 2 directions (assume from IP)

  float E1 = this->getEnergy();
  float E2 = cl2->getEnergy();


  float mass = sqrt(2*E1*E2*(1.-cosTheta)); // the invariant mass

  // and its expected uncertainty
  float stochTerm=GarlicAlgorithmParameters::Instance().GetStochasticTerm();
  float constTerm=GarlicAlgorithmParameters::Instance().GetConstantTerm();

  float dE1 = sqrt ( pow(stochTerm*sqrt(E1) , 2) + pow( constTerm*E1 , 2) );
  float dE2 = sqrt ( pow(stochTerm*sqrt(E2) , 2) + pow( constTerm*E2 , 2) );

  float dMass =
    sqrt( (1.-cosTheta)/2. ) *
    sqrt( (E2/E1)*pow(dE1,2) + (E1/E2)*pow(dE2,2) );

  return std::pair < float, float > (mass, dMass);
}


void GarlicExtendedCluster::changed() {
  _validClusterShapes=false;
  _validLongProj=false;
  _validEMFit=false;
  _validTransClusterShapes=false;
  _validHitEnMeasures=false;
  _validTrackDistances=false;
  _validTubeFracs=false;
  _validLongShape=false;
  _validZone=false;
  _validNNResults=false;
  _validPhotonCutSelResults=false;
  _validElectronCutSelResults=false;
  //    _validHistos=false;
  _validFractal=false;
  _validMC=false;
  _validPlayerHoles=false;
  _validTrackLike=false;
  _validClusterMass=false;
  _validClusterStart=false;
  return;
}

void GarlicExtendedCluster::init() {

  _longProf=NULL;

  _caloHitSimHitRelation=NULL;

  _associatedTrack=NULL;

  changed(); // set all valids to false

  _PreCluster=NULL;
  _Ghosts=NULL;
  _clusterShapes=NULL;
  _biggestNeighbourCluster=NULL;
  _hEnergy=NULL;
  _hHits=NULL;
  _validDirection=false;
  for (int j=0; j<3; j++) {
    _refPoint[j]=0;
    _refDirection[j]=0;
    _xprimeaxis[j]=0;
    _yprimeaxis[j]=0;
    _seededFrom[j]=0;
    _dir[j]=0;
    _seed_dir[j]=0;
  }
  _tracks = NULL;
  _rawEn=0;
  _hEnergy=NULL;
  _hHits=NULL;
  _thisCounter=counter++;

  _radii.clear();
  _radii.push_back(0.4);  // these are (now) Moliere radii
  _radii.push_back(0.8);
  _radii.push_back(1.2);
  _radii.push_back(1.6);
  _radii.push_back(2.0);

  _X0s.clear();
  _X0s.push_back(5.);  // these are X0
  _X0s.push_back(10.);
  _X0s.push_back(15.);
  _X0s.push_back(20.);
  _X0s.push_back(25.);

  _clusterID = -999;

  _nhits=0;
  _energy=0;
  _theta=-999;
  _phi=-999;
  _zone=-999;
  _trackDist_cog=-999;
  _trackDist_proj=-999;
  _trackDist_min=-999;
  _trackDist_first=-999;
  _trackAng_proj=-999;
  _clusterPointingAngle=-999;
  _fitted_eccentricity=0;
  _fitted_width=0;
  _fitted_volume=0;
  _start=0;
  _end=0;
  _meandepth=0;
  _relmeandepth=0;
  for (int j=0; j<ntubes; j++) _tube_eFracs[j]=0;
  for (int j=0; j<ntubes; j++) _tube_nFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_eFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_nFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_rel_eFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_rel_rel_eFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_rel_nFracs[j]=0;
  for (int j=0; j<3; j++)  _centreOfGravity[j]=0;
  for (int j=0; j<3; j++)  _projectedPosition[j]=0;
  _hitMeanEn=-999;
  _hitRMSEn=-999;
  _hitQ1En=-999;
  _hitQ3En=-999;
  _NNval=-999;
  _NNsel=false;
  for (int j=0; j<nfrac; j++) _fractalDimension[j]=0;

  for (int j=0; j<2; j++) {
    _MCCalPdg[j]=-999;
    _MCCalFrac[j]=-999;

    _MCGenPdg[j]=-999;
    _MCGenFrac[j]=-999;
  }

  for (int i=0; i<nlongfitpars; i++)
    _longFitPars[i]=-999;

  _transRmsMin=-999;
  _transRmsMax=-999;

  _minMol60=-999; _maxMol60=-999;
  _minMol80=-999; _maxMol80=-999;
  _minMol90=-999; _maxMol90=-999;
  _minMol95=-999; _maxMol95=-999;

  return;
}


float GarlicExtendedCluster::getClusterProperty( std::string variable ) {

  // get cluster properties by variable name (string)

  float value(-999);
  if      ( variable=="tube0_E"     ) value = this->getTubeEn()[0];
  else if ( variable=="tube0_N"     ) value = this->getTubeN ()[0];
  else if ( variable=="eOnP"        ) {
    value = this->getAssociatedTrack() ?
      this->getEnergy()/this->getAssociatedTrack()->getTotalMomentum() :
      0 ;}
  else if ( variable=="logPointAng" ) value = log(this->getClusterPointAngle());
  else if ( variable=="start"       ) value = this->getStart();
  else if ( variable=="reldepth"    ) value = this->getRelMeanDepth();
  else if ( variable=="relrelLongE0") value = this->getRelRelLongEn()[0];
  else if ( variable=="relrelLongE1") value = this->getRelRelLongEn()[1];
  else if ( variable=="relrelLongE2") value = this->getRelRelLongEn()[2];
  else if ( variable=="hitEnDistr"  ) value = (this->getHitQ3En() - this->getHitQ1En())/this->getHitQ2En();
  else if ( variable=="fracDim"     ) value = this->getFractalDimension()[0];
  else if ( variable=="transRmsA"   ) value = this->getTransverseRMS().first;
  else if ( variable=="transRmsB"   ) value = this->getTransverseRMS().second;
  else if ( variable=="eccen"       ) {
    float transRMSa = this->getTransverseRMS().first;
    float transRMSb = this->getTransverseRMS().second;
    value = (transRMSb-transRMSa)/(transRMSb+transRMSa);}
  else if ( variable=="molA"        ) value = this->get1dMol90().first;
  else if ( variable=="molB"        ) value = this->get1dMol90().second;
  else if ( variable=="eccenMol"    ) {
    float mol90a = this->get1dMol90().first;
    float mol90b = this->get1dMol90().second;
    value = (mol90b-mol90a)/(mol90a+mol90b) ;}
  else if ( variable=="earlyMolB"   ) value = this->getEarly1dMol90().second;
  else if ( variable=="fracPLay"    ) value = this->getFracPseudoLayers(); 
  else if ( variable=="pLayHoles"   ) value = this->getBiggestPlayHole(); 
  else if ( variable=="clMass"      ) value = this->getClusterMass(); 
  else if ( variable=="nLay"        ) value = float( this->getNLayers() );
  else if ( variable=="nPLay"       ) value = float( this->getNPseudoLayers() );
  else {
    cout << "GarlicExtendedCluster::getClusterProperty: asking for unknown variable! : " << variable << endl;
    assert(0);
  }

  return value;
}


