#include "ECALGarlicExtendedCluster.hh"
#include "ECALGarlicExtendedHit.hh"
#include "ECALGarlicClusterVarsGenericObject.hh"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <assert.h>

#include <UTIL/LCRelationNavigator.h>
#include "EVENT/SimCalorimeterHit.h"

#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryParameters.hh"
#include "ECALGarlicNeuralNetwork.hh"

using std::vector;
using std::min;
using std::max;

#include "TVector3.h"
#include "TF1.h"
#include "TFile.h"

#include "CLHEP/Vector/ThreeVector.h"

using std::cout;
using std::endl;

int ExtendedCluster2::counter=0;

void ExtendedCluster2::makeLongitudinalProjectionGraph() {

  if (_validLongProj) return;
  else _validLongProj=true;

  if (_hitVec.size()>0) {
    // group hits into pseudolayers
    std::map <int, std::vector <ExtendedHit2*> > playerHits;
    for (size_t ih=0; ih<_hitVec.size(); ih++) {
      ExtendedHit2* eh=_hitVec[ih];
      int play = eh->getPseudoLayer();

      if (playerHits.find(play)==playerHits.end()) {
	std::vector <ExtendedHit2*> temp;
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

    for (std::map <int, std::vector <ExtendedHit2*> >::iterator tt=playerHits.begin(); tt!=playerHits.end(); tt++) {

      float etot(0);
      float meanPos[3]={0,0,0};
      for (size_t i=0; i<tt->second.size(); i++) {
	ExtendedHit2* ehit = (tt->second)[i];
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


void ExtendedCluster2::calculateTransverseProjectionShapes() {

  if (_validTransClusterShapes) return;
  else _validTransClusterShapes=true;

  calculateClusterShapes();

  _transProjAxisLengths[0]=-999;
  _transProjAxisLengths[1]=-999;

  if (this->getHits()->size()>0) {
    std::pair < TH2F*, TH2I* > projs = this->GetProjectionHistos(); // use hits in all layers
    TH2F* hh = projs.first;
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
	
      _transProjAxisLengths[0] = min(p1, p2);
      _transProjAxisLengths[1] = max(p1, p2);

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

      _transRmsMin = std::min(h_xp.GetRMS(), h_yp.GetRMS());
      _transRmsMax = std::max(h_xp.GetRMS(), h_yp.GetRMS());

      float mol1 = get1dmoliere(&h_xp, 0.95);
      float mol2 = get1dmoliere(&h_yp, 0.95);
      _minMol95 = std::min( mol1, mol2 );
      _maxMol95 = std::max( mol1, mol2 );

      mol1 = get1dmoliere(&h_xp, 0.9);
      mol2 = get1dmoliere(&h_yp, 0.9);
      _minMol90 = std::min( mol1, mol2 );
      _maxMol90 = std::max( mol1, mol2 );

      mol1 = get1dmoliere(&h_xp, 0.8);
      mol2 = get1dmoliere(&h_yp, 0.8);
      _minMol80 = std::min( mol1, mol2 );
      _maxMol80 = std::max( mol1, mol2 );

      mol1 = get1dmoliere(&h_xp, 0.6);
      mol2 = get1dmoliere(&h_yp, 0.6);
      _minMol60 = std::min( mol1, mol2 );
      _maxMol60 = std::max( mol1, mol2 );

      if (p1<0 || p2<0) {
	cout << "WARNING strange transverse profile: ";
	cout << this->getHits()->size() << " rmsx,y = " << rmsx << " " << rmsy << 
	  " correl = " << cor << " , axis lengths = " << p1 << " " << p2 << endl;
      }

    }
  }

  return;
}


float ExtendedCluster2::get1dmoliere(TH1F* h, float fraction) {
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


void ExtendedCluster2::calculateFractalDimension() {

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

      int M = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["M"];
      int S = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["S-1"];
      int K = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["K-1"];

      int I = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["I"];
      int J = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (hit)["J"];

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

void ExtendedCluster2::calculateEMFit() {

  makeLongitudinalProjectionGraph();

  return;
}


void ExtendedCluster2::calculateClusterShapes() {

  if (_hitVec.size()>MAXHITS) cout << "WARNING, ignoring hits: nhits =" << _hitVec.size() << " MAXHITS=" << MAXHITS << endl;

  if (_hitVec.size()==0) {
    //    cout << "WARNING, cluster with only " << _hitVec.size() << " hits!" << endl;
    return;
  }

  if (_validClusterShapes) return;
  _validClusterShapes=true;

  float hitE[MAXHITS]={0};
  float hitX[MAXHITS]={0};
  float hitY[MAXHITS]={0};
  float hitZ[MAXHITS]={0};
 
  int nhit(0);
  _energy=0;
  _nhits=_hitVec.size();
  for (size_t i=0; i<_hitVec.size(); i++) {
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


  if (_clusterShapes) {delete _clusterShapes; _clusterShapes=NULL;}
  _clusterShapes = new ClusterShapes(nhit, hitE, hitX, hitY, hitZ);
  _clusterPointingAngle = angle(_clusterShapes->getEigenVecInertia(), _clusterShapes->getCentreOfGravity());

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
  // projection onto front ECAL face
  getFrontFaceProj(_clusterShapes->getCentreOfGravity(), _projectedPosition);

  return;
}

void ExtendedCluster2::calculateTrackDistances() {

  if (_validTrackDistances) return;
  _validTrackDistances=true;

  _trackDist_cog=999;
  _trackDist_proj=999;
  _trackDist_min=999;
  _trackDist_first=999;

  ExtendedTrack* closestTrack(0);

  if (_tracks && _tracks->size()>0) {
    for(size_t t_i=0; t_i<_tracks->size(); t_i++) {
      ExtendedTrack *a_track = dynamic_cast<ExtendedTrack*> (_tracks->at(t_i) );
      float dist = a_track->getDistanceToPoint( _centreOfGravity );
      _trackDist_cog=min(_trackDist_cog, dist);
      dist = a_track->getDistanceToPoint( _projectedPosition);
      if (dist<_trackDist_proj) {
	_trackDist_proj=dist;
	closestTrack=a_track;
      }
    }

    for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
      ExtendedHit2*   myHit = _hitVec[ihit];
      CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());
      for(size_t t_i=0; t_i<_tracks->size(); t_i++) {
        ExtendedTrack *a_track = dynamic_cast<ExtendedTrack*> (_tracks->at(t_i) );
        float pos[3];
        for (int i=0; i<3; i++) pos[i]=a_hit->getPosition()[i];
        float dist = a_track->getDistanceToPoint(pos);
        _trackDist_min=min(_trackDist_min, dist);
      }
    }
  }

  if (closestTrack) {
    _trackAng_proj = ECALGarlicGeometryHelpers::angle(closestTrack->getEcalEntryDir(), this->getClusterShape()->getEigenVecInertia() );
  } else 
    _trackAng_proj=-999;

  return;
}



void ExtendedCluster2::calculateTubeFracs() {

  if (_validTubeFracs) return;
  _validTubeFracs=true;

  calculateClusterShapes();

  assert (_radii.size() == ntubes);

  const float origin[3]={0,0,0}; // tube passing through origin and COG

  for (size_t i=0; i<ntubes; i++) {
    _tube_eFracs[i]=0;
    _tube_nFracs[i]=0;
  }

  float etot(0);
  int ntot(0);

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {

    ExtendedHit2*   myHit = _hitVec[ihit];
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (myHit->getCaloHit());

    float hit_en(a_hit->getEnergy());

    etot+=hit_en;
    ntot++;

    float distToLine = GetDistToLine( a_hit->getPosition(), origin, _clusterShapes->getCentreOfGravity() );

    int iring(-1);

    for (size_t ir=0; ir<_radii.size(); ir++) {
      if (distToLine<_radii[ir]) {
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

void ExtendedCluster2::calculateLongShape() {
  // longitudinal shower shape

  if (_validLongShape) return;
  _validLongShape=true;

  assert (_X0s.size() == nlong);

  for (size_t i=0; i<_X0s.size(); i++) {
    _long_eFracs[i]=0;
    _long_nFracs[i]=0;
    _long_rel_eFracs[i]=0;
    _long_rel_nFracs[i]=0;
  }

  float etot(0);
  int ntot(0);

  _start=999;
  _meandepth=0;
  _end=0;

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    ExtendedHit2*   myHit = _hitVec[ihit];
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

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    ExtendedHit2*   myHit = _hitVec[ihit];
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

  } // hit loop

  if (etot>0 && ntot>0) {
    for (size_t ir=0; ir<_X0s.size(); ir++) {
      _long_eFracs[ir]/=etot;
      _long_nFracs[ir]/=ntot;

      _long_rel_eFracs[ir]/=etot;
      _long_rel_nFracs[ir]/=ntot;
    }
  } else {
    cout << "etot or ntot zero ! " << etot << " " << ntot << endl;
  }

  return;
}


void ExtendedCluster2::calculateZone() {

  if (_validZone) return;
  _validZone=true;

  bool hasBarrelHits = false;
  bool hasEndcapHits = false;

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (_hitVec[ihit]->getCaloHit());
    int module = (*ECALGarlicGeometryParameters::Instance().Get_defaultDecoder()) (a_hit)["M"];
    if(0 < module && module < 6) hasBarrelHits=true;
    else                         hasEndcapHits=true;
  }

  if      ( hasBarrelHits && !hasEndcapHits) _zone=CLUS2_LOCATION_BARREL;
  else if (!hasBarrelHits &&  hasEndcapHits) _zone=CLUS2_LOCATION_ENDCAP;
  else if ( hasBarrelHits &&  hasEndcapHits) _zone=CLUS2_LOCATION_OVERLAP;
  else _zone=CLUS2_LOCATION_UNKNOWN;

  return;
}


void ExtendedCluster2::calculateHitEnMeasures() {

  if (_validHitEnMeasures) return;
  _validHitEnMeasures=true;

  std::vector <float> hitenergies;
  float etot(0);
  float ensqsum(0);

  for(size_t ihit=0; ihit<_hitVec.size(); ihit++) {
    if (! _hitVec[ihit]) cout << "ERROR ExtendedCluster2::calculateHitEnMeasures could not get hit " << ihit << " " << _hitVec[ihit] << endl;
    else {
      CalorimeterHit* a_hit = dynamic_cast <CalorimeterHit*> (_hitVec[ihit]->getCaloHit());
      if (!a_hit) cout << "ERROR ExtendedCluster2::calculateHitEnMeasures could not find calohit " << ihit << " " << _hitVec[ihit] << " " << _hitVec[ihit]->getCaloHit() << endl; 
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

    iq=_hitVec.size()-qq;
    if (iq<0) iq=0;
    if (iq>= (int) _hitVec.size()) iq=_hitVec.size()-1;
    _hitQ3En = hitenergies[iq];
  } 

  return;
}

void ExtendedCluster2::calculateNNResults() {

  if (_validNNResults) return;
  _validNNResults=true;

  std::pair <float, bool> nnres = ECALGarlicNeuralNetwork::Instance().getNNoutput(this);

  _NNval=nnres.first;
  _NNsel=nnres.second;

  //  cout << "calculateNNResults: " << _NNval << " " << _NNsel << endl;
  
  return;
}


void ExtendedCluster2::freeHits() {
  for (size_t i=0; i<_hitVec.size(); i++)
    _hitVec[i]->setCluster(NULL);
  _hitVec.clear(); 
  changed();
  return;
}



void ExtendedCluster2::GetDistancesToCluster(ExtendedCluster2 *b, double *distances)
{
  double smallest_dist=9999;
  double firstHitDistance=9999;
  int firstBPSLayer=99;
  int firstAPSLayer=99;

  for(size_t b_hit_i=0; b_hit_i<b->getHits()->size(); b_hit_i++) {
    ExtendedHit2 *myBHit = b->getHits()->at(b_hit_i);

    for(size_t a_hit_i=0; a_hit_i<_hitVec.size(); a_hit_i++) {
      ExtendedHit2 *myAHit = _hitVec[a_hit_i];

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
    
float ExtendedCluster2::getDistToClusterLine(ExtendedCluster2* cl) {
  float origin[3]={0};
  return GetDistToLine(this->getCentreOfGravity(), origin, cl->getCentreOfGravity());
}



void ExtendedCluster2::MakeProjectionHistos(float encut, int pseudolayerCut) {
  // algo changed a little by daniel
  // take plane normal to cluster COG at front face

  //  if (_validHistos) return;
  //  _validHistos=true;

  if (_hEnergy) {delete _hEnergy; _hEnergy=NULL;}
  if (_hHits)   {delete _hHits;   _hHits=NULL;}

  //  const float cogProj[3] = {_Parameters.POSx, _Parameters.POSy, _Parameters.POSz};

  for (int i=0; i<3; i++)
    _origin[i] = this->getProjectedPosition() [i];

  const float zero[3]={0,0,0};
  const float zaxis[3]={0,0,1};

  // define the axes in the plane normal to cluster cog projection to front face
  _xprimeaxis[0] = _origin[1]*zaxis[2] - _origin[2]*zaxis[1];
  _xprimeaxis[1] = _origin[2]*zaxis[0] - _origin[0]*zaxis[2];
  _xprimeaxis[2] = _origin[0]*zaxis[1] - _origin[1]*zaxis[0];
  float mag(0);
  for (int i=0; i<3; i++) mag+=pow(_xprimeaxis[i], 2);
  mag=sqrt(mag);
  for (int i=0; i<3; i++) _xprimeaxis[i]/=mag;

  _yprimeaxis[0] = _origin[1]*_xprimeaxis[2] - _origin[2]*_xprimeaxis[1];
  _yprimeaxis[1] = _origin[2]*_xprimeaxis[0] - _origin[0]*_xprimeaxis[2];
  _yprimeaxis[2] = _origin[0]*_xprimeaxis[1] - _origin[1]*_xprimeaxis[0];
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
    ExtendedHit2 *a_ext_hit=_hitVec[hit_i];

    float engev = a_ext_hit->getCaloHit()->getEnergy();
    float enmip = engev*_enConvFac;

    if ( enmip<encut ) continue;
    if ( a_ext_hit->getPseudoLayer()<=pseudolayerCut ) {

      nforSeed++;
      // project onto plane
      float projection[3] = {0};
      ECALGarlicGeometryHelpers::get3dLinePlaneIntersection(_origin, _origin, zero, a_ext_hit->getCaloHit()->getPosition(), projection);

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
    float cellsize = ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

    // make bins about 1mm size
    //  int nbinsx = int ( (max_xprime-min_xprime)/cellsize ); // this was for one bin per cell
    //  int nbinsy = int ( (max_yprime-min_yprime)/cellsize );
    int nbinsx = int ( (max_xprime-min_xprime) ); // this is for one bin per mm
    int nbinsy = int ( (max_yprime-min_yprime) );

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

    if (_hEnergy) {cout << "deleting th2f " << _hEnergy->GetName() << endl; delete _hEnergy; _hEnergy=NULL;}
    TString hname = "enHist"; hname+=_thisCounter;
    _hEnergy = new TH2F(hname, hname,
			nbinsx, min_xprime-cellsize, max_xprime+cellsize, 
			nbinsy, min_yprime-cellsize, max_yprime+cellsize);  
    
    if (_hHits) {delete _hHits; _hHits=NULL;}
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

  return;
}

void ExtendedCluster2::GetGlobalPositionFromLocal(float localx, float localy, float* pos) {
  for (int i=0; i<3; i++) pos[i] = _origin[i] + localx*_xprimeaxis[i] + localy*_yprimeaxis[i];
  return;
}

IMPL::ClusterImpl ExtendedCluster2::getClusterImpl() {

  ClusterImpl cluster;
  for (size_t ih=0; ih<_hitVec.size(); ih++) {
    CalorimeterHit* calohit = _hitVec[ih]->getCaloHit();
    cluster.addHit(calohit, 1.);
  }

  cluster.setEnergy( getEnergy() );
  cluster.setPosition( getCentreOfGravity() );

  return cluster;
}


void ExtendedCluster2::calculateMCnature() {

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


MCParticle* ExtendedCluster2::lastGeneratedMCAncestor( MCParticle* mcp ) {
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


MCParticle* ExtendedCluster2::firstDecInCaloMCAncestor( MCParticle* mcp ) {
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

MCParticleVec ExtendedCluster2::generatedDecInCaloMCAncestors( MCParticle* mcp ) {
  MCParticleVec ancs = allMCAncestors(mcp);
  MCParticleVec genDecInCalo;
  for (size_t i=0; i<ancs.size(); i++) {
    MCParticle* mc = ancs[i];
    if (!mc->isCreatedInSimulation() && mc->isDecayedInCalorimeter()) genDecInCalo.push_back(mc);
  }
  return genDecInCalo;
}


MCParticleVec ExtendedCluster2::allMCAncestors( MCParticle* mcp ) {
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
	  temp2.push_back(pars[j]);
	}
      }
    }
    temp=temp2;
    if (temp.size()==0) break;
  }
  return ancs;
}



IMPL::LCGenericObjectImpl ExtendedCluster2::getGenericObject() {

  ECALGarlicClusterVarsGenericObject obj;

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
  obj.setHITEN_MEAN          (getHitMeanEn());
  obj.setHITEN_RMS           (getHitRMSEn());
  obj.setHITEN_Q1            (getHitQ1En());
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

  f = getRelLongN();
  obj.setRELLONG_N(f);

  //  obj.setLONGFITPAR(f);

  return obj.getGenericObject();


  // 
  // 
  // const int nint=4;
  // const int nfl=70;
  // 
  // IMPL::LCGenericObjectImpl obj(nint,nfl,0);
  // 
  // int i=0;
  // obj.setIntVal(i++, getID() );
  // obj.setIntVal(i++, getNhits() );
  // obj.setIntVal(i++, getZone() );
  // obj.setIntVal(i++, getNNsel() );
  // assert(i==nint);
  // 
  // i=0;
  // obj.setFloatVal(i++, getEnergy() );
  // obj.setFloatVal(i++, getTheta() );
  // obj.setFloatVal(i++, getPhi() );
  // obj.setFloatVal(i++, getTrackDist_cog() );
  // obj.setFloatVal(i++, getTrackDist_proj() );
  // obj.setFloatVal(i++, getTrackDist_min() );
  // obj.setFloatVal(i++, getTrackDist_first() );
  // obj.setFloatVal(i++, getClusterPointAngle() );
  // obj.setFloatVal(i++, getEccentricity() );
  // obj.setFloatVal(i++, getWidth() );
  // obj.setFloatVal(i++, getVolume() );
  // obj.setFloatVal(i++, getStart() );
  // obj.setFloatVal(i++, getEnd() );
  // obj.setFloatVal(i++, getMeanDepth() );
  // float* f = getTubeEn();
  // for (int j=0; j<ntubes; j++) obj.setFloatVal(i++, f[j]);
  // f = getTubeN();
  // for (int j=0; j<ntubes; j++) obj.setFloatVal(i++, f[j]);
  // f = getLongEn();
  // for (int j=0; j<nlong; j++)  obj.setFloatVal(i++, f[j]);
  // f = getLongN();
  // for (int j=0; j<nlong; j++)  obj.setFloatVal(i++, f[j]);
  // 
  // 
  // f = getFractalDimension();
  // for (int j=0; j<nfrac; j++)  obj.setFloatVal(i++, f[j]);
  // 
  // std::pair <float, float> lengths = getTransverseAxisLengths();
  // obj.setFloatVal(i++, lengths.first);
  // obj.setFloatVal(i++, lengths.second);
  // 
  // float* longf = getLongFitPars();
  // for (int jj=0; jj<5; jj++) 
  //   obj.setFloatVal(i++, longf[jj]);
  // 
  // obj.setFloatVal(i++, getTrackAng_proj() );
  // 
  // std::pair <float, float> rmsses=getTransverseRMS();
  // obj.setFloatVal(i++, rmsses.first);
  // obj.setFloatVal(i++, rmsses.second);
  // 
  // obj.setFloatVal(i++, _minMol60);
  // obj.setFloatVal(i++, _maxMol60);
  // obj.setFloatVal(i++, _minMol80);
  // obj.setFloatVal(i++, _maxMol80);
  // obj.setFloatVal(i++, _minMol90);
  // obj.setFloatVal(i++, _maxMol90);
  // obj.setFloatVal(i++, _minMol95);
  // obj.setFloatVal(i++, _maxMol95);
  // 
  // f = getRelLongEn();
  // for (int j=0; j<nlong; j++)  obj.setFloatVal(i++, f[j]);
  // f = getRelLongN();
  // for (int j=0; j<nlong; j++)  obj.setFloatVal(i++, f[j]);
  // 
  // assert(i==nfl);
  // 
  // return obj;

}




void ExtendedCluster2::init() {

  _longProf=NULL;

  _caloHitSimHitRelation=NULL;

  _validClusterShapes=false;
  _validHitEnMeasures=false;
  _validTrackDistances=false;
  _validTubeFracs=false;
  _validLongShape=false;
  _validZone=false;
  _validNNResults=false;
  _validFractal=false;
  _validMC=false;
  
  _PreCluster=NULL;
  _Ghosts=NULL;
  _clusterShapes=NULL;
  _biggestNeighbourCluster=NULL;
  _hEnergy=NULL;
  _hHits=NULL;
  for (int j=0; j<3; j++) {
    _origin[j]=0;
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
  _enConvFac=1./(41.*1.8e-4);

  _radii.clear();
  _radii.push_back(5.);
  _radii.push_back(10.);
  _radii.push_back(15.);
  _radii.push_back(20.);
  _radii.push_back(30.);

  _X0s.clear();
  _X0s.push_back(5.);
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
  for (int j=0; j<ntubes; j++) _tube_eFracs[j]=0;
  for (int j=0; j<ntubes; j++) _tube_nFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_eFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_nFracs[j]=0;
  for (int j=0; j<nlong; j++)  _long_rel_eFracs[j]=0;
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
