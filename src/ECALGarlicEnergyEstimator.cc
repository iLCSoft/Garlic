#include "ECALGarlicEnergyEstimator.hh"
#include <assert.h>

#include "ECALGarlicConstants.hh"
using namespace ECALGarlicConstants;

#include <vector>
using std::vector;

#include "TVector3.h"

#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryParameters.hh"

#include "ECALGarlicGeometryHelpers.hh"

#include <iostream>
using std::cout;
using std::endl;


ECALGarlicEnergyEstimator::ECALGarlicEnergyEstimator(ECALGarlicAlgorithmParameters* pars, ECALGarlicGeometryParameters* geopars) {

  _algoParams=pars;
  _geomParams=geopars;

  _geomHelper = new ECALGarlicGeometryHelpers(_algoParams, _geomParams);

  _corrTh = new TF1("corrTh","([0]*x*x*x*x + [1]*x*x + 1)",0,0.8);
  _corrThPar0 = new TF1("corrThPar1","pol1",0,500);
  _corrThPar1 = new TF1("corrThPar1","pol1",0,500);
  _corrPhi = new TF1("corrPhi","(1-[0]*exp(-0.5*TMath::Power( (x-[1])/[2] , 2) ))",0,0.8);
  _corrPhiPar0 = new TF1("corrPhiPar0","exp([0]+[1]*pow(x,[5]))+[2]*x*x+[3]*x+[4]",0,500);
  _corrPhiPar1 = new TF1("corrPhiPar1","[0]*(1-exp(-(pow(x,[3]))/[1]))+[2]*x",0,500);
  _corrPhiPar2 = new TF1("corrPhiPar2","[0]*(1-exp(-(pow(x,[3]))/[1]))+[2]*x",0,500);

  for (int i=0; i<2; i++) {
    _corrThPar0->SetParameter(i,_algoParams->GetCorrThetaPars0()[i]);
    _corrThPar1->SetParameter(i,_algoParams->GetCorrThetaPars1()[i]);
  }
  for (int i=0; i<6; i++) _corrPhiPar0->SetParameter(i, _algoParams->GetCorrPhiPars0()[i]);
  for (int i=0; i<4; i++) _corrPhiPar1->SetParameter(i, _algoParams->GetCorrPhiPars1()[i]);
  for (int i=0; i<4; i++) _corrPhiPar2->SetParameter(i, _algoParams->GetCorrPhiPars2()[i]);

  _alp = new TF1("alp","(tanh([0]*x)-[1])*(exp((-x+[2])/[3]+ exp((-x+[6])/[7])*[4]))+[5]",0,400);
  _bet = new TF1("bet","[0]*(1-exp(-(pow(x,[2]))/[1]))",0,400);
  _g = new TF1("g","[0]*pow(x+[2],[1])+[3]*x",0,400);
  _d = new TF1("d","[0]+[2]*log([1]*x+[3])+[4]*x",0,400);
  _lam = new TF1("lam","1-[0]*exp(-(pow(x,[3]))/[1])+[2]*x",0,400);
  
  par1_f = new TF1("par1_f","pol2",0,160);
  par0_f = new TF1("par0_f","pol2",0,160);
  l_corr = new TF1("l_corr","[2]/2*(TMath::TanH(([0]-x)/[1])-1)",0.,30.);

  for (uint i=0; i<_algoParams->GetPar0FPars().size(); i++) {
    par0_f->SetParameter(i, _algoParams->GetPar0FPars()[i]);
    par1_f->SetParameter(i, _algoParams->GetPar1FPars()[i]);
  }
  
  _f1_b = new TF1("F1_b","pol1",0,500);
  _f2_b = new TF1("F2_b","pol1",0,500);
  _f1_e = new TF1("F1_e","pol1",0,500);
  _f2_e = new TF1("F2_e","pol1",0,500);

  for (int i=0; i<2; i++) {
    _f1_b->SetParameter(i, _algoParams->GetF1ParsBarrel()[i]);
    _f2_b->SetParameter(i, _algoParams->GetF2ParsBarrel()[i]);
    _f1_e->SetParameter(i, _algoParams->GetF1ParsEndcap()[i]);
    _f2_e->SetParameter(i, _algoParams->GetF2ParsEndcap()[i]);
  }

  return;
}


double ECALGarlicEnergyEstimator::toGeVfunction(double x, std::vector <float> * parameters) {
  double y = 0;

  // double a = (*parameters)[2];
  // double b = (*parameters)[1];
  // double c = (*parameters)[0];
  // double alpha = (*parameters)[3];

  assert (parameters->size()==4);

  double a     = (*parameters)[2];
  double b     = (*parameters)[1];
  double c     = (*parameters)[0];
  double alpha = (*parameters)[3];

  double beta = 40*a-40*alpha+b;
  double gamma = 400*a+20*b+c-400*alpha-20*beta;
  if(x < 20)  y = a*x*x+b*x+c;
  else        y = alpha*x*x+beta*x+gamma;
  return y;
}

double ECALGarlicEnergyEstimator::getGeV(double en, bool isBarrel){

  std::vector <float> pars = isBarrel ? _algoParams->GetToGeVParsBarrel() : _algoParams->GetToGeVParsEndcap();

  return toGeVfunction(en, &pars);
}

double ECALGarlicEnergyEstimator::correctEnergyTheta(double in_en, double theta) {
  //  double in_en = clusPar->E_GeV;
  double par0 = _corrThPar0->Eval(in_en);
  double par1 = _corrThPar1->Eval(in_en);
  _corrTh->SetParameter(0,par0);
  _corrTh->SetParameter(1,par1);
  double factor = _corrTh->Eval(fabs(cos(theta)));
  double corr_en =  (in_en)/factor;
  return corr_en;
}

double ECALGarlicEnergyEstimator::correctEnergyPhi(double in_en, double phi) {
   double par0 = _corrPhiPar0->Eval(in_en);
   double par1 = _corrPhiPar1->Eval(in_en);
   double par2 = _corrPhiPar2->Eval(in_en);
   _corrPhi->SetParameter(0,par0);
   _corrPhi->SetParameter(1,par1);
   _corrPhi->SetParameter(2,par2);
   double at_phi = (1.0/1000.)*(int(1000*(phi+pi))%int(1000*2*pi/8));
   double factor = _corrPhi->Eval(at_phi);
   double corr_en = (in_en)/factor;
   return corr_en;
}


double ECALGarlicEnergyEstimator::EstimateEnergyByMix(ClusterParameters *clusPar, double in_en) {

  int zone = clusPar->zone;

  //  int iset = (zone==1 || zone==3) ? 0 : 1;

  bool isbarrel =  (zone==1 || zone==3);

  if (isbarrel) {
    for (uint i=0; i<_algoParams->GetAlphaParsBarrel().size(); i++) _alp->SetParameter(i, _algoParams->GetAlphaParsBarrel()[i]);
    for (uint i=0; i<_algoParams->GetBetaParsBarrel().size(); i++) _bet->SetParameter(i, _algoParams->GetBetaParsBarrel()[i]);
    for (uint i=0; i<_algoParams->GetGammaParsBarrel().size(); i++) _g->SetParameter(i, _algoParams->GetGammaParsBarrel()[i]);
    for (uint i=0; i<_algoParams->GetDeltaParsBarrel().size(); i++) _d->SetParameter(i, _algoParams->GetDeltaParsBarrel()[i]);
    for (uint i=0; i<_algoParams->GetLambdaParsBarrel().size(); i++) _lam->SetParameter(i, _algoParams->GetLambdaParsBarrel()[i]);
  } else {
    for (uint i=0; i<_algoParams->GetAlphaParsEndcap().size(); i++) _alp->SetParameter(i, _algoParams->GetAlphaParsEndcap()[i]);
    for (uint i=0; i<_algoParams->GetBetaParsEndcap().size(); i++) _bet->SetParameter(i, _algoParams->GetBetaParsEndcap()[i]);
    for (uint i=0; i<_algoParams->GetGammaParsEndcap().size(); i++) _g->SetParameter(i, _algoParams->GetGammaParsEndcap()[i]);
    for (uint i=0; i<_algoParams->GetDeltaParsEndcap().size(); i++) _d->SetParameter(i, _algoParams->GetDeltaParsEndcap()[i]);
    for (uint i=0; i<_algoParams->GetLambdaParsEndcap().size(); i++) _lam->SetParameter(i, _algoParams->GetLambdaParsEndcap()[i]);
  }

  double energy = 0;
  double energy_en = 0;
  double energy_hit = 0;
  double out_en = 0;
  double en1even = clusPar->En1even;
  double en1odd = clusPar->En1odd;
  double en2even = clusPar->En2even;
  double en2odd = clusPar->En2odd;
  double alpha = 0;
  double beta = 0;
  double f1 = 0;
  double f2 = 0;
  double n1even = clusPar->N1even;
  double n1odd = clusPar->N1odd;
  double n2even = clusPar->N2even;
  double n2odd = clusPar->N2odd;
  double n1 = n1even+n1odd;
  double n2 = n2even+n2odd;
  double gamma = 0;
  double delta = 0;
  double lambda = 0;
  for(int i=0;i<1;i++) {
    if(zone==1 || zone==3) {
      f1 = _f1_b->Eval(in_en);
      f2 = _f2_b->Eval(in_en);
      alpha = _alp->Eval(in_en);
      //alpha = 1;
    } else if(zone==2) {
      f1 = _f1_e->Eval(in_en);
      f2 = _f2_e->Eval(in_en);
      alpha = _alp->Eval(in_en);
      //alpha = 1;
    }
    beta = _bet->Eval(in_en);
    gamma = _g->Eval(in_en);
    if(in_en>0.5) delta = _d->Eval(in_en);
    lambda = _lam->Eval(in_en);
    if(lambda>1) lambda=1;

    energy = lambda*(alpha*(f1*en1even + (1-f1)*en1odd) + beta*(f2*en2even + (1-f2)*en2odd)) + (1-lambda)*(gamma*n1 + delta*n2);
    lambda=1;
    energy_en = lambda*(alpha*(f1*en1even + (1-f1)*en1odd) + beta*(f2*en2even + (1-f2)*en2odd)) + (1-lambda)*(gamma*n1 + delta*n2);
    lambda=0;
    energy_hit = lambda*(alpha*(f1*en1even + (1-f1)*en1odd) + beta*(f2*en2even + (1-f2)*en2odd)) + (1-lambda)*(gamma*n1 + delta*n2);

    clusPar->E_GeV_en = energy_en;
    clusPar->E_GeV_hits = energy_hit;

    if(fabs(in_en-energy)<0.001)
      break;
    in_en = energy;
  }
  out_en = energy;

  return out_en;

}

void ECALGarlicEnergyEstimator::CorrectLeakage(ClusterParameters *clusPar)
{
  double energy = clusPar->E_GeV;
  double depth = clusPar->depth;
  double corr = 0;
  double p0 = par0_f->Eval(energy);
  double p1 = par1_f->Eval(energy);
  TString test;
  test+=energy;
  l_corr->SetParameter(0,p0);
  l_corr->SetParameter(1,p1);
  l_corr->SetParameter(2,energy);
  corr = l_corr->Eval(depth);
  clusPar->Etot_g_noLC = clusPar->Etot_g;
  clusPar->Etot_g = clusPar->Etot_g-(corr*(clusPar->Etot_g));
  clusPar->E_GeV_noLC = ECALGarlicEnergyEstimator::getGeV(clusPar->Etot_g_noLC, true);
  clusPar->E_GeV = ECALGarlicEnergyEstimator::getGeV(clusPar->Etot_g, true);
}



void ECALGarlicEnergyEstimator::ApplyFullGapCorrection(LCEvent *evt, ExtendedCluster *myCluster, GhostCluster *myGhostCluster)
{ 
  vector<ExtendedHit* > *clusterHits = &(myCluster->hitVec);
  int NClusteredHits = clusterHits->size();

  double Etot=0;
  int has_b_hits = 0;
  int has_e_hits = 0;
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    float hit_en=a_hit->getEnergy();
    if(myHit->pseudoLayer==0)
      hit_en=0;
    Etot+=hit_en;
    CalorimeterHitZone hit_zone = (CalorimeterHitZone)myHit->zone;
    if(hit_zone == CALHITZONE_BARREL)
       has_b_hits = 1;
    else 
      if(hit_zone == CALHITZONE_ENDCAP || hit_zone == CALHITZONE_RING || hit_zone == CALHITZONE_OVERLAP )
	has_e_hits = 1;
  }
  double E_GeV = 0;

  bool isbarrel = !( has_e_hits==1 && has_b_hits==0 );
  E_GeV = ECALGarlicEnergyEstimator::getGeV(Etot, isbarrel);

//  if(has_e_hits==1 && has_b_hits==0 ) {
//    //    E_GeV = toGeVFctn_EC->Eval(Etot);
//    E_GeV = ECALGarlicEnergyEstimator::getGeV_EC(Etot);
//  } else {
//    //    E_GeV = toGeVFctn->Eval(Etot);
//    E_GeV = ECALGarlicEnergyEstimator::getGeV(Etot);
//  }

  ExtendedHit *refHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[0]);
  CalorimeterHit *ref_hit = dynamic_cast<CalorimeterHit* > (refHit->hit);
  TVector3 ref(ref_hit->getPosition());

  double theta_mod = 1;
  double en_mod = 1;
  if(fabs(cos(ref.Theta()))<_geomParams->Get_cosOfBarrel()) {
    theta_mod = fabs(pow(sin(ref.Theta()),4));
  }
  TF1 *f = new TF1("f","2.*exp(-x/[0])",0,500);  
  f->SetParameter(0,75);
  en_mod = f->Eval(E_GeV);  
  delete f;
  
  double mod = en_mod*theta_mod;
  if (_algoParams->GetDebug()>1) cout << "Mod " << mod << endl;
  double inter_wafer = mod*(_geomParams->Get_padSizeEcal()[0]*2*_geomParams->Get_guardringSize())/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  double inter_alveola = mod*(_geomParams->Get_padSizeEcal()[0]*(2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize()))/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  double inter_module = mod*(_geomParams->Get_padSizeEcal()[0]*(2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()))/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  //double inter_module = 0.6;
  double inter_wafer_cross = mod*(2*_geomParams->Get_guardringSize()*2*_geomParams->Get_guardringSize())/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  double inter_alveola_cross = mod*((2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize())*(2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize()))/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  double inter_module_cross = mod*((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule())*(2*_geomParams->Get_guardringSize()+_geomParams->Get_fiberSize()+_geomParams->Get_fiberSizeModule()))/(_geomParams->Get_padSizeEcal()[0]*_geomParams->Get_padSizeEcal()[0]);
  //double inter_module_cross = 0.65;
  
  vector<GhostHit*> ghostHits;
  int NGhostHits = 0;
  LCCollection *preClusColl = 0;
  preClusColl = evt->getCollection(_algoParams->GetEcalPreClusterCollectionName());

  if(_algoParams->GetDebug()>2)
    cout << "...with " << NClusteredHits << " hits" << endl;
  CellIDDecoder<CalorimeterHit> decoder(preClusColl);  
  for(int cl_hit_i=0;cl_hit_i<NClusteredHits;cl_hit_i++) {
    ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_i]);
    if(myHit->preShower) //neglect preshower hits
      continue; 
    CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit* > (myHit->hit);
    int cell_i = decoder(a_hit)["I"];
    int cell_j = decoder(a_hit)["J"];
    int module = decoder(a_hit)["M"];
    int stave = decoder(a_hit)["S-1"];
    int layer = decoder(a_hit)["K-1"];
    int cell_z = myHit->zone;
    float en=a_hit->getEnergy();
    vec3 hitPos;
    hitPos.x=a_hit->getPosition()[0];
    hitPos.y=a_hit->getPosition()[1];
    hitPos.z=a_hit->getPosition()[2];
    if(0<module && module<6) { // barrel correction: THE INDEX J IS ALONG Z

      if((cell_j+1)%_geomParams->Get_nCellsPerWafer()==0 || (cell_i+1)%_geomParams->Get_nCellsPerWafer()==0) { // search upper row and right column of each wafer (to avoid double counting)
	for(int cl_hit_j=0;cl_hit_j<NClusteredHits;cl_hit_j++) {
	  ExtendedHit *myOtherHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_j]);
	  CalorimeterHit *another_hit = dynamic_cast<CalorimeterHit* > (myOtherHit->hit);
	  if(a_hit==another_hit) 
	    continue;
	  int another_module = decoder(another_hit)["M"];
	  int another_layer = decoder(another_hit)["K-1"];
	  int another_stave = decoder(another_hit)["S-1"];
	  int another_cell_i = decoder(another_hit)["I"];
	  int another_cell_j = decoder(another_hit)["J"];
	  if(another_stave!=stave || another_layer!=layer) // only correct in same layer, not between staves
	    continue;
	  if((cell_i+1)%_geomParams->Get_nCellsPerWafer()==0) { // search upper row for column match
	    if(another_cell_i==cell_i+1 && another_cell_j==cell_j && (module==another_module)) { // column match
	      float another_en=another_hit->getEnergy();
	      double ghost_en = (en+another_en)/2;
	      vec3 anotherHitPos;
	      anotherHitPos.x=another_hit->getPosition()[0];
	      anotherHitPos.y=another_hit->getPosition()[1];
	      anotherHitPos.z=another_hit->getPosition()[2];
	      vec3 ghostPos;
	      ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	      ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	      ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	      GhostHit *myGhost=new GhostHit();
	      double type_factor = inter_wafer;
	      if(type_factor>=0.7)
		myGhost->Count=1;
	      myGhost->Stave=stave;
	      myGhost->Energy=ghost_en*type_factor;
	      myGhost->Position=ghostPos;
	      myGhost->Type=GAP_TYPE_WAFER;
	      myGhost->Layer = layer;
	      myGhost->PseudoLayer = myOtherHit->pseudoLayer;
	      ghostHits.push_back(myGhost);
	      if(_algoParams->GetDebug()>2) {
		cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
	      }
	    }
	  }
	  if((cell_j+1)%_geomParams->Get_nCellsPerWafer()==0) { // search right column for row match
	    if((another_cell_j==cell_j+1 && another_cell_i==cell_i && (module==another_module)) ||
	       ((cell_j+1)%(_geomParams->Get_nCellsPerWafer()*10)==0 && another_cell_j==0 && another_cell_i==cell_i && (module==another_module-1))) { 
	      // row match, here we have to handle spread over two modules as well
	      float another_en=another_hit->getEnergy();
	      double ghost_en = (en+another_en)/2;
	      vec3 anotherHitPos;
	      anotherHitPos.x=another_hit->getPosition()[0];
	      anotherHitPos.y=another_hit->getPosition()[1];
	      anotherHitPos.z=another_hit->getPosition()[2];
	      vec3 ghostPos;
	      ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	      ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	      ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	      GhostHit *myGhost=new GhostHit();
	      double type_factor = inter_wafer;
	      if(type_factor>=0.7)
		myGhost->Count=1;
	      myGhost->Energy=ghost_en*type_factor;
	      myGhost->Position=ghostPos;
	      myGhost->Type=GAP_TYPE_WAFER;
	      myGhost->Layer = layer;
	      myGhost->Stave=stave;
	      myGhost->PseudoLayer = myOtherHit->pseudoLayer;
	      if((cell_j+1)%(2*_geomParams->Get_nCellsPerWafer())==0) {
		type_factor = inter_alveola;
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Type=GAP_TYPE_ALVEOLA;
		myGhost->Energy=ghost_en*type_factor;
	      }
	      if((cell_j+1)%(10*_geomParams->Get_nCellsPerWafer())==0) {
		type_factor = inter_module;
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Type=GAP_TYPE_MODULE;
		myGhost->Energy=ghost_en*type_factor;
	      }	      
	      ghostHits.push_back(myGhost);
	      if(_algoParams->GetDebug()>2) {
		cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
	      }
	    }
	  }
	  // for the diagonal ones avoid double counting!!!
	  if((cell_j+1)%(_geomParams->Get_nCellsPerWafer())==0 && (cell_i+1)%(_geomParams->Get_nCellsPerWafer())==0 ) { // search diagonal UR->BL,  here we have to handle spread over two modules as well
	    if(((another_cell_j)==cell_j+1 && (another_cell_i)==cell_i+1 && (module==another_module)) ||
	       ((cell_j+1)%(10*_geomParams->Get_nCellsPerWafer())==0 && another_cell_j==0 && (another_cell_i)==cell_i+1 && (module==another_module-1)) ) {
	      float another_en=another_hit->getEnergy();
	      double ghost_en = (en+another_en)/2;
	      vec3 anotherHitPos;
	      anotherHitPos.x=another_hit->getPosition()[0];
	      anotherHitPos.y=another_hit->getPosition()[1];
	      anotherHitPos.z=another_hit->getPosition()[2];
	      vec3 ghostPos;
	      ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	      ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	      ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	      double type_factor=1;
	      if(((cell_j+1)%_geomParams->Get_nCellsPerWafer())==0)
		type_factor=inter_wafer_cross;
	      if(((cell_j+1)%(2*_geomParams->Get_nCellsPerWafer()))==0)
		type_factor=inter_alveola_cross;
	      if(((cell_j+1)%(10*_geomParams->Get_nCellsPerWafer()))==0) {
		type_factor=inter_module_cross;
	      }
	      NGhostHits = ghostHits.size();
	      bool double_count=false;
	      for(int g_i=0;g_i<NGhostHits;g_i++) {
		GhostHit* a_ghost_hit = dynamic_cast<GhostHit*>(ghostHits[g_i]);
		vec3 anotherGhostPos = a_ghost_hit->Position;
		double ghostDist = _geomHelper->Get3dDistance(&ghostPos,&anotherGhostPos);
		if(ghostDist < .2) {
		  double_count = true;
		  double averaged_en = (a_ghost_hit->Energy+ghost_en*type_factor)/2;
		  a_ghost_hit->Energy=averaged_en;
		  break;
		}
	      }
	      if(!double_count) {
		GhostHit *myGhost=new GhostHit();
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Energy=ghost_en*type_factor;
		myGhost->Position=ghostPos;
		myGhost->Type=GAP_TYPE_CROSS;
		myGhost->Layer = layer;
		myGhost->Stave=stave;
		myGhost->PseudoLayer = myOtherHit->pseudoLayer;
		ghostHits.push_back(myGhost);
		if(_algoParams->GetDebug()>2) {
		  cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		}
	      }
	    }
	  }
	  if((cell_i)%_geomParams->Get_nCellsPerWafer()==1 || (cell_j)%_geomParams->Get_nCellsPerWafer()==(_geomParams->Get_nCellsPerWafer()-1)) { // search diagonal BR->UL
	    if(((another_cell_i)==cell_i-1 && (another_cell_j)==cell_j+1 && (module==another_module)) ||
	       ((cell_j+1)%(10*_geomParams->Get_nCellsPerWafer())==0 && another_cell_j==0 && (another_cell_i)==cell_i-1 && (module==another_module-1)) ) {
	      float another_en=another_hit->getEnergy();
	      double ghost_en = (en+another_en)/2;
	      vec3 anotherHitPos;
	      anotherHitPos.x=another_hit->getPosition()[0];
	      anotherHitPos.y=another_hit->getPosition()[1];
	      anotherHitPos.z=another_hit->getPosition()[2];
	      vec3 ghostPos;
	      ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	      ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	      ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	      double type_factor=1;
	      if(((cell_j+1)%_geomParams->Get_nCellsPerWafer())==0)
		type_factor=inter_wafer_cross;
	      if(((cell_j+1)%(2*_geomParams->Get_nCellsPerWafer()))==0)
		type_factor=inter_alveola_cross;
	      if(((cell_j+1)%(10*_geomParams->Get_nCellsPerWafer()))==0) {
		type_factor=inter_module_cross;
	      }
	      NGhostHits = ghostHits.size();
	      bool double_count=false;
	      for(int g_i=0;g_i<NGhostHits;g_i++) {
		GhostHit* a_ghost_hit = dynamic_cast<GhostHit*>(ghostHits[g_i]);
		vec3 anotherGhostPos = a_ghost_hit->Position;
		double ghostDist = _geomHelper->Get3dDistance(&ghostPos,&anotherGhostPos);
		if(ghostDist < .2) {
		  double_count = true;
		  double averaged_en = (a_ghost_hit->Energy+ghost_en*type_factor)/2;
		  a_ghost_hit->Energy=averaged_en;
		  break;
		}
	      }
	      if(!double_count) {
		GhostHit *myGhost=new GhostHit();
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Energy=ghost_en*type_factor;
		myGhost->Position=ghostPos;
		myGhost->Type=GAP_TYPE_CROSS;
		myGhost->Layer = layer;
		myGhost->Stave=stave;
		myGhost->PseudoLayer = myOtherHit->pseudoLayer;
		ghostHits.push_back(myGhost);
		if(_algoParams->GetDebug()>2) {
		  cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		}
	      }
	    }
	  }
	}
      }
    }
    else { // Endcap correction, handle inter stave junctions seperately
      if(cell_z!=CALHITZONE_RING) { //exclude EcalRing
	if((cell_j+1)%_geomParams->Get_nCellsPerWafer()==0 || (cell_i+1)%_geomParams->Get_nCellsPerWafer()==0) { // search upper row and right column of each wafer (to avoid double counting)
	  for(int cl_hit_j=0;cl_hit_j<NClusteredHits;cl_hit_j++) {
	    ExtendedHit *myOtherHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_j]);
	    CalorimeterHit *another_hit = dynamic_cast<CalorimeterHit* > (myOtherHit->hit);
	    if(a_hit==another_hit) 
	      continue;
	    int another_module = decoder(another_hit)["M"];
	    int another_layer = decoder(another_hit)["K-1"];
	    int another_stave = decoder(another_hit)["S-1"];
	    int another_cell_i = decoder(another_hit)["I"];
	    int another_cell_j = decoder(another_hit)["J"];
	    if(another_module!=module || another_layer!=layer) // only correct in same layer, not between modules
	      continue;

	    if((cell_j+1)%_geomParams->Get_nCellsPerWafer()==0) { // search upper row for column match
	      if(another_cell_i==cell_i && another_cell_j==cell_j+1 && stave==another_stave && (module==another_module)) { // column match
		float another_en=another_hit->getEnergy();
		double ghost_en = (en+another_en)/2;
		vec3 anotherHitPos;
		anotherHitPos.x=another_hit->getPosition()[0];
		anotherHitPos.y=another_hit->getPosition()[1];
		anotherHitPos.z=another_hit->getPosition()[2];
		vec3 ghostPos;
		ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
		ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
		ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
		GhostHit *myGhost=new GhostHit();
		double type_factor = inter_wafer;
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Energy=ghost_en*type_factor;
		myGhost->Position=ghostPos;
		myGhost->Type=GAP_TYPE_WAFER;
		myGhost->Layer = layer;
		myGhost->Stave=stave;
		myGhost->PseudoLayer = myOtherHit->pseudoLayer;
		ghostHits.push_back(myGhost);
		//cout << "I added a ghost hit in a endcap gap!"  << endl;
		if(_algoParams->GetDebug()>2) {
		  cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		}
	      }
	    }

	    if((cell_i+1)%_geomParams->Get_nCellsPerWafer()==0) { // search right column for row match
	      if((another_cell_i==cell_i+1 && another_cell_j==cell_j && stave==another_stave && (module==another_module)) ) { // row match
		float another_en=another_hit->getEnergy();
		double ghost_en = (en+another_en)/2;
		vec3 anotherHitPos;
		anotherHitPos.x=another_hit->getPosition()[0];
		anotherHitPos.y=another_hit->getPosition()[1];
		anotherHitPos.z=another_hit->getPosition()[2];
		vec3 ghostPos;
		ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
		ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
		ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
		GhostHit *myGhost=new GhostHit();
		double type_factor = inter_wafer;
		if(type_factor>=0.7)
		  myGhost->Count=1;
		myGhost->Energy=ghost_en*type_factor;
		myGhost->Position=ghostPos;
		myGhost->Type=GAP_TYPE_WAFER;
		myGhost->Layer = layer;
		myGhost->Stave=stave;
		myGhost->PseudoLayer = myOtherHit->pseudoLayer;

		if((cell_i+1)%(2*_geomParams->Get_nCellsPerWafer())==0) {
		  type_factor = inter_alveola;
		  if(type_factor>=0.7)
		    myGhost->Count=1;
		  myGhost->Type=GAP_TYPE_ALVEOLA;
		  myGhost->Energy=ghost_en*type_factor;
		}
		ghostHits.push_back(myGhost);
		if(_algoParams->GetDebug()>2) {
		  cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		}
	      }
	    }
	    // for the diagonal ones avoid double counting!!!
	    if((cell_j+1)%_geomParams->Get_nCellsPerWafer()==0 && (cell_i+1)%_geomParams->Get_nCellsPerWafer()==0 ) { // search diagonal UR->BL
	      if(((another_cell_j)==cell_j+1 && (another_cell_i)==cell_i+1 && stave==another_stave && (module==another_module)) ) {
		float another_en=another_hit->getEnergy();
		double ghost_en = (en+another_en)/2;
		vec3 anotherHitPos;
		anotherHitPos.x=another_hit->getPosition()[0];
		anotherHitPos.y=another_hit->getPosition()[1];
		anotherHitPos.z=another_hit->getPosition()[2];
		vec3 ghostPos;
		ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
		ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
		ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
		double type_factor=1;
		if(((cell_i+1)%_geomParams->Get_nCellsPerWafer())==0)
		  type_factor=inter_wafer_cross;
		if(((cell_i+1)%(2*_geomParams->Get_nCellsPerWafer()))==0)
		  type_factor=inter_alveola_cross;
		NGhostHits = ghostHits.size();
		bool double_count=false;
		for(int g_i=0;g_i<NGhostHits;g_i++) {
		  GhostHit* a_ghost_hit = dynamic_cast<GhostHit*>(ghostHits[g_i]);
		  vec3 anotherGhostPos = a_ghost_hit->Position;
		  double ghostDist = _geomHelper->Get3dDistance(&ghostPos,&anotherGhostPos);
		  if(ghostDist < .2) {
		    double_count = true;
		    double averaged_en = (a_ghost_hit->Energy+ghost_en*type_factor)/2;
		    a_ghost_hit->Energy=averaged_en;
		    break;
		  }
		}
		if(!double_count) {
		  GhostHit *myGhost=new GhostHit();
		  if(type_factor>=0.7)
		    myGhost->Count=1;
		  myGhost->Energy=ghost_en*type_factor;
		  myGhost->Position=ghostPos;
		  myGhost->Type=GAP_TYPE_CROSS;
		  myGhost->Layer = layer;
		  myGhost->Stave=stave;
		  myGhost->PseudoLayer = myOtherHit->pseudoLayer;		  
		  ghostHits.push_back(myGhost);
		  if(_algoParams->GetDebug()>2) {
		    cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		  }
		}
	      }
	    }

	    if((cell_j)%_geomParams->Get_nCellsPerWafer()==0 && (cell_i+1)%_geomParams->Get_nCellsPerWafer()==0) { // search diagonal BR->UL
	      if(((another_cell_j)==cell_j-1 && (another_cell_i)==cell_i+1 && stave==another_stave && (module==another_module)) ) {
		float another_en=another_hit->getEnergy();
		double ghost_en = (en+another_en)/2;
		vec3 anotherHitPos;
		anotherHitPos.x=another_hit->getPosition()[0];
		anotherHitPos.y=another_hit->getPosition()[1];
		anotherHitPos.z=another_hit->getPosition()[2];
		vec3 ghostPos;
		ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
		ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
		ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
		double type_factor=1;
		if(((cell_i+1)%_geomParams->Get_nCellsPerWafer())==0)
		  type_factor=inter_wafer_cross;
		if(((cell_i+1)%(2*_geomParams->Get_nCellsPerWafer()))==0)
		  type_factor=inter_alveola_cross;
		NGhostHits = ghostHits.size();
		bool double_count=false;
		for(int g_i=0;g_i<NGhostHits;g_i++) {
		  GhostHit* a_ghost_hit = dynamic_cast<GhostHit*>(ghostHits[g_i]);
		  vec3 anotherGhostPos = a_ghost_hit->Position;
		  double ghostDist = _geomHelper->Get3dDistance(&ghostPos,&anotherGhostPos);
		  if(ghostDist < .2) {
		    double_count = true;
		    double averaged_en = (a_ghost_hit->Energy+ghost_en*type_factor)/2;
		    a_ghost_hit->Energy=averaged_en;
		    break;
		  }
		}
		if(!double_count) {
		  GhostHit *myGhost=new GhostHit();
		  if(type_factor>=0.7)
		    myGhost->Count=1;
		  myGhost->Energy=ghost_en*type_factor;
		  myGhost->Position=ghostPos;
		  myGhost->Type=GAP_TYPE_CROSS;
		  myGhost->Layer = layer;
		  myGhost->Stave=stave;
		  myGhost->PseudoLayer = myOtherHit->pseudoLayer;
		  ghostHits.push_back(myGhost);
		  if(_algoParams->GetDebug()>2) {
		    cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
		  }
		}
	      }
	    }
	  }
	}
	//additional ENDCAP inter-stave correction search along buttom row and left column 
	if(cell_j==0) {
	  for(int cl_hit_j=0;cl_hit_j<NClusteredHits;cl_hit_j++) {
	    ExtendedHit *myOtherHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_j]);
	    CalorimeterHit *another_hit = dynamic_cast<CalorimeterHit* > (myOtherHit->hit);
	    if(a_hit==another_hit) 
	      continue;
	    int another_module = decoder(another_hit)["M"];
	    int another_layer = decoder(another_hit)["K-1"];
	    int another_stave = decoder(another_hit)["S-1"];
	    int another_cell_i = decoder(another_hit)["I"];
	    //int another_cell_j = decoder(another_hit)["J"];
	    if(another_layer!=layer || another_module!=module)
	      continue;
	    if(another_stave==3)
	      another_stave=-1;
	    if(another_stave!=stave-1) // only correct in same layer, not between modules
	      continue;
	    if(another_cell_i!=0 ) // look at left column
	      continue;
	    vec3 anotherHitPos;
	    anotherHitPos.x=another_hit->getPosition()[0];
	    anotherHitPos.y=another_hit->getPosition()[1];
	    anotherHitPos.z=another_hit->getPosition()[2];
	    double dist=_geomHelper->Get3dDistance(&hitPos,&anotherHitPos);

	    double maxdist=sqrt((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+_geomParams->Get_padSizeEcal()[layer])*
				(2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+_geomParams->Get_padSizeEcal()[layer])+
				(_geomParams->Get_padSizeEcal()[layer]*_geomParams->Get_padSizeEcal()[layer]));
	    if(dist<maxdist) {     
	      float another_en=another_hit->getEnergy();
	      double ghost_en = (en+another_en)/2;
	      vec3 ghostPos;
	      ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	      ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	      ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	      GhostHit *myGhost=new GhostHit();
	      double type_factor = inter_module;
	      if(type_factor>=0.7)
		myGhost->Count=1;
	      myGhost->Energy=ghost_en*type_factor;
	      myGhost->Position=ghostPos;
	      myGhost->Type=GAP_TYPE_MODULE;
	      myGhost->Layer = layer;
	      myGhost->Stave=stave;
	      myGhost->PseudoLayer = myOtherHit->pseudoLayer;
	      ghostHits.push_back(myGhost);
	      if(_algoParams->GetDebug()>2) {
		cout << "I added a ghost hit in a module gap!"  << endl;
		cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
	      }
	    }
	  }
	}
      }
      if(cell_z==CALHITZONE_RING) {// Additional correction for Endcap Ring
	for(int cl_hit_j=0;cl_hit_j<NClusteredHits;cl_hit_j++) {
	  ExtendedHit *myOtherHit = dynamic_cast<ExtendedHit* > ((*clusterHits)[cl_hit_j]);
	  CalorimeterHit *another_hit = dynamic_cast<CalorimeterHit* > (myOtherHit->hit);
	  if(a_hit==another_hit) 
	    continue;
	  int another_module = decoder(another_hit)["M"];
	  int another_layer = decoder(another_hit)["K-1"];
	  int another_cell_i = decoder(another_hit)["I"];
	  int another_cell_z = myOtherHit->zone;
	  if(another_layer!=layer || another_module!=module || another_cell_z==CALHITZONE_ENDCAP) // only correct in same layer, not between modules and only gap between Ring and big encap pieces
	    continue;
	  if(another_cell_i!=0)
	    continue; // only look at innermost endcap row
	  vec3 anotherHitPos;
	  anotherHitPos.x=another_hit->getPosition()[0];
	  anotherHitPos.y=another_hit->getPosition()[1];
	  anotherHitPos.z=another_hit->getPosition()[2];
	  double dist=_geomHelper->Get3dDistance(&hitPos,&anotherHitPos);
	  double maxdist=sqrt((2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+_geomParams->Get_padSizeEcal()[layer])*(2*_geomParams->Get_guardringSize()+2*_geomParams->Get_fiberSizeModule()+_geomParams->Get_padSizeEcal()[layer])+(_geomParams->Get_padSizeEcal()[layer]/2*_geomParams->Get_padSizeEcal()[layer]/2))+1;
	  if(dist<maxdist) {     
	    float another_en=another_hit->getEnergy();
	    double ghost_en = (en+another_en)/2;
	    vec3 ghostPos;
	    ghostPos.x=(anotherHitPos.x+hitPos.x)/2;
	    ghostPos.y=(anotherHitPos.y+hitPos.y)/2;
	    ghostPos.z=(anotherHitPos.z+hitPos.z)/2;
	    GhostHit *myGhost=new GhostHit();
	    double type_factor = inter_module;
	    if(type_factor>=0.7)
	      myGhost->Count=1;
	    myGhost->Energy=ghost_en*type_factor;
	    myGhost->Position=ghostPos;
	    myGhost->Type=GAP_TYPE_MODULE;
	    myGhost->Layer = layer;
	    myGhost->Stave=stave;
	    myGhost->PseudoLayer = myOtherHit->pseudoLayer;
	    ghostHits.push_back(myGhost);
	    if(_algoParams->GetDebug()>2) {
	      cout << "I added a ghost hit in a module gap Ring!"  << endl;
	      cout << "Added GhostHit at " << ghostPos.x << " ," << ghostPos.y << " ," << ghostPos.z << endl;
	    }
	  }
	}
      }
    }
  }
  myGhostCluster->ghostHitVec=ghostHits;
  if(_algoParams->GetDebug()>2) {
    if(ghostHits.size()>0)
      cout << ghostHits.size() << " ghost hits added" << endl;
  }
}

