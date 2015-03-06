#include "photonSelector.hh"

#include <iostream>

using std::cout;
using std::endl;

int photonSelector::photon_select(ExtendedCluster2* ecl) {

  int plong = photon_longProfile(ecl);
  int trans = photon_transProfile(ecl);
  int hiten = photon_hitEnergies(ecl);
  int point = photon_pointing(ecl);

  int sel(0);
  if ( plong>1 && trans>1 && hiten>1 ) sel=2;
  else if ( plong>0 && trans>0 && hiten>0 ) sel=1;

  return sel;
}



int photonSelector::photon_longProfile(ExtendedCluster2* ecl) {
  int npassT(0);
  int npassL(0);

  int sel_start =          photon_select_start(ecl);
  int sel_reldepth =       photon_select_reldepth(ecl);
  int sel_relrelLongE0 =   photon_select_relrelLongE0(ecl);
  int sel_relrelLongE1 =   photon_select_relrelLongE1(ecl);
  int sel_relrelLongE2 =   photon_select_relrelLongE2(ecl);
  int sel_fracPLay =       photon_select_fracPLay(ecl);

  if (verbose)
    cout << "photon_longProfile " << sel_start <<  " " << sel_reldepth << " " <<
      sel_relrelLongE0 << " " << sel_relrelLongE1 << " " << sel_relrelLongE2 << endl;

  if ( sel_start       >1 ) npassT++;
  if ( sel_reldepth    >1 ) npassT++;
  if ( sel_relrelLongE0>1 ) npassT++;
  if ( sel_relrelLongE1>1 ) npassT++;
  if ( sel_relrelLongE2>1 ) npassT++;
  if ( sel_fracPLay    >1 ) npassT++;

  if ( sel_start       >0 ) npassL++;
  if ( sel_reldepth    >0 ) npassL++;
  if ( sel_relrelLongE0>0 ) npassL++;
  if ( sel_relrelLongE1>0 ) npassL++;
  if ( sel_relrelLongE2>0 ) npassL++;
  if ( sel_fracPLay    >0 ) npassL++;

  if (verbose)
    cout << "    npass (T/L) = " << npassT << " " << npassL << endl;

  int sel=0;
  if       ( npassT>=6 ) sel=2; // tight
  else if  ( (npassT>=3 && npassL>=6) || (npassT>=4 && npassL>=5) ) sel=1; // loose

  return sel;
}

int photonSelector::photon_transProfile(ExtendedCluster2* ecl) {
  int npassT(0);
  int npassL(0);

  int sel_molA      = photon_select_molA(ecl)     ;
  int sel_eccenMol  = photon_select_eccenMol(ecl) ;
  int sel_earlyMolB = photon_select_earlyMolB(ecl);

  if (verbose)
    cout << "photon_transProfile " << sel_molA << " " << sel_eccenMol << " " << sel_earlyMolB << endl;

  if ( sel_molA      >1 ) npassT++;
  if ( sel_eccenMol  >1 ) npassT++;
  if ( sel_earlyMolB >1 ) npassT++;

  if ( sel_molA      >0 ) npassL++;
  if ( sel_eccenMol  >0 ) npassL++;
  if ( sel_earlyMolB >0 ) npassL++;

  int sel(0);
  if ( npassT>=3) sel=2; // tight
  else if ( (npassL>=3 && npassT>=1) || (npassL>=2 && npassT>=2) ) sel=1; // loose;
  
  return sel;
}

int photonSelector::photon_pointing(ExtendedCluster2* ecl) {
  return photon_select_pointAng(ecl)  ;
}

int photonSelector::photon_hitEnergies(ExtendedCluster2* ecl) {
  return photon_select_hitEnDistr(ecl);
}

//-----------------------------------
// these are the detailed functions
//-----------------------------------

int photonSelector::photon_select_pointAng(ExtendedCluster2* ecl) {
  int sel=0;

  float e = log10( ecl->getEnergy() );
  float val = log(ecl->getClusterPointAngle());

  float min = -4.5566808831991166 + e*-1.091225137987069 + e*e*0.058322579109377021;
  float max = -0.86028180711597024 + e*-1.7397242513243989 + e*e*0.27511987435604651;

  if      ( val < max ) sel=2;
  else if ( val < 0.8*max ) sel=1; // loose
 
  if (verbose)
    if ( sel<2 ) 
      cout << "photon_select_pointAng " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_start(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getStart();

  float min = -1.1130427836546725e-06 + e*-2.0183867880026647e-06 + e*e*1.5833064284359502e-06;
  float max = 5.2873879406674158 + e*-0.066861781687773417 + e*e*-0.48909671653109177;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_start " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_reldepth(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getRelMeanDepth();

  float min = 2.4563891414211008 + e*2.423142985994279 + e*e*-0.015760211816292804;
  float max = 7.4911635634927505 + e*1.6152719458052429 + e*e*0.50461205373449458;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_reldepth " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_relrelLongE0(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val=ecl->getRelRelLongEn()[0];

  float min = 0.062614370953640405 + e*-0.0096959797900057886 + e*e*-0.0080761444015274818;
  float max = 0.48371799123885817 + e*-0.27969383346517207 + e*e*0.066576649580260439;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_relrelLongE0 " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_relrelLongE1(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val=ecl->getRelRelLongEn()[1];

  float min = 0.10710876444629448 + e*0.1934285892866138 + e*e*-0.032910055086320299;
  float max = 0.58990342924321426 + e*-0.10758761529120128 + e*e*0.048168069849473236;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_relrelLongE1 " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_relrelLongE2(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val=ecl->getRelRelLongEn()[2];

  float min = 0.052187338459342751 + e*0.11313428013185459 + e*e*-0.01029458182900249;
  float max = 0.5037736948686875 + e*-0.15381339270537767 + e*e*0.060324775681572053;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_relrelLongE2 " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_hitEnDistr(ExtendedCluster2* ecl) {

  float e = log10( ecl->getEnergy() );
  float val = (ecl->getHitQ3En() - ecl->getHitQ1En())/ecl->getHitQ2En();

  float min = 0.82524344911377645 + e*0.53903618136562992 + e*e*-0.041399789341445926;
  float max = 2.5646623649985836 + e*-0.5378568076332465 + e*e*0.22435386700304119;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/3. &&  val < max+(max-min)/3. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_hitEnDistr " << val << " : " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_fracDim(ExtendedCluster2* ecl) {

  float e = log10( ecl->getEnergy() );
  float val = ecl->getFractalDimension()[0];

  float min = 0.16417635147188189 + e*0.38756993930927741 + e*e*-0.040371498542121576;
  float max = 0.7438543155838111 + e*0.037809975087085446 + e*e*0.036320227697249664;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_fracDim " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_transRmsA(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getTransverseRMS().first;

  float min = 2.5281085479623751 + e*2.3036575959619414 + e*e*-0.23305718888122864;
  float max = 10.227172030677941 + e*-0.50498884359193674 + e*e*-0.055506128000675621;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_transRmsA " << val  << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_transRmsB(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getTransverseRMS().second;

  float min = 3.5084538425024543 + e*2.0543155071950205 + e*e*-0.22929956860403264;
  float max = 16.754636511463051 + e*-8.0737047106533524 + e*e*2.6310750675908277;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_transRmsB " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_eccen(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float transRMSa = ecl->getTransverseRMS().first;
  float transRMSb = ecl->getTransverseRMS().second;

  float val = (transRMSb-transRMSa)/(transRMSb+transRMSa);

  float min = 0.0099988867225685157 + e*-2.018406971846097e-06 + e*e*1.5833222614783173e-06;
  float max = 0.48326403401845441 + e*-0.41588732322577188 + e*e*0.1142589520645127;

  int sel=0;
  if      ( val < max ) sel=2; 
  else if ( val < 1.1*max ) sel=1; 

  if (verbose)
    if ( sel<2 ) 
      cout << "photon_select_eccen " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_molA(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->get1dMol90().first;

  float min = 7.5633836433795514 + e*7.1605190220100212 + e*e*-1.0182203237129848;
  float max = 27.318146779955114 + e*-4.5059491039324833 + e*e*1.5011379621761167;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  //  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 
  else if ( val > min-(max-min)/5. &&  val < max+(max-min)/10. ) sel=1;  // be more generous in narrow direction

  if (verbose)
    if (sel<2) 
      cout << "photon_select_molA " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_molB(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->get1dMol90().second;

  float min = 9.5467759937305168 + e*7.0384938598778355 + e*e*-1.0232711903020026;
  float max = 45.664157075803722 + e*-17.798623785229221 + e*e*4.3940451526558579;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_molB " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_eccenMol(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float mol90a = ecl->get1dMol90().first;
  float mol90b = ecl->get1dMol90().second;

  float val = (mol90b-mol90a)/(mol90a+mol90b) ;

  float min = 0.037498888405593407 + e*-2.0184624774090563e-06 + e*e*1.5833658023480255e-06;
  float max = 0.53604829616936078 + e*-0.40139192225563119 + e*e*0.095421964113882929;

  int sel=0;
  if      ( val < max ) sel=2; 
  else if ( val < max*1.5 ) sel=1; 

  if (verbose)
    if (sel<2)
      cout << "photon_select_eccenMol " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_earlyMolB(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getEarly1dMol90().second;

  float min = 5.9724736109172198 + e*1.4341364173082367 + e*e*0.062002036088697074;
  float max = 36.105975005210432 + e*-15.372495032231821 + e*e*10.269868025466799;

  int sel=0;
  if      ( val > min && val < max ) sel=2; 
  else if ( val > min-(max-min)/10. &&  val < max+(max-min)/10. ) sel=1; 

  if (verbose)
    if (sel<2)
      cout << "photon_select_earlyMolB " << val << " " << min << " " << max <<  " : " << sel << endl;

  return sel;
}

int photonSelector::photon_select_fracPLay(ExtendedCluster2* ecl) {
  float e = log10( ecl->getEnergy() );
  float val = ecl->getFracPseudoLayers();

  float min = 0.53534781157521982 + e*0.34178374980157655 + e*e*-0.081957404728803238;

  int sel=0;
  if ( val > 0.66 && val > min ) sel=2; 
  else if ( val > 0.9*min ) sel=1; 

  if (verbose)
    if (sel<2) 
      cout << "photon_select_fracPLay " << val << " " << min <<  " : " << sel << endl;

  return sel;
}

