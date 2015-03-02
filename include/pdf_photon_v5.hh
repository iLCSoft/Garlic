#ifndef pdf_photon_H_ 
#define pdf_photon_H_ 
#include <TMath.h>
namespace pdf_photon { 

  float pdf_point_parM(float x) {
     return -2.14617 + x*-0.984742 ; 
  } 

  float pdf_start_parM(float x) {
     return -0.74311 + x*-0.0744177 + x*x*-0.100932 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.66467 + x*2.37509 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.243303 + x*-0.0988994 + x*x*0.0140968 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.34393 + x*0.0573117 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.227261 + x*0.0489167 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.26243 + x*0.379458 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.411095 + x*0.302701 + x*x*-0.0344615 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 3.77264 + x*2.47121 + x*x*-0.419385 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 5.09618 + x*1.73452 + x*x*-0.254474 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.654964 + x*0.424965 + x*x*-0.0298706 ; 
  } 

  float pdf_Mol90A_parM(float x) {
     return 10.5442 + x*8.07728 + x*x*-1.55205 ; 
  } 

  float pdf_Mol90B_parM(float x) {
     return 14.1711 + x*7.1429 + x*x*-1.60465 ; 
  } 

  float pdf_cylinMol_parM(float x) {
     return -6.18998 + x*-5.94764 + x*x*-3.43679 ; 
  } 

  float pdf_earlyMol90B_parM(float x) {
     return 9.82705 + x*3.15555 + x*x*-1.17559 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.841114 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.550546 ; 
  } 

  float pdf_long0_parW1(float x) {
     return exp ( -2.46868 + x*-0.572849 ) ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.100967 + x*-0.0629721 + x*x*0.0159025 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.0926679 + x*-0.0561439 + x*x*0.0147749 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return exp ( -1.59625 + x*-0.315382 ) ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return exp ( -2.18384 + x*-0.612312 ) ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return exp ( -0.533114 + x*-0.202659 ) ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( -0.258632 + x*-0.359214 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.122053 ; 
  } 

  float pdf_Mol90A_parW1(float x) {
     return 1.65936 ; 
  } 

  float pdf_Mol90B_parW1(float x) {
     return 1.81862 ; 
  } 

  float pdf_earlyMol90B_parW1(float x) {
     return 1.76195 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.510664 + x*-0.265163 + x*x*0.108586 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.22614 + x*-0.286949 + x*x*0.181098 ; 
  } 

  float pdf_hitenwidth_parW2(float x) {
     return 0.481084 + x*-0.36474 + x*x*0.100567 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return exp ( 0.559077 + x*-0.81394 ) ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 0.986644 + x*-1.08369 ) ; 
  } 

  float pdf_Mol90A_parW2(float x) {
     return 4.99923 + x*-2.80633 + x*x*0.450254 ; 
  } 

  float pdf_Mol90B_parW2(float x) {
     return 8.20686 + x*-6.64248 + x*x*1.62444 ; 
  } 

  float pdf_earlyMol90B_parW2(float x) {
     return 6.20223 + x*-4.04286 + x*x*1.06547 ; 
  } 

// v = log(pointAng)
  float PDF_point ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_point_parM(x);
     float width1 = pdf_point_parW1(x);
     float width2 = pdf_point_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = start
  float PDF_start ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float coeff = pdf_start_parM(x);
     return -coeff*exp ( coeff*v ); 
  } 

// v = reldepth
  float PDF_relmeanlong ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_relmeanlong_parM(x);
     float width1 = pdf_relmeanlong_parW1(x);
     float width2 = pdf_relmeanlong_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = relrelLongE[0]
  float PDF_long0 ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_long0_parM(x);
     float width = pdf_long0_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = relrelLongE[1]
  float PDF_long1 ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_long1_parM(x);
     float width = pdf_long1_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = relrelLongE[2]
  float PDF_long2 ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_long2_parM(x);
     float width = pdf_long2_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = (hitEnQ3-hitEnQ1)/hitEnQ2
  float PDF_hitenwidth ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_hitenwidth_parM(x);
     float width1 = pdf_hitenwidth_parW1(x);
     float width2 = pdf_hitenwidth_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = fracDim[0]
  float PDF_fracdim ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_fracdim_parM(x);
     float width = pdf_fracdim_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = transRMSa
  float PDF_widthMin ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_widthMin_parM(x);
     float width1 = pdf_widthMin_parW1(x);
     float width2 = pdf_widthMin_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = transRMSb
  float PDF_widthMax ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_widthMax_parM(x);
     float width1 = pdf_widthMax_parW1(x);
     float width2 = pdf_widthMax_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = -log10( 0.0001+ (transRMSb-transRMSa)/(transRMSb+transRMSa) )
  float PDF_cylin ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mpv = pdf_cylin_parM(x);
     float width = pdf_cylin_parW1(x);
     return TMath::Landau(v, mpv, width, true); 
  } 

// v = mol90a
  float PDF_Mol90A ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_Mol90A_parM(x);
     float width1 = pdf_Mol90A_parW1(x);
     float width2 = pdf_Mol90A_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = mol90b
  float PDF_Mol90B ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_Mol90B_parM(x);
     float width1 = pdf_Mol90B_parW1(x);
     float width2 = pdf_Mol90B_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = (mol90b-mol90a)/(mol90a+mol90b)
  float PDF_cylinMol ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float coeff = pdf_cylinMol_parM(x);
     return -coeff*exp ( coeff*v ); 
  } 

// v = earlymol90b
  float PDF_earlyMol90B ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_earlyMol90B_parM(x);
     float width1 = pdf_earlyMol90B_parW1(x);
     float width2 = pdf_earlyMol90B_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

}
#endif
