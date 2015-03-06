#ifndef pdf_photon_H_ 
#define pdf_photon_H_ 
#include <TMath.h>
namespace pdf_photon { 

  float pdf_point_parM(float x) {
     return -2.14888 + x*-0.986555 ; 
  } 

  float pdf_start_parM(float x) {
     return -0.743445 + x*-0.0749005 + x*x*-0.100745 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.68003 + x*2.36357 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.243568 + x*-0.0996342 + x*x*0.0144107 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.340895 + x*0.0595911 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.227032 + x*0.0490843 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.25817 + x*0.383173 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.412252 + x*0.301675 + x*x*-0.0342265 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 3.77785 + x*2.48208 + x*x*-0.427106 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 5.0987 + x*1.73391 + x*x*-0.255539 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.654518 + x*0.428625 + x*x*-0.0316169 ; 
  } 

  float pdf_Mol90A_parM(float x) {
     return 10.5772 + x*8.16522 + x*x*-1.61247 ; 
  } 

  float pdf_Mol90B_parM(float x) {
     return 14.2284 + x*7.11058 + x*x*-1.60522 ; 
  } 

  float pdf_cylinMol_parM(float x) {
     return -6.39238 + x*-5.18978 + x*x*-3.85833 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.845769 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.555916 ; 
  } 

  float pdf_long0_parW1(float x) {
     return exp ( -2.4729 + x*-0.573008 ) ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.102889 + x*-0.0689299 + x*x*0.0187016 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.09542 + x*-0.0631449 + x*x*0.0178921 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return exp ( -1.59698 + x*-0.312043 ) ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return exp ( -2.17031 + x*-0.626878 ) ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return exp ( -0.529295 + x*-0.203002 ) ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( -0.254077 + x*-0.364408 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.120792 + x*0.00367806 ; 
  } 

  float pdf_Mol90A_parW1(float x) {
     return 1.70431 ; 
  } 

  float pdf_Mol90B_parW1(float x) {
     return 1.87741 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.50823 + x*-0.25471 + x*x*0.103899 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.22359 + x*-0.285933 + x*x*0.181682 ; 
  } 

  float pdf_hitenwidth_parW2(float x) {
     return 0.476996 + x*-0.360945 + x*x*0.0997478 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return exp ( 0.559899 + x*-0.820864 ) ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 0.986685 + x*-1.08402 ) ; 
  } 

  float pdf_Mol90A_parW2(float x) {
     return 4.98054 + x*-2.91774 + x*x*0.518754 ; 
  } 

  float pdf_Mol90B_parW2(float x) {
     return 8.14156 + x*-6.58232 + x*x*1.61179 ; 
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

}
#endif
