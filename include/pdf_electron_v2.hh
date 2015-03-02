#ifndef pdf_electron_H_ 
#define pdf_electron_H_ 

#include "TMath.h"

namespace pdf_electron { 

  float pdf_point_parM(float x) {
     return -2.1082 + x*-1.0177 + x*x*0.0114869 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.50567 + x*2.30138 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.248932 + x*-0.0892252 + x*x*0.0117142 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.340424 + x*0.0500509 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.224735 + x*0.0439793 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.42564 + x*0.284248 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.408185 + x*0.298487 + x*x*-0.0313042 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 4.59071 + x*1.50167 + x*x*-0.108325 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 15.7441 + x*-18.8572 + x*x*12.8825 + x*x*x*-2.7411 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.0404586 + x*1.17556 + x*x*-0.26009 ; 
  } 

  float pdf_EonP_parM(float x) {
     return 0.901011 + x*0.0562354 + x*x*-0.00676626 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.815343 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.736453 + x*-0.163775 ; 
  } 

  float pdf_long0_parW1(float x) {
     return 0.0949935 + x*-0.088323 + x*x*0.0461426 + x*x*x*-0.0093343 ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.100425 + x*-0.0601007 + x*x*0.0129074 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.100166 + x*-0.0737128 + x*x*0.0216683 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.313575 + x*-0.161895 + x*x*0.0303413 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return 0.112441 + x*-0.0732514 + x*x*0.0178705 ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.797915 + x*-0.335932 + x*x*0.054893 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return 2.13997 + x*-1.97741 + x*x*0.552067 ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.00471417 + x*0.163999 + x*x*-0.0541417 ; 
  } 

  float pdf_EonP_parW1(float x) {
     return 0.135139 + x*-0.103645 + x*x*0.0235199 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.771427 + x*-0.536616 + x*x*0.168699 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.07597 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 2.36815 + x*-2.20174 + x*x*0.606956 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 1.19534 + x*-1.19279 ) ; 
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

     if ( clusEn<1.) return 1;

     float x = std::log10(clusEn);
     float mean = pdf_long1_parM(x);
     float width = pdf_long1_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = relrelLongE[2]
  float PDF_long2 ( float clusEn, float v ) { 

     if ( clusEn<1.) return 1;

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
     float width = pdf_hitenwidth_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

// v = fracDim[0]
  float PDF_fracdim ( float clusEn, float v ) { 

     if ( clusEn<1.) return 1;

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

// v = clusEn/trackMom
  float PDF_EonP ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_EonP_parM(x);
     float width = pdf_EonP_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

}
#endif
