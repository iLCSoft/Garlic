#ifndef pdf_photon_H_ 
#define pdf_photon_H_ 

#include "TMath.h"

namespace pdf_photon { 

  float pdf_point_parM(float x) {
     return -2.15912 + x*-0.993181 ; 
  } 

  float pdf_start_parM(float x) {
     return -0.737787 + x*-0.085326 + x*x*-0.0762639 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.62854 + x*2.3693 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.242011 + x*-0.0971392 + x*x*0.01368 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.34359 + x*0.0598989 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.226649 + x*0.0496268 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.26131 + x*0.382073 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.413172 + x*0.298564 + x*x*-0.0322242 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 3.81077 + x*2.42759 + x*x*-0.420401 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 5.07107 + x*1.79877 + x*x*-0.314245 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.669621 + x*0.416929 + x*x*-0.0270468 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.805177 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.532456 ; 
  } 

  float pdf_long0_parW1(float x) {
     return 0.0846516 + x*-0.0707912 + x*x*0.0400158 + x*x*x*-0.00885945 ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.1023 + x*-0.0715461 + x*x*0.02118 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.0957958 + x*-0.0643409 + x*x*0.0189069 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.206517 + x*-0.0399897 + x*x*-0.0119912 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return exp ( -2.1731 + x*-0.634009 ) ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.612938 + x*0.0400606 + x*x*-0.195146 + x*x*x*0.0547099 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( -0.286214 + x*-0.380176 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.127337 + x*0.00114921 + x*x*0.00275766 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.525478 + x*-0.230479 + x*x*0.0850748 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.19685 + x*-0.258038 + x*x*0.210522 ; 
  } 

  float pdf_hitenwidth_parW2(float x) {
     return 0.486388 + x*-0.424366 + x*x*0.142391 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 1.6618 + x*-1.18635 + x*x*0.268299 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 0.956238 + x*-1.06616 ) ; 
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
     float width1 = pdf_hitenwidth_parW1(x);
     float width2 = pdf_hitenwidth_parW2(x);
     float width = v<mean ? width1 : width2 ; 
     float norm = 1./sqrt(2.*acos(-1))/ ( (width1+width2)/2. );
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


}
#endif
