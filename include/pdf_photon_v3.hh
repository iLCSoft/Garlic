#ifndef pdf_photon_H_ 
#define pdf_photon_H_ 
#include <TMath.h>
namespace pdf_photon { 

  float pdf_point_parM(float x) {
     return -2.15918 + x*-0.987434 ; 
  } 

  float pdf_start_parM(float x) {
     return -0.726121 + x*-0.0742235 + x*x*-0.107614 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.62578 + x*2.37664 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.241798 + x*-0.0941728 + x*x*0.0118328 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.347734 + x*0.0549888 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.227247 + x*0.0488966 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.2654 + x*0.374043 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.411108 + x*0.305236 + x*x*-0.0359019 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 3.78799 + x*2.36838 + x*x*-0.36437 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 5.06038 + x*1.66824 + x*x*-0.204248 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.672252 + x*0.420593 + x*x*-0.0333164 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.796661 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.532708 ; 
  } 

  float pdf_long0_parW1(float x) {
     return 0.0848722 + x*-0.0701725 + x*x*0.0385635 + x*x*x*-0.00841105 ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.100372 + x*-0.0641046 + x*x*0.0167569 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.0940306 + x*-0.0584957 + x*x*0.0156629 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.20552 + x*-0.0604208 + x*x*0.0049226 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return exp ( -2.1878 + x*-0.611222 ) ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.630018 + x*0.0343713 + x*x*-0.252387 + x*x*x*0.0867921 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( -0.303536 + x*-0.344436 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.129939 + x*0.00541208 + x*x*-0.00510106 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.526202 + x*-0.230584 + x*x*0.0820722 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.22354 + x*-0.216535 + x*x*0.139406 ; 
  } 

  float pdf_hitenwidth_parW2(float x) {
     return 0.486594 + x*-0.378793 + x*x*0.106706 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 1.66023 + x*-1.18145 + x*x*0.268195 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 0.954467 + x*-1.05448 ) ; 
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

}
#endif
