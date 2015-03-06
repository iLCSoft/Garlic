#ifndef pdf_photon_H_ 
#define pdf_photon_H_ 

#include "TMath.h"

namespace pdf_photon { 

  float pdf_point_parM(float x) {
     return -2.17484 + x*-0.97848 ; 
  } 

  float pdf_start_parM(float x) {
     return -0.758381 + x*-0.0066306 + x*x*-0.122814 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.71844 + x*2.28076 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.227525 + x*0.0472342 + x*x*0.00123124 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.30705 + x*0.346405 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.413172 + x*0.298564 + x*x*-0.0322242 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 3.82361 + x*2.41784 + x*x*-0.41998 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 5.19274 + x*1.53515 + x*x*-0.193024 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.645201 + x*0.500108 + x*x*-0.0730608 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.793919 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.527727 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.096158 + x*-0.0666364 + x*x*0.0219754 + x*x*x*-0.00111725 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.201013 + x*-0.0241121 + x*x*-0.0202137 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return exp ( -2.1731 + x*-0.634009 ) ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.591753 + x*0.14842 + x*x*-0.322998 + x*x*x*0.0974327 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return 0.812224 + x*-0.225291 + x*x*-0.163726 + x*x*x*0.0902511 ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.116818 + x*0.0400788 + x*x*-0.0190351 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.488798 + x*-0.116531 + x*x*0.0242741 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.24868 + x*-0.425031 + x*x*0.302358 ; 
  } 

  float pdf_hitenwidth_parW2(float x) {
     return 0.472513 + x*-0.388024 + x*x*0.124233 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return exp ( 0.56049 + x*-0.859873 ) ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 0.963023 + x*-1.07495 ) ; 
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

// v = clusEn/trackMom
//  float PDF_EonP ( float clusEn, float v ) { 
//     float x = std::log10(clusEn);
//  } 

}
#endif
