#ifndef pdf_electron_H_ 
#define pdf_electron_H_ 
#include <TMath.h>
namespace pdf_electron { 

  float pdf_point_parM(float x) {
     return -1.82173 + x*-1.69491 + x*x*0.320869 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 2.90292 + x*2.78066 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.257896 + x*-0.107367 + x*x*0.019037 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.334671 + x*0.0545814 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.224726 + x*0.0441268 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.4331 + x*0.281688 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.402055 + x*0.307049 + x*x*-0.0343074 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 5.12027 + x*-0.0988621 + x*x*1.25168 + x*x*x*-0.348418 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 13.6398 + x*-13.2767 + x*x*8.44927 + x*x*x*-1.64966 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.271722 + x*0.320087 + x*x*0.586778 + x*x*x*-0.24332 ; 
  } 

  float pdf_EonP_parM(float x) {
     return 0.849967 + x*0.212952 + x*x*-0.14159 + x*x*x*0.034969 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.828436 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.584918 + x*-0.0458606 ; 
  } 

  float pdf_long0_parW1(float x) {
     return 0.103056 + x*-0.114794 + x*x*0.0698915 + x*x*x*-0.0157246 ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.10343 + x*-0.0654347 + x*x*0.0147756 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.102231 + x*-0.0783584 + x*x*0.0236573 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.34346 + x*-0.289695 + x*x*0.157578 + x*x*x*-0.0355707 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return 0.113971 + x*-0.0753153 + x*x*0.0184896 ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.855021 + x*-0.426104 + x*x*0.0873848 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return 2.26066 + x*-2.1605 + x*x*0.61721 ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.0702091 + x*-0.0274492 + x*x*0.0380753 ; 
  } 

  float pdf_EonP_parW1(float x) {
     return 0.147859 + x*-0.121289 + x*x*0.029215 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.819359 + x*-0.586699 + x*x*0.178697 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.05151 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 2.77711 + x*-2.84402 + x*x*0.835329 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 1.24808 + x*-1.23385 ) ; 
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
     float width = pdf_hitenwidth_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
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
  float PDF_EonP ( float clusEn, float v ) { 
     float x = std::log10(clusEn);
     float mean = pdf_EonP_parM(x);
     float width = pdf_EonP_parW1(x);
     float norm = 1./sqrt(2.*acos(-1))/width;
     return norm*exp( -pow( (v-mean)/width, 2 )/2. );
  } 

}
#endif
