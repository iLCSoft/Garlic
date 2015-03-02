#ifndef pdf_electron_H_ 
#define pdf_electron_H_ 
#include <TMath.h>
namespace pdf_electron { 

  float pdf_point_parM(float x) {
     return -1.92171 + x*-1.56956 + x*x*0.284581 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.1016 + x*2.70785 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.261623 + x*-0.119157 + x*x*0.0246445 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.324013 + x*0.0619685 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.225144 + x*0.0442744 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.43982 + x*0.278024 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.399444 + x*0.317312 + x*x*-0.039444 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 4.86543 + x*0.251484 + x*x*1.1484 + x*x*x*-0.355179 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 14.6762 + x*-14.4776 + x*x*8.64909 + x*x*x*-1.57089 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.260004 + x*0.260356 + x*x*0.657141 + x*x*x*-0.260693 ; 
  } 

  float pdf_EonP_parM(float x) {
     return 0.860213 + x*0.207702 + x*x*-0.148323 + x*x*x*0.0385285 ; 
  } 

  float pdf_Mol90A_parM(float x) {
     return 13.0045 + x*4.64537 + x*x*-0.286872 ; 
  } 

  float pdf_Mol90B_parM(float x) {
     return 37.1916 + x*-21.5491 + x*x*7.24992 ; 
  } 

  float pdf_cylinMol_parM(float x) {
     return -2.46245 + x*-2.17496 + x*x*-4.58057 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.836091 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.803485 + x*-0.213376 ; 
  } 

  float pdf_long0_parW1(float x) {
     return exp ( -2.39595 + x*-0.695079 ) ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.10524 + x*-0.074404 + x*x*0.0193132 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.101604 + x*-0.0793905 + x*x*0.0245262 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return exp ( -1.15072 + x*-0.54483 ) ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return 0.115038 + x*-0.0767981 + x*x*0.0189845 ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return exp ( -0.269873 + x*-0.401633 ) ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( 1.07547 + x*-1.28778 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.0636518 + x*-0.0465857 + x*x*0.0543971 ; 
  } 

  float pdf_EonP_parW1(float x) {
     return exp ( -1.91005 + x*-1.01237 ) ; 
  } 

  float pdf_Mol90A_parW1(float x) {
     return 2.35231 + x*-0.61343 ; 
  } 

  float pdf_Mol90B_parW1(float x) {
     return 7.79422 + x*-7.39087 + x*x*2.0759 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.882177 + x*-0.895506 + x*x*0.342582 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.11628 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 2.96125 + x*-3.3933 + x*x*1.09244 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 1.48699 + x*-1.44915 ) ; 
  } 

  float pdf_Mol90A_parW2(float x) {
     return 10.8014 + x*-18.5985 + x*x*13.3309 + x*x*x*-3.30286 ; 
  } 

  float pdf_Mol90B_parW2(float x) {
     return 19.1614 + x*-38.4884 + x*x*28.4465 + x*x*x*-6.89563 ; 
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
