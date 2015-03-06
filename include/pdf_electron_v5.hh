#ifndef pdf_electron_H_ 
#define pdf_electron_H_ 
#include <TMath.h>
namespace pdf_electron { 

  float pdf_point_parM(float x) {
     return -1.92231 + x*-1.56505 + x*x*0.282565 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.09487 + x*2.68544 ; 
  } 

  float pdf_long0_parM(float x) {
     return 0.261722 + x*-0.119488 + x*x*0.024791 ; 
  } 

  float pdf_long1_parM(float x) {
     return 0.323382 + x*0.0620912 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.225269 + x*0.0441979 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.44001 + x*0.278196 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.399373 + x*0.317841 + x*x*-0.0397175 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 4.88382 + x*0.155354 + x*x*1.23285 + x*x*x*-0.375347 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 14.4302 + x*-13.2909 + x*x*7.53906 + x*x*x*-1.28124 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.264581 + x*0.250401 + x*x*0.628659 + x*x*x*-0.242588 ; 
  } 

  float pdf_EonP_parM(float x) {
     return 0.859984 + x*0.212122 + x*x*-0.152529 + x*x*x*0.0395287 ; 
  } 

  float pdf_Mol90A_parM(float x) {
     return 13.057 + x*4.36179 + x*x*-0.13818 ; 
  } 

  float pdf_Mol90B_parM(float x) {
     return 37.0102 + x*-20.2835 + x*x*6.58966 ; 
  } 

  float pdf_cylinMol_parM(float x) {
     return -2.47557 + x*-2.0637 + x*x*-4.3483 ; 
  } 

  float pdf_earlyMol90B_parM(float x) {
     return 22.3545 + x*-13.0612 + x*x*4.0412 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.840689 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.806302 + x*-0.214456 ; 
  } 

  float pdf_long0_parW1(float x) {
     return exp ( -2.39561 + x*-0.695103 ) ; 
  } 

  float pdf_long1_parW1(float x) {
     return 0.105433 + x*-0.0752649 + x*x*0.0197223 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.101738 + x*-0.08003 + x*x*0.0248425 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return exp ( -1.15056 + x*-0.546427 ) ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return 0.114456 + x*-0.0732882 + x*x*0.0172019 ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return exp ( -0.271323 + x*-0.401999 ) ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return exp ( 1.08056 + x*-1.26749 ) ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.0643894 + x*-0.0503097 + x*x*0.054072 ; 
  } 

  float pdf_EonP_parW1(float x) {
     return exp ( -1.91011 + x*-1.01632 ) ; 
  } 

  float pdf_Mol90A_parW1(float x) {
     return 2.34669 + x*-0.611637 ; 
  } 

  float pdf_Mol90B_parW1(float x) {
     return 7.75231 + x*-7.2037 + x*x*1.98298 ; 
  } 

  float pdf_earlyMol90B_parW1(float x) {
     return 4.50043 + x*-3.2197 + x*x*0.755979 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.891618 + x*-0.944902 + x*x*0.366036 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.11983 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return 2.96327 + x*-3.39189 + x*x*1.09089 ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 1.49098 + x*-1.43659 ) ; 
  } 

  float pdf_Mol90A_parW2(float x) {
     return 10.8485 + x*-18.8175 + x*x*13.4327 + x*x*x*-3.29696 ; 
  } 

  float pdf_Mol90B_parW2(float x) {
     return 19.6674 + x*-40.585 + x*x*30.1938 + x*x*x*-7.30241 ; 
  } 

  float pdf_earlyMol90B_parW2(float x) {
     return 18.6235 + x*-39.7077 + x*x*30.4326 + x*x*x*-7.40409 ; 
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
