#ifndef pdf_electron_H_ 
#define pdf_electron_H_ 
namespace pdf_electron { 

  float pdf_point_parM(float x) {
     return -2.08682 + x*-1.05989 + x*x*0.0291189 ; 
  } 

  float pdf_relmeanlong_parM(float x) {
     return 3.42297 + x*2.36034 ; 
  } 

  float pdf_long2_parM(float x) {
     return 0.223459 + x*0.0443988 + x*x*0.000262647 ; 
  } 

  float pdf_hitenwidth_parM(float x) {
     return 1.42377 + x*0.285431 ; 
  } 

  float pdf_fracdim_parM(float x) {
     return 0.401472 + x*0.30908 + x*x*-0.0351402 ; 
  } 

  float pdf_widthMin_parM(float x) {
     return 4.62587 + x*1.44778 + x*x*-0.0892179 ; 
  } 

  float pdf_widthMax_parM(float x) {
     return 11.3072 + x*-5.98786 + x*x*2.12259 ; 
  } 

  float pdf_cylin_parM(float x) {
     return 0.245501 + x*0.62143 ; 
  } 

  float pdf_EonP_parM(float x) {
     return 0.887581 + x*0.0754252 + x*x*-0.0132496 ; 
  } 

  float pdf_point_parW1(float x) {
     return 0.814062 ; 
  } 

  float pdf_relmeanlong_parW1(float x) {
     return 0.564074 ; 
  } 

  float pdf_long2_parW1(float x) {
     return 0.10274 + x*-0.0832399 + x*x*0.0307958 + x*x*x*-0.00255668 ; 
  } 

  float pdf_hitenwidth_parW1(float x) {
     return 0.314748 + x*-0.163827 + x*x*0.0310593 ; 
  } 

  float pdf_fracdim_parW1(float x) {
     return 0.113192 + x*-0.0744504 + x*x*0.018309 ; 
  } 

  float pdf_widthMin_parW1(float x) {
     return 0.810652 + x*-0.356368 + x*x*0.0623423 ; 
  } 

  float pdf_widthMax_parW1(float x) {
     return 2.12738 + x*-1.95997 + x*x*0.546311 ; 
  } 

  float pdf_cylin_parW1(float x) {
     return 0.040497 + x*0.0652689 + x*x*-0.00830354 ; 
  } 

  float pdf_EonP_parW1(float x) {
     return 0.140251 + x*-0.111135 + x*x*0.0260952 ; 
  } 

  float pdf_point_parW2(float x) {
     return 0.940561 + x*-0.880211 + x*x*0.31412 ; 
  } 

  float pdf_relmeanlong_parW2(float x) {
     return 1.08189 ; 
  } 

  float pdf_widthMin_parW2(float x) {
     return exp ( 1.00415 + x*-1.20581 ) ; 
  } 

  float pdf_widthMax_parW2(float x) {
     return exp ( 1.24056 + x*-1.22608 ) ; 
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
//  float PDF_start ( float clusEn, float v ) { 
//     float x = std::log10(clusEn);
//  } 

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
