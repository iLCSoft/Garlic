#include "photonPDFs.h"

#include <math.h>
#include <TMath.h>
#include <iostream>

using std::cout;
using std::endl;

Double_t photonPDFs_splitGaus(Double_t* x, Double_t* par) {
  // gaussian with different widths on left and right sides
  double mean = par[0];
  Double_t width = x[0]<mean ? par[1] : par[2];
  double norm = 1./sqrt(2.*acos(-1))/ ( (par[1]+par[2])/2. );
  return norm*exp ( -pow( (x[0]-mean)/width, 2)/2. );
}

Double_t photonPDFs_singleGaus(Double_t* x, Double_t* par) {
  // gaussian with different widths on left and right sides
  double mean = par[0];
  Double_t width = par[1];
  float norm = 1./sqrt(2.*acos(-1))/width;
  return norm*exp ( -pow( (x[0]-mean)/width, 2)/2. );
}

Double_t photonPDFs_myExp(Double_t* x, Double_t* par) {
  double coeff = par[0];
  float norm = -coeff;
  return norm*exp ( coeff*x[0]);
}

Double_t photonPDFs_myLandau(Double_t* x, Double_t* par) {
  double mpv = par[0];
  double width = par[1];
  return TMath::Landau(x[0], mpv, width, true);
}

photonPDFs::photonPDFs() {

  _fsplitGaus = new TF1("splitGaus",photonPDFs_splitGaus, -10, 15, 3);
  _fsingleGaus = new TF1("singleGaus",photonPDFs_singleGaus, -10, 10, 2);
  _fExp = new TF1("myExp",photonPDFs_myExp,0,10,1);
  _fLandau = new TF1("myLandau",photonPDFs_myLandau,-10,10,2);

}

photonPDFs::~photonPDFs() {

  delete _fsplitGaus;
  delete _fsingleGaus;
  delete _fExp;
  delete _fLandau;

}

TF1* photonPDFs::get_pdffn_point( float clEn ) {

  float mean = 
    -2.15216 
    + -0.420425*std::log(clEn);

  float width = 
    0.861474
    + -0.00962259*pow( std::log(clEn), 1 );

  float width2 = 
    0.482209
    +   -0.106342*pow( std::log(clEn), 1 )
    +     0.03398*pow( std::log(clEn), 2 )
    + -0.00416087*pow( std::log(clEn), 3 );

  _fsplitGaus->SetParameter(0, mean);
  _fsplitGaus->SetParameter(1, width);
  _fsplitGaus->SetParameter(2, width2);
  return _fsplitGaus;
}

float photonPDFs::get_pdf_point( float clEn, float pointAng ) {
  return get_pdffn_point( clEn )->Eval( std::log(pointAng) );
}


TF1* photonPDFs::get_pdffn_start( float clEn ) {

  float coeff = 
    -1.18138
    + -0.0918138*pow( std::log(clEn), 1 )
    +  0.0233499*pow( std::log(clEn), 2 )
    + -0.0135623*pow( std::log(clEn), 3 );

  _fExp->SetParameter(0, coeff);
  return _fExp;
}

float photonPDFs::get_pdf_start( float clEn, float start ) {
  return get_pdffn_start( clEn )->Eval( start );
}

TF1* photonPDFs::get_pdffn_relmeanlong( float clEn ) {

  float mean = 
    3.65621
    + 1.03553*std::log(clEn);

  float width = 
    0.524113;

  float width2 = 
    1.23401
    + -0.157015*pow( std::log(clEn), 1 )
    + 0.0493261*pow( std::log(clEn), 2 );

  _fsplitGaus->SetParameter(0, mean);
  _fsplitGaus->SetParameter(1, width);
  _fsplitGaus->SetParameter(2, width2);
  return _fsplitGaus;
}

float photonPDFs::get_pdf_relmeanlong( float clEn, float relmeanlong) {
  return get_pdffn_relmeanlong(clEn)->Eval( relmeanlong );
}


TF1* photonPDFs::get_pdffn_long2( float clEn ) {

  float mean = 
    0.21713
    +   0.0387728*pow( std::log(clEn), 1 )
    + -0.00736207*pow( std::log(clEn), 2 )
    + 0.000917662*pow( std::log(clEn), 3 );

  float width = 
    0.106584
    +  -0.0471935*pow( std::log(clEn), 1 )
    +   0.0125436*pow( std::log(clEn), 2 )
    + -0.00118624*pow( std::log(clEn), 3 );
  


  _fsingleGaus->SetParameter(0, mean);
  _fsingleGaus->SetParameter(1, width);
  return _fsingleGaus;
}

float photonPDFs::get_pdf_long2( float clEn, float relrelLongE2) {
  return get_pdffn_long2( clEn )->Eval( relrelLongE2 );
}


TF1* photonPDFs::get_pdffn_hitenmean( float clEn , bool bb ) {

  float mean = exp (
		    -3.85794e+00
		    + 2.88352e-01*std::log(clEn) );

  float width = 
    0.002476
    + -0.000182825*pow( std::log(clEn), 1 )
    + -6.82256e-05*pow( std::log(clEn), 2 )
    +  7.23914e-05*pow( std::log(clEn), 3 );

  float width2 = 
    0.00414569
    + -0.00109089*pow( std::log(clEn), 1 )
    + 0.000185521*pow( std::log(clEn), 2 )
    + 7.68093e-05*pow( std::log(clEn), 3 );


  if ( bb ) cout << "get_pdffn_hitenmean " << clEn << " " << mean << " " << width << " " << width2 << endl;

  _fsplitGaus->SetParameter(0, mean);
  _fsplitGaus->SetParameter(1, width);
  _fsplitGaus->SetParameter(2, width2);

  return _fsplitGaus;
}

float photonPDFs::get_pdf_hitenmean( float clEn, float hitEnMean) {
  return get_pdffn_hitenmean(clEn)->Eval( hitEnMean );
}

TF1* photonPDFs::get_pdffn_hitenwidth( float clEn, bool pp ) {

  float mean = 
    1.00549
    + -0.0538259*pow( std::log(clEn), 1 )
    + -0.0159273*pow( std::log(clEn), 2 )
    + 0.00113385*pow( std::log(clEn), 3 );

  float width = exp (
		     -1.68936e+00
		     + -3.60028e-01*pow( std::log(clEn), 1 ) );



  if (pp) cout << "get_pdffn_hitenwidth " << clEn << " " << mean << " " << width << endl;

  _fsingleGaus->SetParameter(0, mean);
  _fsingleGaus->SetParameter(1, width);
  return _fsingleGaus;
}

float photonPDFs::get_pdf_hitenwidth( float clEn, float hitEnQ1, float hitEnMean, float hitEnQ3) {
  return get_pdffn_hitenwidth( clEn )->Eval( (hitEnQ3 - hitEnQ1)/hitEnMean );
}

TF1* photonPDFs::get_pdffn_fracdim( float clEn) {

  float mean = 
    0.412255
    +    0.134072*pow( std::log(clEn), 1 )
    + -0.00699831*pow( std::log(clEn), 2 );

  float width = exp (
		     -2.13728e+00
		     + -2.88384e-01*pow( std::log(clEn), 1 ) );

  _fsingleGaus->SetParameter(0, mean);
  _fsingleGaus->SetParameter(1, width);
  return _fsingleGaus;
}

float photonPDFs::get_pdf_fracdim( float clEn, float fracDim) {
  return get_pdffn_fracdim(clEn)->Eval( fracDim );
}

TF1* photonPDFs::get_pdffn_widthMax( float clEn ) {

  float mean = 
    +    4.98541
    +   0.801503*pow( std::log(clEn), 1 )
    + -0.0512686*pow( std::log(clEn), 2 );

  float width = 
    0.748867
    + -0.0355234*pow( std::log(clEn), 1 )
    + -0.0363861*pow( std::log(clEn), 2 )
    + 0.00600245*pow( std::log(clEn), 3 );

  float width2 = exp (
		      9.77032e-01
		      + -4.40066e-01*pow( std::log(clEn), 1 ) );


  _fsplitGaus->SetParameter(0, mean);
  _fsplitGaus->SetParameter(1, width);
  _fsplitGaus->SetParameter(2, width2);

  return _fsplitGaus;
}

float photonPDFs::get_pdf_widthMax( float clEn, float transRMSb) {
  return get_pdffn_widthMax(clEn)->Eval( transRMSb );
}

TF1* photonPDFs::get_pdffn_cylin( float clEn ) {

  float mpv = 
    0.66779
    +    0.185063*pow( std::log(clEn), 1 )
    + -0.00721781*pow( std::log(clEn), 2 );

  float width = 
    0.129364
    +  0.00022425*pow( std::log(clEn), 1 )
    + -0.00122763*pow( std::log(clEn), 2 );


  _fLandau->SetParameter(0, mpv);
  _fLandau->SetParameter(1, width);

  return _fLandau;
}

float photonPDFs::get_pdf_cylin( float clEn, float transRMSa, float transRMSb) {
  return get_pdffn_cylin(clEn)->Eval( -log10( 0.0001+ (transRMSb-transRMSa)/(transRMSb+transRMSa) ) );
}
