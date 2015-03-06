#ifndef photonPDFs_hh
#define photonPDFs_hh

#include "TF1.h"

// some parameterised photon PDFs
// for various observables
// everything is nastily hard-coded...
// 
// daniel jeans, aug'14

class photonPDFs {

 public:
  photonPDFs();
  ~photonPDFs();

  float get_pdf_point( float clEn, float pointAng );
  float get_pdf_start( float clEn, float start );
  float get_pdf_relmeanlong( float clEn, float relmeanlong);
  float get_pdf_long2( float clEn, float relrelLongE2);
  float get_pdf_hitenmean( float clEn, float hitEnMean);
  float get_pdf_hitenwidth( float clEn, float hitEnQ1, float hitEnMean, float hitEnQ3);
  float get_pdf_fracdim( float clEn, float fracDim);
  float get_pdf_widthMax( float clEn, float transRMSb);
  float get_pdf_cylin( float clEn, float transRMSa, float transRMSb);

 private:

  TF1* get_pdffn_point( float clEn );
  TF1* get_pdffn_start( float clEn );
  TF1* get_pdffn_relmeanlong( float clEn );
  TF1* get_pdffn_long2( float clEn );
  TF1* get_pdffn_hitenmean( float clEn, bool pp=false );
  TF1* get_pdffn_hitenwidth( float clEn, bool pp=false );
  TF1* get_pdffn_fracdim( float clEn );
  TF1* get_pdffn_widthMax( float clEn );
  TF1* get_pdffn_cylin( float clEn );

  TF1* _fsplitGaus;
  TF1* _fsingleGaus;
  TF1* _fExp;
  TF1* _fLandau;

};

#endif
