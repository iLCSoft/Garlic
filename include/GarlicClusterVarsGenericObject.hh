#ifndef GarlicClusterVarsGenericObject_hh_
#define GarlicClusterVarsGenericObject_hh_

#include <IMPL/LCGenericObjectImpl.h>
#include <assert.h>
#include <iostream>
using std::cout;
using std::endl;

class GarlicClusterVarsGenericObject {

public:
  GarlicClusterVarsGenericObject() {
    _genObj = new LCGenericObjectImpl(NINT, NFLOAT, 0);
    for (int i=0; i<NINT; i++) _genObj->setIntVal(i, -999);
    for (int i=0; i<NFLOAT; i++) _genObj->setFloatVal(i, -999);
  }

  GarlicClusterVarsGenericObject(LCGenericObject* obj) {
    assert(obj->getNInt()==NINT);
    assert(obj->getNFloat()==NFLOAT);
    assert(obj->getNDouble()==0);
    _genObj = new LCGenericObjectImpl(NINT, NFLOAT, 0);
    for (int i=0; i<NINT; i++) _genObj->setIntVal( i, obj->getIntVal(i) );
    for (int i=0; i<NFLOAT; i++) _genObj->setFloatVal( i, obj->getFloatVal(i) );
  }

  ~GarlicClusterVarsGenericObject() {}

  IMPL::LCGenericObjectImpl getGenericObject() {return *_genObj;}

  void setID    (int ij) {_genObj->setIntVal( ID     , ij);}
  void setNHITS (int ij) {_genObj->setIntVal( NHITS  , ij);}
  void setZONE  (int ij) {_genObj->setIntVal( ZONE   , ij);}
  void setNNSEL (int ij) {_genObj->setIntVal( NNSEL  , ij);}

  void setENERGY              (float ff) {_genObj->setFloatVal( ENERGY              , ff);}
  void setTHETA               (float ff) {_genObj->setFloatVal( THETA 		    , ff);}
  void setPHI                 (float ff) {_genObj->setFloatVal( PHI 		    , ff);}
  void setTRKDIST_COG         (float ff) {_genObj->setFloatVal( TRKDIST_COG 	    , ff);}
  void setTRKDIST_PROJ        (float ff) {_genObj->setFloatVal( TRKDIST_PROJ 	    , ff);}
  void setTRKDIST_MIN         (float ff) {_genObj->setFloatVal( TRKDIST_MIN 	    , ff);}
  void setTRKDIST_FIRST       (float ff) {_genObj->setFloatVal( TRKDIST_FIRST	    , ff);}
  void setPOINTING_ANGLE      (float ff) {_genObj->setFloatVal( POINTING_ANGLE	    , ff);}
  void setECCENTRICITY        (float ff) {_genObj->setFloatVal( ECCENTRICITY 	    , ff);}
  void setWIDTH               (float ff) {_genObj->setFloatVal( WIDTH 		    , ff);}
  void setVOLUME              (float ff) {_genObj->setFloatVal( VOLUME 	            , ff);}
  void setSTART               (float ff) {_genObj->setFloatVal( START 		    , ff);}
  void setEND                 (float ff) {_genObj->setFloatVal( END 		    , ff);}
  void setMEAN_DEPTH          (float ff) {_genObj->setFloatVal( MEAN_DEPTH	    , ff);}
  void setREL_MEAN_DEPTH      (float ff) {_genObj->setFloatVal( REL_MEAN_DEPTH	    , ff);}
  void setHITEN_MEAN          (float ff) {_genObj->setFloatVal( HITEN_MEAN 	    , ff);}
  void setHITEN_RMS           (float ff) {_genObj->setFloatVal( HITEN_RMS 	    , ff);}
  void setHITEN_Q1            (float ff) {_genObj->setFloatVal( HITEN_Q1 	    , ff);}
  void setHITEN_Q2            (float ff) {_genObj->setFloatVal( HITEN_Q2 	    , ff);}
  void setHITEN_Q3            (float ff) {_genObj->setFloatVal( HITEN_Q3	    , ff);}
  void setNNOUT               (float ff) {_genObj->setFloatVal( NNOUT		    , ff);}
  void setFRACDIM_2           (float ff) {_genObj->setFloatVal( FRACDIM_2 	    , ff);}
  void setFRACDIM_4           (float ff) {_genObj->setFloatVal( FRACDIM_4 	    , ff);}
  void setFRACDIM_8           (float ff) {_genObj->setFloatVal( FRACDIM_8	    , ff);}
  void setTRANSAXISLENGTH_MIN (float ff) {_genObj->setFloatVal( TRANSAXISLENGTH_MIN , ff);}
  void setTRANSAXISLENGTH_MAX (float ff) {_genObj->setFloatVal( TRANSAXISLENGTH_MAX , ff);}
  void setTRACKANGLE_PROJ     (float ff) {_genObj->setFloatVal( TRACKANGLE_PROJ     , ff);}
  void setTRANSVERSERMS_MIN   (float ff) {_genObj->setFloatVal( TRANSVERSERMS_MIN   , ff);}
  void setTRANSVERSERMS_MAX   (float ff) {_genObj->setFloatVal( TRANSVERSERMS_MAX   , ff);}
  void setMOL60_MIN           (float ff) {_genObj->setFloatVal( MOL60_MIN 	    , ff);}
  void setMOL60_MAX           (float ff) {_genObj->setFloatVal( MOL60_MAX	    , ff);}
  void setMOL80_MIN           (float ff) {_genObj->setFloatVal( MOL80_MIN 	    , ff);}
  void setMOL80_MAX           (float ff) {_genObj->setFloatVal( MOL80_MAX	    , ff);}
  void setMOL90_MIN           (float ff) {_genObj->setFloatVal( MOL90_MIN 	    , ff);}
  void setMOL90_MAX           (float ff) {_genObj->setFloatVal( MOL90_MAX	    , ff);}
  void setMOL95_MIN           (float ff) {_genObj->setFloatVal( MOL95_MIN 	    , ff);}
  void setMOL95_MAX           (float ff) {_genObj->setFloatVal( MOL95_MAX           , ff);}


  void setPDF_POINT         (float ff) {_genObj->setFloatVal( PDF_POINT       ,ff); }
  void setPDF_START 	    (float ff) {_genObj->setFloatVal( PDF_START       ,ff); }
  void setPDF_RELMEANLONG   (float ff) {_genObj->setFloatVal( PDF_RELMEANLONG ,ff); }
  void setPDF_LONG2	    (float ff) {_genObj->setFloatVal( PDF_LONG2       ,ff); }
  void setPDF_HITENMEAN     (float ff) {_genObj->setFloatVal( PDF_HITENMEAN   ,ff); }
  void setPDF_HITENWIDTH    (float ff) {_genObj->setFloatVal( PDF_HITENWIDTH  ,ff); }
  void setPDF_FRACDIM	    (float ff) {_genObj->setFloatVal( PDF_FRACDIM     ,ff); }
  void setPDF_WIDTHMAX      (float ff) {_genObj->setFloatVal( PDF_WIDTHMAX    ,ff); }
  void setPDF_WIDTHMIN      (float ff) {_genObj->setFloatVal( PDF_WIDTHMIN    ,ff); }
  void setPDF_CYLIN	    (float ff) {_genObj->setFloatVal( PDF_CYLIN       ,ff); }
  void setPDF_LOWEN	    (float ff) {_genObj->setFloatVal( PDF_LOWEN       ,ff); }
  void setPDF_HIGHEN        (float ff) {_genObj->setFloatVal( PDF_HIGHEN      ,ff); }

  void setPDF_ELECTRON_POINT        (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_POINT      ,ff); }
  void setPDF_ELECTRON_RELMEANLONG  (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_RELMEANLONG,ff); }
  void setPDF_ELECTRON_LONG2        (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_LONG2      ,ff); }
  void setPDF_ELECTRON_HITENWIDTH   (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_HITENWIDTH ,ff); }
  void setPDF_ELECTRON_FRACDIM      (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_FRACDIM    ,ff); }
  void setPDF_ELECTRON_WIDTHMAX     (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_WIDTHMAX   ,ff); }
  void setPDF_ELECTRON_WIDTHMIN     (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_WIDTHMIN   ,ff); }
  void setPDF_ELECTRON_CYLIN        (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_CYLIN      ,ff); }
  void setPDF_ELECTRON_EONP         (float ff) {_genObj->setFloatVal(   PDF_ELECTRON_EONP       ,ff); }


  void setELECTRON_EONP     (float ff) {_genObj->setFloatVal( ELECTRON_EONP    ,ff); }
  void setEARLY1DMOL90_MIN  (float ff) {_genObj->setFloatVal( EARLY1DMOL90_MIN ,ff); }
  void setEARLY1DMOL90_MAX  (float ff) {_genObj->setFloatVal( EARLY1DMOL90_MAX ,ff); }
  void setFRACPLAYERS	    (float ff) {_genObj->setFloatVal( FRACPLAYERS      ,ff); }
  void setMAXPLAYERHOLE     (float ff) {_genObj->setFloatVal( MAXPLAYERHOLE    ,ff); }
  void setCLMASS	    (float ff) {_genObj->setFloatVal( CLMASS           ,ff); }

  void setNLAYERS	    (int ff) {_genObj->setIntVal( NLAYERS          ,ff); }
  void setNPLAYERS	    (int ff) {_genObj->setIntVal( NPLAYERS         ,ff); }


  void setTUBE_EN(float* ff) {
    _genObj->setFloatVal(TUBE_EN_0, ff[0]);
    _genObj->setFloatVal(TUBE_EN_1, ff[1]);
    _genObj->setFloatVal(TUBE_EN_2, ff[2]);
    _genObj->setFloatVal(TUBE_EN_3, ff[3]);
    _genObj->setFloatVal(TUBE_EN_4, ff[4]);
  }

  void setTUBE_N(float* ff) {
    _genObj->setFloatVal(TUBE_N_0, ff[0]);
    _genObj->setFloatVal(TUBE_N_1, ff[1]);
    _genObj->setFloatVal(TUBE_N_2, ff[2]);
    _genObj->setFloatVal(TUBE_N_3, ff[3]);
    _genObj->setFloatVal(TUBE_N_4, ff[4]);
  }

  void setLONG_EN(float* ff) {
    _genObj->setFloatVal(LONG_EN_0, ff[0]);
    _genObj->setFloatVal(LONG_EN_1, ff[1]);
    _genObj->setFloatVal(LONG_EN_2, ff[2]);
    _genObj->setFloatVal(LONG_EN_3, ff[3]);
    _genObj->setFloatVal(LONG_EN_4, ff[4]);
  }

  void setLONG_N(float* ff) {
    _genObj->setFloatVal(LONG_N_0, ff[0]);
    _genObj->setFloatVal(LONG_N_1, ff[1]);
    _genObj->setFloatVal(LONG_N_2, ff[2]);
    _genObj->setFloatVal(LONG_N_3, ff[3]);
    _genObj->setFloatVal(LONG_N_4, ff[4]);
  }

  void setRELLONG_EN(float* ff) {
    _genObj->setFloatVal(RELLONG_EN_0, ff[0]);
    _genObj->setFloatVal(RELLONG_EN_1, ff[1]);
    _genObj->setFloatVal(RELLONG_EN_2, ff[2]);
    _genObj->setFloatVal(RELLONG_EN_3, ff[3]);
    _genObj->setFloatVal(RELLONG_EN_4, ff[4]);
  }

  void setRELRELLONG_EN(float* ff) {
    _genObj->setFloatVal(RELRELLONG_EN_0, ff[0]);
    _genObj->setFloatVal(RELRELLONG_EN_1, ff[1]);
    _genObj->setFloatVal(RELRELLONG_EN_2, ff[2]);
    _genObj->setFloatVal(RELRELLONG_EN_3, ff[3]);
    _genObj->setFloatVal(RELRELLONG_EN_4, ff[4]);
  }

  void setRELLONG_N(float* ff) {
    _genObj->setFloatVal(RELLONG_N_0, ff[0]);
    _genObj->setFloatVal(RELLONG_N_1, ff[1]);
    _genObj->setFloatVal(RELLONG_N_2, ff[2]);
    _genObj->setFloatVal(RELLONG_N_3, ff[3]);
    _genObj->setFloatVal(RELLONG_N_4, ff[4]);
  }

  void setLONGFITPAR(float* ff) {
    _genObj->setFloatVal(LONGFITPAR_0, ff[0]);
    _genObj->setFloatVal(LONGFITPAR_1, ff[1]);
    _genObj->setFloatVal(LONGFITPAR_2, ff[2]);
    _genObj->setFloatVal(LONGFITPAR_3, ff[3]);
    _genObj->setFloatVal(LONGFITPAR_4, ff[4]);
  }


  int getID    () {return _genObj->getIntVal( ID    );}
  int getNHITS () {return _genObj->getIntVal( NHITS );}
  int getZONE  () {return _genObj->getIntVal( ZONE  );}
  int getNNSEL () {return _genObj->getIntVal( NNSEL );}

  float getENERGY              () {return _genObj->getFloatVal( ENERGY             );}
  float getTHETA               () {return _genObj->getFloatVal( THETA 		    );}
  float getPHI                 () {return _genObj->getFloatVal( PHI 		    );}
  float getTRKDIST_COG         () {return _genObj->getFloatVal( TRKDIST_COG 	    );}
  float getTRKDIST_PROJ        () {return _genObj->getFloatVal( TRKDIST_PROJ 	    );}
  float getTRKDIST_MIN         () {return _genObj->getFloatVal( TRKDIST_MIN 	    );}
  float getTRKDIST_FIRST       () {return _genObj->getFloatVal( TRKDIST_FIRST	    );}
  float getPOINTING_ANGLE      () {return _genObj->getFloatVal( POINTING_ANGLE	    );}
  float getECCENTRICITY        () {return _genObj->getFloatVal( ECCENTRICITY 	    );}
  float getWIDTH               () {return _genObj->getFloatVal( WIDTH 		    );}
  float getVOLUME              () {return _genObj->getFloatVal( VOLUME 	    );}
  float getSTART               () {return _genObj->getFloatVal( START 		    );}
  float getEND                 () {return _genObj->getFloatVal( END 		    );}
  float getMEAN_DEPTH          () {return _genObj->getFloatVal( MEAN_DEPTH	    );}
  float getREL_MEAN_DEPTH      () {return _genObj->getFloatVal( REL_MEAN_DEPTH	    );}
  float getHITEN_MEAN          () {return _genObj->getFloatVal( HITEN_MEAN 	    );}
  float getHITEN_RMS           () {return _genObj->getFloatVal( HITEN_RMS 	    );}
  float getHITEN_Q1            () {return _genObj->getFloatVal( HITEN_Q1 	    );}
  float getHITEN_Q2            () {return _genObj->getFloatVal( HITEN_Q2 	    );}
  float getHITEN_Q3            () {return _genObj->getFloatVal( HITEN_Q3	    );}
  float getNNOUT               () {return _genObj->getFloatVal( NNOUT		    );}
  float getFRACDIM_2           () {return _genObj->getFloatVal( FRACDIM_2 	    );}
  float getFRACDIM_4           () {return _genObj->getFloatVal( FRACDIM_4 	    );}
  float getFRACDIM_8           () {return _genObj->getFloatVal( FRACDIM_8	    );}
  float getTRANSAXISLENGTH_MIN () {return _genObj->getFloatVal( TRANSAXISLENGTH_MIN);}
  float getTRANSAXISLENGTH_MAX () {return _genObj->getFloatVal( TRANSAXISLENGTH_MAX);}
  float getTRACKANGLE_PROJ     () {return _genObj->getFloatVal( TRACKANGLE_PROJ    );}
  float getTRANSVERSERMS_MIN   () {return _genObj->getFloatVal( TRANSVERSERMS_MIN  );}
  float getTRANSVERSERMS_MAX   () {return _genObj->getFloatVal( TRANSVERSERMS_MAX  );}
  float getMOL60_MIN           () {return _genObj->getFloatVal( MOL60_MIN 	    );}
  float getMOL60_MAX           () {return _genObj->getFloatVal( MOL60_MAX	    );}
  float getMOL80_MIN           () {return _genObj->getFloatVal( MOL80_MIN 	    );}
  float getMOL80_MAX           () {return _genObj->getFloatVal( MOL80_MAX	    );}
  float getMOL90_MIN           () {return _genObj->getFloatVal( MOL90_MIN 	    );}
  float getMOL90_MAX           () {return _genObj->getFloatVal( MOL90_MAX	    );}
  float getMOL95_MIN           () {return _genObj->getFloatVal( MOL95_MIN 	    );}
  float getMOL95_MAX           () {return _genObj->getFloatVal( MOL95_MAX          );}

  float getPDF_POINT         () {return _genObj->getFloatVal( PDF_POINT       ); }
  float getPDF_START 	     () {return _genObj->getFloatVal( PDF_START       ); }
  float getPDF_RELMEANLONG   () {return _genObj->getFloatVal( PDF_RELMEANLONG ); }
  float getPDF_LONG2	     () {return _genObj->getFloatVal( PDF_LONG2       ); }
  float getPDF_HITENMEAN     () {return _genObj->getFloatVal( PDF_HITENMEAN   ); }
  float getPDF_HITENWIDTH    () {return _genObj->getFloatVal( PDF_HITENWIDTH  ); }
  float getPDF_FRACDIM	     () {return _genObj->getFloatVal( PDF_FRACDIM     ); }
  float getPDF_WIDTHMAX      () {return _genObj->getFloatVal( PDF_WIDTHMAX    ); }
  float getPDF_WIDTHMIN      () {return _genObj->getFloatVal( PDF_WIDTHMIN    ); }
  float getPDF_CYLIN	     () {return _genObj->getFloatVal( PDF_CYLIN       ); }
  float getPDF_LOWEN	     () {return _genObj->getFloatVal( PDF_LOWEN       ); }
  float getPDF_HIGHEN        () {return _genObj->getFloatVal( PDF_HIGHEN      ); }

  float getPDF_ELECTRON_POINT       () {return _genObj->getFloatVal(   PDF_ELECTRON_POINT      ); }
  float getPDF_ELECTRON_RELMEANLONG () {return _genObj->getFloatVal(   PDF_ELECTRON_RELMEANLONG); }
  float getPDF_ELECTRON_LONG2       () {return _genObj->getFloatVal(   PDF_ELECTRON_LONG2      ); }
  float getPDF_ELECTRON_HITENWIDTH  () {return _genObj->getFloatVal(   PDF_ELECTRON_HITENWIDTH ); }
  float getPDF_ELECTRON_FRACDIM     () {return _genObj->getFloatVal(   PDF_ELECTRON_FRACDIM    ); }
  float getPDF_ELECTRON_WIDTHMAX    () {return _genObj->getFloatVal(   PDF_ELECTRON_WIDTHMAX   ); }
  float getPDF_ELECTRON_WIDTHMIN    () {return _genObj->getFloatVal(   PDF_ELECTRON_WIDTHMIN   ); }
  float getPDF_ELECTRON_CYLIN       () {return _genObj->getFloatVal(   PDF_ELECTRON_CYLIN      ); }
  float getPDF_ELECTRON_EONP        () {return _genObj->getFloatVal(   PDF_ELECTRON_EONP       ); }

  float getELECTRON_EONP     () {return _genObj->getFloatVal( ELECTRON_EONP    ); }
  float getEARLY1DMOL90_MIN  () {return _genObj->getFloatVal( EARLY1DMOL90_MIN ); }
  float getEARLY1DMOL90_MAX  () {return _genObj->getFloatVal( EARLY1DMOL90_MAX ); }
  float getFRACPLAYERS	     () {return _genObj->getFloatVal( FRACPLAYERS      ); }
  float getMAXPLAYERHOLE     () {return _genObj->getFloatVal( MAXPLAYERHOLE    ); }
  float getCLMASS	     () {return _genObj->getFloatVal( CLMASS           ); }
  int getNLAYERS	   () {return _genObj->getIntVal( NLAYERS          ); }
  int getNPLAYERS	   () {return _genObj->getIntVal( NPLAYERS         ); }

  
  float getTUBE_EN(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( TUBE_EN_0 ); break;
    case 1: return _genObj->getFloatVal( TUBE_EN_1 ); break;
    case 2: return _genObj->getFloatVal( TUBE_EN_2 ); break;
    case 3: return _genObj->getFloatVal( TUBE_EN_3 ); break;
    case 4: return _genObj->getFloatVal( TUBE_EN_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getTUBE_N(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( TUBE_N_0 ); break;
    case 1: return _genObj->getFloatVal( TUBE_N_1 ); break;
    case 2: return _genObj->getFloatVal( TUBE_N_2 ); break;
    case 3: return _genObj->getFloatVal( TUBE_N_3 ); break;
    case 4: return _genObj->getFloatVal( TUBE_N_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getLONG_EN(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( LONG_EN_0 ); break;
    case 1: return _genObj->getFloatVal( LONG_EN_1 ); break;
    case 2: return _genObj->getFloatVal( LONG_EN_2 ); break;
    case 3: return _genObj->getFloatVal( LONG_EN_3 ); break;
    case 4: return _genObj->getFloatVal( LONG_EN_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getLONG_N(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( LONG_N_0 ); break;
    case 1: return _genObj->getFloatVal( LONG_N_1 ); break;
    case 2: return _genObj->getFloatVal( LONG_N_2 ); break;
    case 3: return _genObj->getFloatVal( LONG_N_3 ); break;
    case 4: return _genObj->getFloatVal( LONG_N_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getRELLONG_EN(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( RELLONG_EN_0 ); break;
    case 1: return _genObj->getFloatVal( RELLONG_EN_1 ); break;
    case 2: return _genObj->getFloatVal( RELLONG_EN_2 ); break;
    case 3: return _genObj->getFloatVal( RELLONG_EN_3 ); break;
    case 4: return _genObj->getFloatVal( RELLONG_EN_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getRELRELLONG_EN(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( RELRELLONG_EN_0 ); break;
    case 1: return _genObj->getFloatVal( RELRELLONG_EN_1 ); break;
    case 2: return _genObj->getFloatVal( RELRELLONG_EN_2 ); break;
    case 3: return _genObj->getFloatVal( RELRELLONG_EN_3 ); break;
    case 4: return _genObj->getFloatVal( RELRELLONG_EN_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getRELLONG_N(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( RELLONG_N_0 ); break;
    case 1: return _genObj->getFloatVal( RELLONG_N_1 ); break;
    case 2: return _genObj->getFloatVal( RELLONG_N_2 ); break;
    case 3: return _genObj->getFloatVal( RELLONG_N_3 ); break;
    case 4: return _genObj->getFloatVal( RELLONG_N_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }

  float getLONGFITPAR(int i) {
    switch (i) {
    case 0: return _genObj->getFloatVal( LONGFITPAR_0 ); break;
    case 1: return _genObj->getFloatVal( LONGFITPAR_1 ); break;
    case 2: return _genObj->getFloatVal( LONGFITPAR_2 ); break;
    case 3: return _genObj->getFloatVal( LONGFITPAR_3 ); break;
    case 4: return _genObj->getFloatVal( LONGFITPAR_4 ); break;
    default:
      cout << "GarlicClusterVarsGenericObject ERROR invalid index requested: " << i << endl;
      return -999;
    }
  }


private:

  // integer variables
  enum {
    ID=0, NHITS, ZONE, NNSEL,
    NLAYERS,
    NPLAYERS,    
    NINT
  };

  // and the floats
  enum {
    ENERGY=0, THETA, PHI, 
    TRKDIST_COG, TRKDIST_PROJ, TRKDIST_MIN, TRKDIST_FIRST,
    POINTING_ANGLE,
    ECCENTRICITY, WIDTH, VOLUME, 
    START, END, MEAN_DEPTH, REL_MEAN_DEPTH,
    TUBE_EN_0, TUBE_EN_1, TUBE_EN_2, TUBE_EN_3, TUBE_EN_4, 
    TUBE_N_0, TUBE_N_1, TUBE_N_2, TUBE_N_3, TUBE_N_4, 
    LONG_EN_0, LONG_EN_1, LONG_EN_2, LONG_EN_3, LONG_EN_4, 
    LONG_N_0, LONG_N_1, LONG_N_2, LONG_N_3, LONG_N_4, 
    HITEN_MEAN, HITEN_RMS, HITEN_Q1, HITEN_Q3,
    NNOUT,
    FRACDIM_2, FRACDIM_4, FRACDIM_8,
    TRANSAXISLENGTH_MIN, TRANSAXISLENGTH_MAX,
    LONGFITPAR_0, LONGFITPAR_1, LONGFITPAR_2, LONGFITPAR_3, LONGFITPAR_4, 
    TRACKANGLE_PROJ,
    TRANSVERSERMS_MIN, TRANSVERSERMS_MAX,
    MOL60_MIN, MOL60_MAX,
    MOL80_MIN, MOL80_MAX,
    MOL90_MIN, MOL90_MAX,
    MOL95_MIN, MOL95_MAX,
    RELLONG_EN_0, RELLONG_EN_1, RELLONG_EN_2, RELLONG_EN_3, RELLONG_EN_4, 
    RELRELLONG_EN_0, RELRELLONG_EN_1, RELRELLONG_EN_2, RELRELLONG_EN_3, RELRELLONG_EN_4, 
    RELLONG_N_0, RELLONG_N_1, RELLONG_N_2, RELLONG_N_3, RELLONG_N_4, 
    PDF_POINT, 
    PDF_START, 
    PDF_RELMEANLONG, 
    PDF_LONG2, 
    PDF_HITENMEAN, 
    PDF_HITENWIDTH, 
    PDF_FRACDIM, 
    PDF_WIDTHMAX, 
    PDF_WIDTHMIN, 
    PDF_CYLIN,
    PDF_LOWEN, PDF_HIGHEN,
    PDF_ELECTRON_POINT      ,
    PDF_ELECTRON_RELMEANLONG,
    PDF_ELECTRON_LONG2      ,
    PDF_ELECTRON_HITENWIDTH ,
    PDF_ELECTRON_FRACDIM    ,
    PDF_ELECTRON_WIDTHMAX   ,
    PDF_ELECTRON_WIDTHMIN   ,
    PDF_ELECTRON_CYLIN      ,
    PDF_ELECTRON_EONP      ,

    ELECTRON_EONP,
    EARLY1DMOL90_MIN,
    EARLY1DMOL90_MAX,
    FRACPLAYERS,
    MAXPLAYERHOLE,
    CLMASS,
    HITEN_Q2,
    
    NFLOAT
  };

  IMPL::LCGenericObjectImpl* _genObj;


};


#endif
