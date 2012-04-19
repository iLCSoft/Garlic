#include "ECALGarlicTrainer.hh"
#include "ECALGarlicClusterVarsGenericObject.hh"
#include <EVENT/LCCollection.h>
#include <EVENT/LCGenericObject.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/Cluster.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include <iostream>
#include <assert.h>
#include <cmath>

using std::cout;
using std::endl;

ECALGarlicTrainer anECALGarlicTrainer;

ECALGarlicTrainer::ECALGarlicTrainer() : Processor("ECALGarlicTrainer") {
  std::cout << "hello from ECALGarlicTrainer constructor " << this << std::endl;


  // input pre-clusters
  registerInputCollection( LCIO::CLUSTER,
                           "GARLICClusterCollName" ,
                           "Name of GARLIC cluster collection",
                           _GARLICClusterCollectionName,
                           std::string("GARLICClusters") );

  registerInputCollection( LCIO::LCGENERICOBJECT,
                           "GARLICClusterParametersCollName" ,
                           "Name of GARLIC cluster parameters collection",
                           _GARLICClusterParametersCollectionName,
                           std::string("GARLICClusterParameters") );

  registerInputCollection( LCIO::LCRELATION,
                           "GARLICClusterParLinkCollName" ,
                           "Name of GARLIC cluster to parameter relations collection",
                           _GARLICClusterParameterRelationsCollectionName,
                           std::string("GARLICClusterParameterLinks") );

  registerInputCollection( LCIO::LCRELATION,
                           "caloHitRelationsCollName" ,
                           "Name of calohit eelations collection",
                           _calohitRelationsCollectionName,
                           std::string("RelationCaloHit") );

  registerProcessorParameter("outputRootFile",
                             "name of output root file",
                             _outRootFile,
                             std::string("clusterTree.root") );

}


ECALGarlicTrainer::~ECALGarlicTrainer() {
  std::cout << "hello from ECALGarlicTrainer destructor " << this << std::endl;
}

void ECALGarlicTrainer::init()
{

  _fout = new TFile(_outRootFile.c_str(),"recreate");

  _tree = new TTree("ClusterParameters","clusPar");

  _tree->Branch("energy",&_energy, "energy/F");

  _tree->Branch("costheta", &_costheta, "costheta/F");
  _tree->Branch("phi", &_phi, "phi/F");  

  _tree->Branch("nhits", &_nhits, "nhits/I");

  _tree->Branch("enFrac_tube5", &_enFrac_tube5, "enFrac_tube5/F");
  _tree->Branch("enFrac_tube10", &_enFrac_tube10, "enFrac_tube10/F");
  _tree->Branch("enFrac_tube15", &_enFrac_tube15, "enFrac_tube15/F");
  _tree->Branch("enFrac_tube20", &_enFrac_tube20, "enFrac_tube20/F");
  _tree->Branch("enFrac_tube30", &_enFrac_tube30, "enFrac_tube30/F");

  _tree->Branch("nFrac_tube5",  &_nFrac_tube5,  "nFrac_tube5/F");
  _tree->Branch("nFrac_tube10", &_nFrac_tube10, "nFrac_tube10/F");
  _tree->Branch("nFrac_tube15", &_nFrac_tube15, "nFrac_tube15/F");
  _tree->Branch("nFrac_tube20", &_nFrac_tube20, "nFrac_tube20/F");
  _tree->Branch("nFrac_tube30", &_nFrac_tube30, "nFrac_tube30/F");

  _tree->Branch("enFrac_x0_5", &_enFrac_x0_5, "enFrac_x0_5/F");
  _tree->Branch("enFrac_x0_10", &_enFrac_x0_10, "enFrac_x0_10/F");
  _tree->Branch("enFrac_x0_15", &_enFrac_x0_15, "enFrac_x0_15/F");
  _tree->Branch("enFrac_x0_20", &_enFrac_x0_20, "enFrac_x0_20/F");
  _tree->Branch("enFrac_x0_25", &_enFrac_x0_25, "enFrac_x0_25/F");

  _tree->Branch("nFrac_x0_5",  &_nFrac_x0_5,  "nFrac_x0_5/F");
  _tree->Branch("nFrac_x0_10", &_nFrac_x0_10, "nFrac_x0_10/F");
  _tree->Branch("nFrac_x0_15", &_nFrac_x0_15, "nFrac_x0_15/F");
  _tree->Branch("nFrac_x0_20", &_nFrac_x0_20, "nFrac_x0_20/F");
  _tree->Branch("nFrac_x0_25", &_nFrac_x0_25, "nFrac_x0_25/F");


  _tree->Branch("enFracRel_x0_5", &_enFracRel_x0_5, "enFracRel_x0_5/F");
  _tree->Branch("enFracRel_x0_10", &_enFracRel_x0_10, "enFracRel_x0_10/F");
  _tree->Branch("enFracRel_x0_15", &_enFracRel_x0_15, "enFracRel_x0_15/F");
  _tree->Branch("enFracRel_x0_20", &_enFracRel_x0_20, "enFracRel_x0_20/F");
  _tree->Branch("enFracRel_x0_25", &_enFracRel_x0_25, "enFracRel_x0_25/F");

  _tree->Branch("nFracRel_x0_5",  &_nFracRel_x0_5,  "nFracRel_x0_5/F");
  _tree->Branch("nFracRel_x0_10", &_nFracRel_x0_10, "nFracRel_x0_10/F");
  _tree->Branch("nFracRel_x0_15", &_nFracRel_x0_15, "nFracRel_x0_15/F");
  _tree->Branch("nFracRel_x0_20", &_nFracRel_x0_20, "nFracRel_x0_20/F");
  _tree->Branch("nFracRel_x0_25", &_nFracRel_x0_25, "nFracRel_x0_25/F");


  _tree->Branch("start", &_start, "start/F");
  _tree->Branch("end",   &_end,   "end/F");
  _tree->Branch("depth", &_depth, "depth/F");
  _tree->Branch("halfLength", &_halfLength, "halfLength/F");
  _tree->Branch("fullLength", &_fullLength, "fullLength/F");

  _tree->Branch("distToTrackCOG",  &_distToTrackCOG,  "distToTrackCOG/F");
  _tree->Branch("distToTrackPOS",  &_distToTrackPOS,  "distToTrackPOS/F");
  _tree->Branch("dirErr"	,  &_dirErr       ,   "dirErr/F");
  _tree->Branch("angleToTrackPOS",  &_angToTrackPOS, "angleToTrackPOS/F");

  _tree->Branch("Width"		,  &_Width        ,   "Width/F");
  _tree->Branch("Eccentricity"	,  &_Eccentricity ,   "Eccentricity/F");
  _tree->Branch("Volume"        ,  &_Volume       ,   "Volume/F");

  _tree->Branch("hitsMeanEn",&_hitMeanEn, "hitsMeanEn/F");
  _tree->Branch("hitsRMSEn", &_hitRMSEn,  "hitsRMSEn/F");
  _tree->Branch("hitsQ1En",  &_hitQ1En,   "hitsQ1En/F");
  _tree->Branch("hitsQ3En",  &_hitQ3En,   "hitsQ3En/F");

  _tree->Branch("fracDim2",  &_fracDim2,   "fracDim2/F");
  _tree->Branch("fracDim4",  &_fracDim4,   "fracDim4/F");
  _tree->Branch("fracDim8",  &_fracDim8,   "fracDim8/F");

  _tree->Branch("genpdg1",     &_genpdg1,     "genpdg1/I");
  _tree->Branch("genpdg1Frac", &_genpdg1Frac, "genpdg1Frac/F");
  _tree->Branch("genpdg2",     &_genpdg2,     "genpdg2/I");
  _tree->Branch("genpdg2Frac", &_genpdg2Frac, "genpdg2Frac/F");

  _tree->Branch("ecalpdg1",     &_ecalpdg1,     "ecalpdg1/I");
  _tree->Branch("ecalpdg1Frac", &_ecalpdg1Frac, "ecalpdg1Frac/F");
  _tree->Branch("ecalpdg2",     &_ecalpdg2,     "ecalpdg2/I");
  _tree->Branch("ecalpdg2Frac", &_ecalpdg2Frac, "ecalpdg2Frac/F");


  _tree->Branch("genP1",     &_genP1,     "genP1/I");
  _tree->Branch("genP1Frac", &_genP1Frac, "genP1Frac/F");
  _tree->Branch("genP2",     &_genP2,     "genP2/I");
  _tree->Branch("genP2Frac", &_genP2Frac, "genP2Frac/F");

  _tree->Branch("ecalP1",     &_ecalP1,     "ecalP1/I");
  _tree->Branch("ecalP1Frac", &_ecalP1Frac, "ecalP1Frac/F");
  _tree->Branch("ecalP2",     &_ecalP2,     "ecalP2/I");
  _tree->Branch("ecalP2Frac", &_ecalP2Frac, "ecalP2Frac/F");


  _tree->Branch("transProjAx1", &_transAx1, "transProjAx1/F");
  _tree->Branch("transProjAx2", &_transAx2, "transProjAx2/F");

//  _tree->Branch("emFit_a", &_emFit_a   , "emFit_a/F");    
//  _tree->Branch("emFit_b", &_emFit_b   , "emFit_b/F");
//  _tree->Branch("emFit_c", &_emFit_c   , "emFit_c/F");
//  _tree->Branch("emFit_d", &_emFit_d   , "emFit_d/F");
//  _tree->Branch("emFit_xl0", &_emFit_xl0 , "emFit_xl0/F");
//  _tree->Branch("emFit_chi2", &_emFit_chi2, "emFit_chi2/F");

  _tree->Branch("longFit_E", &_longFit_E     , "longFit_E/F");
  _tree->Branch("longFit_a", &_longFit_alpha , "longFit_a/F");
  _tree->Branch("longFit_T", &_longFit_T     , "longFit_T/F");
  _tree->Branch("longFit_off", &_longFit_off     , "longFit_off/F");
  _tree->Branch("longFit_p", &_longFit_prob  , "longFit_p/F");

  _tree->Branch("transRMSmin", &_transRmsMin, "transRMSmin/F");
  _tree->Branch("transRMSmax", &_transRmsMax, "transRMSmax/F");

  _tree->Branch("mol60min", &_mol60min, "mol60min/F");
  _tree->Branch("mol60max", &_mol60max, "mol60max/F");
  _tree->Branch("mol80min", &_mol80min, "mol80min/F");
  _tree->Branch("mol80max", &_mol80max, "mol80max/F");
  _tree->Branch("mol90min", &_mol90min, "mol90min/F");
  _tree->Branch("mol90max", &_mol90max, "mol90max/F");
  _tree->Branch("mol95min", &_mol95min, "mol95min/F");
  _tree->Branch("mol95max", &_mol95max, "mol95max/F");

  _tree->Branch("enFracRel_x0_5", &_enFracRel_x0_5, "enFracRel_x0_5/F");
  _tree->Branch("enFracRel_x0_10", &_enFracRel_x0_10, "enFracRel_x0_10/F");
  _tree->Branch("enFracRel_x0_15", &_enFracRel_x0_15, "enFracRel_x0_15/F");
  _tree->Branch("enFracRel_x0_20", &_enFracRel_x0_20, "enFracRel_x0_20/F");
  _tree->Branch("enFracRel_x0_25", &_enFracRel_x0_25, "enFracRel_x0_25/F");

  _tree->Branch("nFracRel_x0_5",  &_nFracRel_x0_5,  "nFracRel_x0_5/F");
  _tree->Branch("nFracRel_x0_10", &_nFracRel_x0_10, "nFracRel_x0_10/F");
  _tree->Branch("nFracRel_x0_15", &_nFracRel_x0_15, "nFracRel_x0_15/F");
  _tree->Branch("nFracRel_x0_20", &_nFracRel_x0_20, "nFracRel_x0_20/F");
  _tree->Branch("nFracRel_x0_25", &_nFracRel_x0_25, "nFracRel_x0_25/F");

  _tree->Branch("nnOut", &_nnout, "nnOut/F");

}

void ECALGarlicTrainer::processRunHeader(LCRunHeader * run) {
}

void ECALGarlicTrainer::processEvent(LCEvent * evt) {

  LCCollection* GARLICClusters = 0;
  LCCollection* GARLICClusterParLinks = 0;
  LCCollection* GARLICClusterParams = 0;
  LCCollection* calohitrels = 0;
  initTreeVars();

  try {

    GARLICClusters = evt->getCollection(_GARLICClusterCollectionName);
    GARLICClusterParams = evt->getCollection(_GARLICClusterParametersCollectionName);
    GARLICClusterParLinks = evt->getCollection(_GARLICClusterParameterRelationsCollectionName);
    calohitrels = evt->getCollection(_calohitRelationsCollectionName);

    LCRelationNavigator clparnav(GARLICClusterParLinks);
    LCRelationNavigator calhitnav(calohitrels);

    if (GARLICClusterParams->getNumberOfElements()>0) {

      for (int i=0; i<GARLICClusterParams->getNumberOfElements(); i++) {

	LCGenericObject* clparams = dynamic_cast< LCGenericObject* > (GARLICClusterParams->getElementAt(i) );

	ECALGarlicClusterVarsGenericObject vars(clparams);

	_nhits = vars.getNHITS();
	_energy = vars.getENERGY();
	_distToTrackCOG = vars.getTRKDIST_COG();
	_distToTrackPOS = vars.getTRKDIST_PROJ();
	_dirErr = vars.getPOINTING_ANGLE();
	_Eccentricity = vars.getECCENTRICITY();
	_Width = vars.getWIDTH();
	_Volume = vars.getVOLUME();
	_start = vars.getSTART();
	_end = vars.getEND();
	_depth = vars.getMEAN_DEPTH();

	_halfLength = _depth-_start;
	_fullLength = _end-_start;

	_hitMeanEn = vars.getHITEN_MEAN();
	_hitRMSEn = vars.getHITEN_RMS();
	if (_hitRMSEn!=_hitRMSEn) _hitRMSEn=-999; // check for nan...
	_hitQ1En = vars.getHITEN_Q1();
	_hitQ3En = vars.getHITEN_Q3();

	_fracDim2 = vars.getFRACDIM_2();
	_fracDim4 = vars.getFRACDIM_4();
	_fracDim8 = vars.getFRACDIM_8();

	_transAx1 = vars.getTRANSAXISLENGTH_MIN();
	_transAx2 = vars.getTRANSAXISLENGTH_MAX();

	_angToTrackPOS = vars.getTRACKANGLE_PROJ();

	_transRmsMin = vars.getTRANSVERSERMS_MIN();
	_transRmsMax = vars.getTRANSVERSERMS_MAX();

	_mol60min = vars.getMOL60_MIN();
	_mol60max = vars.getMOL60_MAX();
	_mol80min = vars.getMOL80_MIN();
	_mol80max = vars.getMOL80_MAX();
	_mol90min = vars.getMOL90_MIN();
	_mol90max = vars.getMOL90_MAX();
	_mol95min = vars.getMOL95_MIN();
	_mol95max = vars.getMOL95_MAX();

	_enFrac_tube5  = vars.getTUBE_EN(0);
	_enFrac_tube10 = vars.getTUBE_EN(1);
	_enFrac_tube15 = vars.getTUBE_EN(2);
	_enFrac_tube20 = vars.getTUBE_EN(3);
	_enFrac_tube30 = vars.getTUBE_EN(4);

	_nFrac_tube5  = vars.getTUBE_N(0);
	_nFrac_tube10 = vars.getTUBE_N(1);
	_nFrac_tube15 = vars.getTUBE_N(2);
	_nFrac_tube20 = vars.getTUBE_N(3);
	_nFrac_tube30 = vars.getTUBE_N(4);

	_enFrac_x0_5  = vars.getLONG_EN(0);
	_enFrac_x0_10 = vars.getLONG_EN(1);
	_enFrac_x0_15 = vars.getLONG_EN(2);
	_enFrac_x0_20 = vars.getLONG_EN(3);
	_enFrac_x0_25 = vars.getLONG_EN(4);

	_nFrac_x0_5  = vars.getLONG_N(0);
	_nFrac_x0_10 = vars.getLONG_N(1);
	_nFrac_x0_15 = vars.getLONG_N(2);
	_nFrac_x0_20 = vars.getLONG_N(3);
	_nFrac_x0_25 = vars.getLONG_N(4);

	_enFracRel_x0_5  = vars.getRELLONG_EN(0);
	_enFracRel_x0_10 = vars.getRELLONG_EN(1);
	_enFracRel_x0_15 = vars.getRELLONG_EN(2);
	_enFracRel_x0_20 = vars.getRELLONG_EN(3);
	_enFracRel_x0_25 = vars.getRELLONG_EN(4);

	_nFracRel_x0_5  = vars.getRELLONG_N(0);
	_nFracRel_x0_10 = vars.getRELLONG_N(1);
	_nFracRel_x0_15 = vars.getRELLONG_N(2);
	_nFracRel_x0_20 = vars.getRELLONG_N(3);
	_nFracRel_x0_25 = vars.getRELLONG_N(4);

	_nnout = vars.getNNOUT();

	std::map < int, float > ecalEntryContribs;
	std::map < int, float > genContribs;

	std::map < MCParticle*, float > ecalEntryContribsP;
	std::map < MCParticle*, float > genContribsP;

	LCObjectVec cls = clparnav.getRelatedFromObjects( clparams );
	if (cls.size()==1) {
	  Cluster* garClus = dynamic_cast <Cluster*> (cls[0]);


	  float mm(0);
	  for (int kl=0; kl<3; kl++) mm+=pow( garClus->getPosition()[kl], 2);
	  mm=sqrt(mm);

	  _costheta = garClus->getPosition()[2]/mm;
	  _phi = atan2( garClus->getPosition()[1], garClus->getPosition()[0] );

	  CalorimeterHitVec clhits = garClus->getCalorimeterHits();
	  
	  for (size_t ih=0; ih<clhits.size(); ih++) {
	    CalorimeterHit* calhit = clhits[ih];
	    LCObjectVec simhits = calhitnav.getRelatedToObjects(calhit);
	    if (simhits.size()==1) {
	      SimCalorimeterHit* simhit = dynamic_cast <SimCalorimeterHit*> (simhits[0]);

	      for (int ic=0; ic<simhit->getNMCContributions(); ic++) {

		MCParticle* mcp = simhit->getParticleCont(ic);

		std::pair <MCParticle*, MCParticle*> parents = findGeneratedAndECALEnterParent(mcp);

		MCParticle* enterEcalParent = parents.second;
		MCParticle* genParent = parents.first;

		if (enterEcalParent) {
		  if ( ecalEntryContribsP.find( enterEcalParent ) != ecalEntryContribsP.end() )
		    ecalEntryContribsP[enterEcalParent] += calhit->getEnergy();
		  else
		    ecalEntryContribsP[enterEcalParent] = calhit->getEnergy();

		  int pdg = enterEcalParent->getPDG();
		  if ( ecalEntryContribs.find( pdg ) != ecalEntryContribs.end() ) 
		    ecalEntryContribs[pdg] += calhit->getEnergy();
		  else 
		    ecalEntryContribs[pdg] = calhit->getEnergy();
		}

		if (genParent) {

		  if ( genContribsP.find( genParent )!=genContribsP.end() )
		    genContribsP[genParent] += calhit->getEnergy();
		  else
		    genContribsP[genParent] = calhit->getEnergy();

		  int pdg = genParent->getPDG();
		  if ( genContribs.find( pdg )!=genContribs.end() )
		    genContribs[pdg] += calhit->getEnergy();
		  else
		    genContribs[pdg] = calhit->getEnergy();
		}
	      }
	      
	    } else {
	      cout << "warning, found " << simhits.size() << " simcalohits for a caloriemter hit" << endl;
	    }

	  }


	  float maxen(0);
	  float nexten(0);
	  int maxpdg(0);
	  int nextpdg(0);
	  float toten(0);
	  for ( std::map < int, float >::iterator itt=ecalEntryContribs.begin(); itt!=ecalEntryContribs.end(); itt++ ) {
	    //cout << " -- " << itt->first << " " << itt->second << endl;
	    toten+=itt->second;
	    if (itt->second>maxen) {
	      nexten=maxen;
	      nextpdg=maxpdg;
	      maxen=itt->second;
	      maxpdg=itt->first;
	    } else if (itt->second>nexten) {
	      nexten=itt->second;
	      nextpdg=itt->first;
	    }
	  }

	  _ecalpdg1=maxpdg;
	  _ecalpdg1Frac=maxen/toten;
	  _ecalpdg2=nextpdg;
	  _ecalpdg2Frac=nexten/toten;


	  maxen=nexten=0;
	  maxpdg=nextpdg=0;
	  toten=0;
	  for ( std::map < int, float >::iterator itt=genContribs.begin(); itt!=genContribs.end(); itt++ ) {
	    //cout << " -- " << itt->first << " " << itt->second << endl;
	    toten+=itt->second;
	    if (itt->second>maxen) {
	      nexten=maxen;
	      nextpdg=maxpdg;
	      maxen=itt->second;
	      maxpdg=itt->first;
	    } else if (itt->second>nexten) {
	      nexten=itt->second;
	      nextpdg=itt->first;
	    }
	  }

	  _genpdg1=maxpdg;
	  _genpdg1Frac=maxen/toten;
	  _genpdg2=nextpdg;
	  _genpdg2Frac=nexten/toten;


          maxen=nexten=0;
	  MCParticle* maxMCP(0);
	  MCParticle* nextMCP(0);
	  toten=0;
	  for ( std::map < MCParticle*, float >::iterator itt=genContribsP.begin(); itt!=genContribsP.end(); itt++ ) {
	    toten+=itt->second;
	    if (itt->second>maxen) {
	      nexten=maxen;
	      nextMCP=maxMCP;
	      maxen=itt->second;
	      maxMCP=itt->first;
	    } else if (itt->second>nexten) {
	      nexten=itt->second;
	      nextMCP=itt->first;
	    }
	  }
	  if (maxMCP) _genP1=maxMCP->getPDG();
	  _genP1Frac=maxen/toten;
	  if (nextMCP) _genP2=nextMCP->getPDG();
	  _genP2Frac=nexten/toten;

	  maxen=nexten=0;
	  maxMCP=nextMCP=0;
          toten=0;
          for ( std::map < MCParticle*, float >::iterator itt=ecalEntryContribsP.begin(); itt!=ecalEntryContribsP.end(); itt++ ) {
            toten+=itt->second;
            if (itt->second>maxen) {
              nexten=maxen;
              nextMCP=maxMCP;
              maxen=itt->second;
              maxMCP=itt->first;
            } else if (itt->second>nexten) {
              nexten=itt->second;
              nextMCP=itt->first;
            }
          }
	  if (maxMCP) _ecalP1=maxMCP->getPDG();
	  _ecalP1Frac=maxen/toten;
	  if (nextMCP) _ecalP2=nextMCP->getPDG();
	  _ecalP2Frac=nexten/toten;

	} else {
	  cout << "WARNING could not get the corresponding cluster!" << endl;
	}

	_tree->Fill();
	initTreeVars();

      }


    }

  }
  catch(DataNotAvailableException err) {
    streamlog_out ( WARNING )  << "could not find all collections!" << endl;
  };


}


void ECALGarlicTrainer::initTreeVars() {

  _Ns1=0;  
  _Ns2=0;  
  _Ns3=0;  

  _energy=0; 
  _costheta=0;
  _phi=0;
  _EC1=0; 
  _EC4=0; 
  _EC9=0; 
  _start=0; 
  _end=0; 
  _depth=0; 
  _hitDensity=0; 
  _enDensity=0;
  _Es1=0;  
  _Es2=0;  
  _Es3=0;  
  _distToTrackCOG=0;  
  _distToTrackPOS=0;  
  _angToTrackPOS=0;
  _dirErr=0;  
  _Width=0;  
  _Eccentricity=0;  
  _Volume=0;

  _fullLength=0;
  _halfLength=0;
  
  _genpdg1=0; 
  _genpdg2=0; 
  _ecalpdg1=0; 
  _ecalpdg2=0;
  _genpdg1Frac=0;
  _genpdg2Frac=0; 
  _ecalpdg1Frac=0; 
  _ecalpdg2Frac=0;

  _ecalP1=0;
  _ecalP2=0;
  _ecalP1Frac=0;
  _ecalP2Frac=0;

  _genP1=0;
  _genP2=0;
  _genP1Frac=0;
  _genP2Frac=0;

  _hitMeanEn=0;
  _hitRMSEn=0;
  _hitQ1En=0;
  _hitQ3En=0;

  _enFrac_tube5=0;
  _enFrac_tube10=0;
  _enFrac_tube15=0;
  _enFrac_tube20=0;
  _enFrac_tube30=0;

  _fracDim2=0;
  _fracDim4=0;
  _fracDim8=0;

  _transAx1=0;
  _transAx2=0;

//  _emFit_a=0;
//  _emFit_b=0;
//  _emFit_c=0;
//  _emFit_d=0;
//  _emFit_xl0=0;
//  _emFit_chi2=0;

  _longFit_E=0;
  _longFit_alpha=0;
  _longFit_T=0;
  _longFit_off=0;
  _longFit_prob=0;

  _mol80min=0; _mol80max=0;
  _mol90min=0; _mol90max=0;
  _mol95min=0; _mol95max=0;



  return;
}

void ECALGarlicTrainer::check(LCEvent * evt) {
}

void ECALGarlicTrainer::end() {

  _fout->Write();
  _fout->Close();

}



std::pair <MCParticle*, MCParticle*> ECALGarlicTrainer::findGeneratedAndECALEnterParent(MCParticle* mcp) {

  std::pair <MCParticle*, MCParticle*> pars(NULL, NULL);
  std::vector <MCParticle*> generations;

  if (!mcp) return pars;
  
  generations.push_back(mcp);

  MCParticle* mctemp=mcp;
  while ( mctemp->getParents().size()>0 && mctemp->isCreatedInSimulation()) {
    if (mctemp->getParents().size()>1) 
      cout << "warning " << mctemp->getParents().size() << " parents, considering only first one" << endl;
    mctemp = mctemp->getParents()[0];
    generations.push_back(mctemp);
  }

  // find last one not created in simulation
  MCParticle* genpar=NULL;
  for (int i=generations.size()-1; i>=0; i--) {
    if ( !generations[i]->isCreatedInSimulation() ) {
      genpar=generations[i];
    }
  }

  MCParticle* reachECAL=NULL;
  for (int i=generations.size()-1; i>=0; i--) {
    if ( generations[i]->isDecayedInCalorimeter() ) {
      reachECAL=generations[i];
      break;
    }
  }

  pars.first = genpar;
  pars.second = reachECAL;

  return pars;
}
  

