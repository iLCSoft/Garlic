#include "ECALGarlicCluster.hh"

#include <algorithm>
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
#include <assert.h>

#include "TVector3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TGraph.h"

#include <IMPL/CalorimeterHitImpl.h>

//#include "ECALGarlicClusterHelpers.hh"
#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryParameters.hh"
#include "ECALGarlicGeometryHelpers.hh"

#include "ECALGarlicExtendedCluster.hh"
#include "ECALGarlicExtendedHit.hh"

//std::vector <CalorimeterHit*> ECALGarlicCluster::getSeeds(ExtendedCluster2* preClus) {
std::map <CalorimeterHit*, bool> ECALGarlicCluster::getSeeds(ExtendedCluster2* preClus) {

  float energyCutMip = ECALGarlicAlgorithmParameters::Instance().GetSeedHitEnergyCut();
  float seedThresholdMIP = ECALGarlicAlgorithmParameters::Instance().GetSeedEnergyCut();
  float cutoffdist = 
    ECALGarlicAlgorithmParameters::Instance().GetSeedDistanceCut()*
    ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  if ( cutoffdist>ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius() ) {
    cout << "WARNING, requested scale for seeding is larger than Moliere radius...setting to Moliere radius: " << 
      " distance (cells)= " << ECALGarlicAlgorithmParameters::Instance().GetSeedDistanceCut() << 
      " cell size = " << ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0] << 
      " distance (mm) = " << ECALGarlicAlgorithmParameters::Instance().GetSeedDistanceCut()*
      ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0] << 
      " Moliere: " << ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius() << endl;
    cutoffdist=ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();
  }

  int maxpseulayer = ECALGarlicAlgorithmParameters::Instance().GetSeedNLayers();

  //  const int maxseed=300;

  //  std::vector <CalorimeterHit*> seeds;
  std::map <CalorimeterHit*, bool> seeds;

  if (_fhistos) _fhistos->cd();


  std::pair < TH2F*, TH2I* > histos = preClus->GetProjectionHistos( energyCutMip, maxpseulayer );
  TH2F* henergy = histos.first;
  TH2I* hhits = histos.second;

  if (!henergy) {
    return seeds;
  }

  TString blah;

  if (_hnn!="blah") blah="_"+_hnn+"_";
  else blah="_n_";
  blah+="nhits";
  blah+=preClus->getHits()->size();
  blah+="_";
  blah+=_nSaveHist+1;
    
  // need to clone the histos and rebin them to cell size (originally they were 1mm)
  float cellsize = ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  int ics = int(cellsize);

  henergy->Rebin2D(ics, ics);
  hhits->Rebin2D(ics, ics);

  std::vector < std::pair < int, int > > seedCells;

  seedCells.clear();

  // find highest energy bin
  int maxbin = henergy->GetMaximumBin();
  float maxen = henergy->GetBinContent( maxbin );

  //  int iseed=0;
  //float seedXvec[maxseed]={0};
  //float seedYvec[maxseed]={0};
  std::vector <float> seedXvec;
  std::vector <float> seedYvec;

  //  int ibadseed=0;
  //float badseedXvec[maxseed]={0};
  //float badseedYvec[maxseed]={0};
  std::vector <float> badseedXvec;
  std::vector <float> badseedYvec;
  
  //  while (maxen>seedThresholdMIP ) {
  while (maxen>=energyCutMip) {

    // look at 3x3 neighbouring bins, get COG position
    int binx(0); int biny(0);
    int dummy;
    henergy->GetBinXYZ(maxbin, binx, biny, dummy);

    float cogx(0);
    float cogy(0);
    float toten(0);

    for (int ix=-1; ix<=1; ix++) {
      for (int iy=-1; iy<=1; iy++) {
	float en = henergy->GetBinContent(binx+ix, biny+iy);
	cogx+= en*henergy->GetXaxis()->GetBinCenter(binx+ix);
	cogy+= en*henergy->GetYaxis()->GetBinCenter(biny+iy);
	toten+=en;
      }
    }
    cogx/=toten;
    cogy/=toten;

    // look for cells within some distance of this position
    int windowsize = int ( 1.5 * cutoffdist / henergy->GetXaxis()->GetBinWidth(1) );

    float seed_en(0);
    float seed_rms(0);
    float seed_rmsX(0);
    float seed_rmsY(0);
    float seed_x(0);
    float seed_y(0);
    int seed_nhits(0);
    
    for (int ix=-windowsize; ix<=windowsize; ix++) {
      for (int iy=-windowsize; iy<=windowsize; iy++) {

	float x = henergy->GetXaxis()->GetBinCenter(binx+ix);
	float y = henergy->GetYaxis()->GetBinCenter(biny+iy);

	float en = henergy->GetBinContent(binx+ix, biny+iy);

	float dx = cogx - x;
	float dy = cogy - y;

	float dist = sqrt( pow(dx, 2) + pow (dy, 2) );
	if (dist>cutoffdist) continue;

	seed_en+=en;
	seed_nhits+= (int) hhits->GetBinContent(binx+ix, biny+iy);

	seed_x+=x*en;
	seed_y+=y*en;

	seed_rms+=pow(en*dist, 2);
	seed_rmsX+=pow(en*dx, 2);
	seed_rmsY+=pow(en*dy, 2);

	seedCells.push_back ( std::pair <int, int> (binx+ix, biny+iy) );

      }
    }
    seed_rms=sqrt(seed_rms);
    seed_rms/=seed_en;

    seed_rmsX=sqrt(seed_rmsX);
    seed_rmsX/=seed_en;

    seed_rmsY=sqrt(seed_rmsY);
    seed_rmsY/=seed_en;

    seed_x/=seed_en;
    seed_y/=seed_en;

    bool goodSeed = seed_nhits>=ECALGarlicAlgorithmParameters::Instance().GetSeedMinHits() && toten>=seedThresholdMIP ;

    float pos[3];
    preClus->GetGlobalPositionFromLocal(seed_x, seed_y, pos);

    CalorimeterHitImpl* a_seed=new CalorimeterHitImpl();
    a_seed->setEnergy(seed_en);
    a_seed->setPosition(pos);
    seeds[a_seed] = goodSeed;

    if (goodSeed) {
      //if (iseed<maxseed) {
      //seedXvec[iseed]=seed_x;
      //seedYvec[iseed]=seed_y;
      seedXvec.push_back(seed_x);
      seedYvec.push_back(seed_y);
      //      iseed++;
      // } else {
      //	cout << "WARNING, too many seeds, only keeping first " << maxseed << endl;
      //	cout << " increase maxseed in ECALGarlicCluster::getSeeds" << endl;
      // }
    } else {
      //      if (ibadseed<maxseed) {
      //badseedXvec[ibadseed]=seed_x;
      //badseedYvec[ibadseed]=seed_y;
      //ibadseed++;
      badseedXvec.push_back(seed_x);
      badseedYvec.push_back(seed_y);
      //      }
    }
    
    for (size_t i=0; i<seedCells.size(); i++) {
      henergy->SetBinContent(seedCells[i].first, seedCells[i].second, 0);
      hhits  ->SetBinContent(seedCells[i].first, seedCells[i].second, 0);
    }
    seedCells.clear();
    maxbin = henergy->GetMaximumBin();
    maxen = henergy->GetBinContent( maxbin ); 
  }


  if (henergy && hhits && _fhistos) {

    _fhistos->cd();
    //    _fhistos->ls();

    float x[100];
    float y[100];
    int np=std::min(100, int(seedXvec.size()));
    for (int i=0; i<np; i++) {
      x[i]=seedXvec[i];
      y[i]=seedYvec[i];
    }

    TGraph* gr = new TGraph(np, x, y);
    gr->SetNameTitle("seeds"+blah,"seeds"+blah);
    gr->Write();

    np=std::min(100, int(badseedXvec.size()));
    for (int i=0; i<np; i++) {
      x[i]=badseedXvec[i];
      y[i]=badseedYvec[i];
    }

    TGraph* gr2 = new TGraph(np, x, y);
    gr2->SetNameTitle("badseeds"+blah,"badseeds"+blah);
    gr2->Write();
    _nSaveHist++;
  }

  if (!_fhistos) {
    if (henergy) delete henergy;
    if (hhits) delete hhits;
  }

  return seeds;
}

vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > ECALGarlicCluster::getElectrons( ExtendedCluster2* preCluster, vector <ExtendedTrack* > trks ) {

  if (_verbose) 
    cout << "hello from ECALGarlicCluster::getElectrons "<< trks.size() << " " << preCluster->getNhits() << endl;
  
  vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > electrons;

  float addCut = 2.4*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // be a bit looser...

  //  float maxDistCut = ECALGarlicAlgorithmParameters::Instance().GetClusterMaxDist()*ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();
  float maxDistCut = ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();


  //  float maxAddDist = ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist()*ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();
  float maxAddDist = ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();
  //ECALGarlicAlgorithmParameters::Instance().Get_MaxTransMergeDist();
  
  for ( size_t it=0; it<trks.size(); it++ ) {

    if (_verbose) 
      cout << "looking at track " << it << endl;

    // position and direction of track at ECAL front face
    const float* position  = trks[it]->getEcalEntryPos();
    const float* direction = trks[it]->getEcalEntryDir();

    // build a "seed"
    CalorimeterHitImpl seed;
    seed.setEnergy(1);
    seed.setPosition(position);

    // make the core
    ExtendedCluster2* electronCore = BuildCore( &seed, preCluster, direction);

    electronCore->setAssociatedTrack(trks[it]);

    if (_verbose) 
      cout << "core hits: " << electronCore->getHits()->size() << endl;

    if ( electronCore->getHits()->size() ==0 ) {
      delete electronCore;
      continue;
    }

    if (_verbose) {
      cout << "electron core selected ? " << electronCore->getElectronCutSel() << " : ";
      cout << electronCore->getElectronCutSel(ExtendedCluster2::TRANS ) << " ";
      cout << electronCore->getElectronCutSel(ExtendedCluster2::LONG  ) << " ";
      cout << electronCore->getElectronCutSel(ExtendedCluster2::HITEN ) << " ";
      cout << electronCore->getElectronCutSel(ExtendedCluster2::POINT ) << " ";
      cout << endl;
    }

    //    photon_longProfileTight

    // find hits not too far from core axis
    float point2[3]; // second point along shower axis
    for (int i=0; i<3; i++)
      point2[i]=electronCore->getCentreOfGravity()[i]-direction[i];

    std::vector < ExtendedHit2* > possibleHits;
    for ( size_t ij=0; ij<preCluster->getHits()->size(); ij++) {
      ExtendedHit2* hh = preCluster->getHits()->at(ij);
      if ( find(possibleHits.begin(), possibleHits.end(), hh)==possibleHits.end()  && 
	   find(electronCore->getHits()->begin(), electronCore->getHits()->end(), hh)==electronCore->getHits()->end() 
	   ) {
	// this assumed projective to IP...
	//	float dist = ECALGarlicGeometryHelpers::GetDistToLine(hh->getCaloHit()->getPosition(), origin, electronCore->getCentreOfGravity()); 
	float dist = ECALGarlicGeometryHelpers::GetDistToLine(hh->getCaloHit()->getPosition(), point2, electronCore->getCentreOfGravity());
	if (dist<maxAddDist) possibleHits.push_back(hh);
      }
    }

    if (_verbose) 
      cout << " got " << possibleHits.size() << " possible hits within a distance " << maxAddDist << " from axis " << endl;



    // at each iteration, we should check that we don;t screw up electron properties.

    int addedThisRound=999;
    int nrounds(0);
    while ( addedThisRound>0) {
      addedThisRound=0;
      std::vector < ExtendedHit2* > hitsToAdd;
      for ( int ij=possibleHits.size()-1; ij>=0; ij--) {
	float distToCore = electronCore->getDistToClusterAxis( possibleHits[ij]->getPosition(), 1);
	distToCore = std::min( distToCore, electronCore->getDistToClusterAxis( possibleHits[ij]->getPosition(), 2) );
	if ( distToCore > maxDistCut ) continue; // too far from core
	std::vector < ExtendedHit2* > * corehits = electronCore->getHits();
	for (size_t ik=0; ik<corehits->size(); ik++) {
	  ExtendedHit2* chit = corehits->at(ik);
	  float dist(0);
	  for (int i=0; i<3; i++) 
	    dist += pow ( possibleHits[ij]->getPosition()[i] - chit->getCaloHit()->getPosition()[i], 2);
	  dist = sqrt(dist);
	  if (dist<addCut) {
	    addedThisRound++;
	    //	    hitsToAdd.push_back(chit);
	    hitsToAdd.push_back(possibleHits[ij]);
	    possibleHits.erase( possibleHits.begin()+ij );
	    break;
	  }
	}
      } // possible hits loop
      if (_verbose) 
	cout << "adding to electron core: " << hitsToAdd.size() << " from " << possibleHits.size() << endl;
      if (hitsToAdd.size()>0) {
	electronCore->addHits(hitsToAdd);
      }
    } // while loop





    float EonP = electronCore->getEnergy()/electronCore->getAssociatedTrack()->getTotalMomentum();

    bool lowenergy = electronCore->getEnergy()<1.;

    float PDFcut_loose = lowenergy ? 
      ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_lowE() : 
      ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_hiE() ;
    float PDFcut_tight = lowenergy ? 
      ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_lowE() : 
      ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_hiE() ;

    if (_verbose) {
      cout << " -- get_pdf_point       " << electronCore->get_pdf_point      () << endl;
      cout << " -- get_pdf_start       " << electronCore->get_pdf_start      () << endl;
      cout << " -- get_pdf_relmeanlong " << electronCore->get_pdf_relmeanlong() << endl;
      cout << " -- get_pdf_long2       " << electronCore->get_pdf_long2      () << endl;
      //    cout << " -- get_pdf_hitenmean   " << electronCore->get_pdf_hitenmean  () << endl;
      cout << " -- get_pdf_hitenwidth  " << electronCore->get_pdf_hitenwidth () << endl;
      cout << " -- get_pdf_fracdim     " << electronCore->get_pdf_fracdim    () << endl;
      cout << " -- get_pdf_widthMax    " << electronCore->get_pdf_widthMax   () << endl;
      cout << " -- get_pdf_widthMin    " << electronCore->get_pdf_widthMin   () << endl;
      cout << " -- get_pdf_cylin       " << electronCore->get_pdf_cylin      () << endl;
      
      cout << "  == get_electronpdf_point       " <<  electronCore->get_electronpdf_point      () << endl;
      cout << "  == get_electronpdf_relmeanlong " <<  electronCore->get_electronpdf_relmeanlong() << endl;
      cout << "  == get_electronpdf_long2       " <<  electronCore->get_electronpdf_long2      () << endl;
      cout << "  == get_electronpdf_hitenwidth  " <<  electronCore->get_electronpdf_hitenwidth () << endl;
      cout << "  == get_electronpdf_fracdim     " <<  electronCore->get_electronpdf_fracdim    () << endl;
      cout << "  == get_electronpdf_widthMax    " <<  electronCore->get_electronpdf_widthMax   () << endl;
      cout << "  == get_electronpdf_widthMin    " <<  electronCore->get_electronpdf_widthMin   () << endl;
      cout << "  == get_electronpdf_cylin       " <<  electronCore->get_electronpdf_cylin      () << endl;
      cout << "  == get_electronpdf_eonp        " <<  electronCore->get_electronpdf_eonp       () << endl;
      cout << "  == get_total_electronpdf       " << 
	electronCore->get_total_electronpdf(true) << " " << electronCore->get_total_electronpdf(false) << endl;
      
      //
      cout << " widths: " << electronCore->getTransverseRMS().first << " " << electronCore->getTransverseRMS().second << endl;
      cout << " fracdim: " << electronCore->getFractalDimension()[0] << endl;
      
      cout << "electron candidate? " << " energy " << electronCore->getEnergy();
      cout << " ; eonp = " << EonP << endl;
      
      cout << "  == get_total_electronpdf " << electronCore->get_total_electronpdf( false ) << " " << 
	electronCore->get_total_electronpdf( true ) << " pdf cuts = " << PDFcut_loose << " " << PDFcut_tight << endl;


      cout << "electron candidate: " << electronCore->getNhits() << " " << electronCore->getStart() << " " << electronCore->get_total_electronpdf( lowenergy ) << " " << EonP << endl;
    }
    // make cuts to find electron clusters
    if (electronCore->getNhits()>0  && electronCore->getStart()<2. 
	//	&& electronCore->get_total_electronpdf( lowenergy ) > PDFcut_loose 
	&& EonP > 0.4 && EonP < 1.4
	)  {
      if (_verbose) 
	cout << "accepting electron with e/p = " << EonP << endl;

      // if (EonP<0.25 || EonP>2) cout << "WARNING selectong electron with crazy E/p = " << EonP << endl;

      electrons.push_back( std::pair < ExtendedTrack*, ExtendedCluster2* > ( trks[it], electronCore ) );
      trks[it]->setElectronSel(1);
    } else {
      // delete the new cluster
      delete electronCore; electronCore=NULL;
    }
  } // loop over tracks

  return electrons;
}




std::map < CalorimeterHit*, ExtendedCluster2* > ECALGarlicCluster::getClusters( ExtendedCluster2* preCluster, std::map < CalorimeterHit*, ExtendedCluster2* > cores ) {
  // assign unassigned hits to "best" core

  //  float addCut = 1.9*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  float addCut = 2.4*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  //  int nIterations = ECALGarlicAlgorithmParameters::Instance().GetClusterNIterations();
  float maxDistCut = ECALGarlicAlgorithmParameters::Instance().GetClusterMaxDist()*ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();

  // duplicate the cores into the clusters
  std::map < CalorimeterHit*, ExtendedCluster2* > clusters;
  for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator clit = cores.begin(); clit!=cores.end(); clit++ ) {
    ExtendedCluster2* cluster = new ExtendedCluster2();
    (*cluster) = *(clit->second);
    clusters[clit->first] = cluster;
  }

  // first extract the unassigned hits within precluster
  std::vector < ExtendedHit2* > unassignedHits = *(preCluster->getHits());
  for ( int ihh=unassignedHits.size()-1; ihh>=0; ihh-- ) {
    for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator clit = clusters.begin(); clit!=clusters.end(); clit++ ) {
      if ( find( clit->second->getHits()->begin(), clit->second->getHits()->end(), unassignedHits[ihh] ) != clit->second->getHits()->end() ) {
	// cout << "erasing hit " << ihh << " at " << unassignedHits[ihh] << endl;
	unassignedHits.erase( unassignedHits.begin()+ihh );
	break;
      }
    }
  }

  //  for (int iiter=0; iiter<nIterations; iiter++) {
  int addedthisround=999;
  while ( addedthisround>0 ) {
    addedthisround=0;
    std::map < CalorimeterHit*, std::vector < ExtendedHit2* > > preassignedHits;
    // find distance to closest hit in core
    for ( int ihh=unassignedHits.size()-1; ihh>=0; ihh-- ) {
      float mindist=99999;
      CalorimeterHit* closestCore=0;
      for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator clit = clusters.begin(); clit!=clusters.end(); clit++ ) {

	// distance from hit to core axis
	float distToCore = clit->second->getDistToClusterAxis( unassignedHits[ihh]->getPosition(), 0 ); // this assumes pointing photon
	if ( distToCore>maxDistCut ) continue; // do not consider this combination

	std::vector < ExtendedHit2* >* corehits = clit->second->getHits();
	for (size_t ich=0; ich<corehits->size(); ich++) {
	  float dist(0);
	  for (int i=0; i<3; i++) dist+=pow( (*corehits)[ich]->getPosition()[i] - unassignedHits[ihh]->getPosition()[i] , 2. );
	  dist = sqrt(dist);
	  if (dist<mindist) {
	    mindist=dist;
	    closestCore=clit->first;
	  }
	}
      }
      if (mindist<addCut) {
	addedthisround++;
	if ( preassignedHits.find(closestCore)!=preassignedHits.end() ) {
	  preassignedHits[closestCore].push_back(unassignedHits[ihh]);
	} else {
	  std::vector < ExtendedHit2* > temp;
	  temp.push_back( unassignedHits[ihh] );
	  preassignedHits[closestCore]=temp;
	}
	unassignedHits.erase(unassignedHits.begin()+ihh);
      }
    }
    for ( std::map < CalorimeterHit*, std::vector < ExtendedHit2* > >::iterator preass = preassignedHits.begin(); preass!=preassignedHits.end(); preass++) {
      //      cout << "adding preassigned hits to core" << endl;
      clusters[preass->first]->addHits(preass->second);
    }
  } // iterations

  return clusters;
}


std::map < CalorimeterHit*, ExtendedCluster2* > ECALGarlicCluster::getCores( ExtendedCluster2* preCluster, std::vector <CalorimeterHit*> seeds ) {
  std::map < CalorimeterHit*, ExtendedCluster2* > cores;
  if (seeds.size()==0) {
    cout << "WARNING getCores: no seeds given!" << endl;
    return cores;
  }

  for (size_t iseed = 0; iseed<seeds.size(); iseed++) {
    CalorimeterHit* seed = seeds[iseed];
    const float* clusterDirection = seed->getPosition(); // look for pointing photons
    cores[seed] = BuildCore( seed, preCluster, clusterDirection);
    //    cout << " core " << iseed << " has " << cores[seed]->getHits()->size() << " hits" << endl;
  }

  // check for shared hits...
  // give to most energetic
  std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt;
  std::map < CalorimeterHit*, ExtendedCluster2* >::iterator jtt;

  for (itt=cores.begin(); itt!=cores.end(); itt++) {
    std::vector<ExtendedHit2*>* ih = itt->second->getHits();
    if (ih->size()==0) continue;
    for (jtt=itt; jtt!=cores.end(); jtt++) {
      if (jtt==itt) continue;
      if (jtt==cores.end()) continue;
      for (int ii= (int) ih->size()-1; ii>=0; ii--) {
	std::vector<ExtendedHit2*>* jh = jtt->second->getHits();
	if (jh->size()==0) continue;
	for (int jj= (int) jh->size()-1; jj>=0; jj--) {
	  if ( (*ih)[ii] == (*jh)[jj] ) {
	    float e_i = itt->second->getHits()->size()==0 ? 0 : itt->second->getEnergy();
	    float e_j = jtt->second->getHits()->size()==0 ? 0 : jtt->second->getEnergy();
	    if ( e_i > e_j ) {
	      jtt->second->removeHit( (*ih)[ii] );
	    } else {
	      itt->second->removeHit( (*ih)[ii] );
	    }
	  }
	}
      }
    }
  }

  // remove any cores which have been completely eaten by neighbours (shared hits)...
  bool changed=true;
  while (changed) {
    changed=false;
    for (itt=cores.begin(); itt!=cores.end(); itt++) {
      if (itt->second->getHits()->size()==0) {
	cores.erase(itt);
	changed=true;
	break;
      }
    }
  }

  return cores;
}


ExtendedCluster2* ECALGarlicCluster::BuildCore(CalorimeterHit* mySeed, ExtendedCluster2* preCluster, const float *clusterDir) {

  int nl_1         = ECALGarlicAlgorithmParameters::Instance().GetCoreLayersSection1();
  int allowedGap_1 = ECALGarlicAlgorithmParameters::Instance().GetCoreMaxHoleSection1();
  int allowedGap_2 = ECALGarlicAlgorithmParameters::Instance().GetCoreMaxHoleSection2();

  ExtendedHit2 seedHit(mySeed);

  std::vector <ExtendedHit2*> coreHits;

  std::vector <ExtendedHit2*> allPreclusterHits = *(preCluster->getHits());

  std::map <int, int> hitPseudoLayersBarrel;
  std::map <int, int> hitPseudoLayersEndcap;

  // define the distance window to add hits to the core
  //  float distCut = 1.5*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  float distCut = 
    ECALGarlicAlgorithmParameters::Instance().GetCoreDistanceCut()*
    ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

  //  cout << "cell size, dist cut = " << ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0] << " " << distCut << endl;

  for (size_t ih=0; ih<allPreclusterHits.size(); ih++) {

    ExtendedHit2* thisHit = allPreclusterHits[ih];

    // project all hit onto plane defined by seed position and clusterDir
    const float* hitpos = thisHit->getPosition();
    float proj[3];
    ECALGarlicGeometryHelpers::GetGeneralPointPlaneProjection(hitpos, mySeed->getPosition(), mySeed->getPosition(), clusterDir, proj);

    // distance from seed position to projected hit
    float dist(0);
    for (int i=0; i<3; i++) dist+= pow( proj[i] - mySeed->getPosition()[i], 2);
    dist=sqrt(dist);

    // if (dist>1000) {
    //   cout << "ECALGarlicCluster::BuildCore : probable ERROR in projection measurement! too large distance " << dist << endl;
    //   cout << "seed pos = " << mySeed->getPosition()[0] << " " << mySeed->getPosition()[1] << " " << mySeed->getPosition()[2] << endl;
    //   cout << "hit pos = " << hitpos[0] << " " << hitpos[1] << " " << hitpos[2] << endl;
    //   cout << "projection = " << proj[0] << " " << proj[1] << " " << proj[2] << endl;
    //   cout << "cluster dir = " << clusterDir[0] << " " << clusterDir[1] << " " << clusterDir[2] << endl;
    // }

    // take account of non-normal detector cells
    float inflationFactor(1.);
    float mag(0);
    for (int i=0; i<3; i++) mag+=pow( thisHit->getPosition()[i], 2. );
    mag = sqrt(mag);
    float r(0);
    for (int i=0; i<2; i++) r+=pow( thisHit->getPosition()[i], 2. );
    r=sqrt(r);
    float theta = acos ( fabs( thisHit->getPosition()[2])/mag );
    if ( thisHit->getZone()==CALHITZONE_BARREL ) {
      float alpha = acos(-1)/2. - theta;
      inflationFactor = 1./cos(alpha);
    } else if ( thisHit->getZone()==CALHITZONE_ENDCAP || thisHit->getZone()==CALHITZONE_RING) {
      inflationFactor = 1./cos(theta);
    } else {
      cout << "hit " << ih << " in strange detector region " << thisHit->getZone() << " position " <<
	thisHit->getPosition()[0] << " " << thisHit->getPosition()[1] << " " << thisHit->getPosition()[2] << endl;
    }

    float thisHitDistCut = distCut*inflationFactor;

    //    cout << thisHit->getZone() << " " << mag << " " << r << " " << theta << " " << inflationFactor << " " << distCut << " " << thisHitDistCut << endl;

    if ( seedHit.getZone()==CALHITZONE_BARREL && thisHit->getZone()==CALHITZONE_ENDCAP ) 
      thisHitDistCut*=3;

    if (dist<thisHitDistCut) {
      coreHits.push_back(thisHit);
      int player = thisHit->getPseudoLayer();
      if ( thisHit->getZone()==CALHITZONE_BARREL ) {
	if ( hitPseudoLayersBarrel.find(player)==hitPseudoLayersBarrel.end() ) hitPseudoLayersBarrel[player]=1;
	else hitPseudoLayersBarrel[player]++;
      } else if ( thisHit->getZone()==CALHITZONE_ENDCAP || thisHit->getZone()==CALHITZONE_RING) {
	if ( hitPseudoLayersEndcap.find(player)==hitPseudoLayersEndcap.end() ) hitPseudoLayersEndcap[player]=1;
	else hitPseudoLayersEndcap[player]++;
      }

      //      cout << "added hit to core: player, dist , cut " << player << " " << dist << " " << distCut << " " << thisHitDistCut << endl;

    }
  }

  // now check for missed pseudo-layers
  if ( (hitPseudoLayersBarrel.size()==0 && hitPseudoLayersEndcap.size()>0) ||
       (hitPseudoLayersBarrel.size()>0 && hitPseudoLayersEndcap.size()==0) ) { // all barrel or all endcap
    std::map <int, int>* pls = hitPseudoLayersBarrel.size() ? &hitPseudoLayersBarrel : &hitPseudoLayersEndcap;

    int firstHit(-1);
    int consecutiveEmpty(0);
    int lastEmpty(-1);
    int lastFilled(-1);

    bool stopit(false);

    for (int i=0; i< ECALGarlicGeometryParameters::Instance().Get_nPseudoLayers(); i++) {
      if ( pls->find(i) != pls->end() ) { // this player hit
	if (firstHit<0) firstHit=i;
	lastFilled=i;
	lastEmpty=-1;
	consecutiveEmpty=0;
      } else { // no hit
	if (firstHit>=0) {
	  consecutiveEmpty++;
	  lastEmpty=i;
	}
      }
      // cout << "player " << i << " fHit=" << firstHit << " lFill=" << lastFilled << " lEmp=" << lastEmpty << " consEm=" << consecutiveEmpty << endl;
      if ( lastFilled<=nl_1 ) {
	if ( consecutiveEmpty > allowedGap_1 ) stopit = true;
      } else {
	if ( consecutiveEmpty > allowedGap_2 ) stopit = true;
      }
      if (stopit) {
	// cout << "deciding to stop this core building, gapSize=" << consecutiveEmpty << " lastFilledLayer=" << lastFilled << endl;
	break;
      }
    }

    // remove the hits beyond the break point
    for (int ihh=coreHits.size()-1; ihh>=0; ihh--) {
      if (coreHits[ihh]->getPseudoLayer()>lastFilled) coreHits.erase(coreHits.begin()+ihh);
    }
    
  } else if (hitPseudoLayersBarrel.size()==0 && hitPseudoLayersEndcap.size()==0) {

    // cout << "STRANGE (probably an error), found no hits in core!" << endl;

  } else { // got mixed barrel endcap case

    cout << "mixed barrel endcap case!  nBarrel=" << hitPseudoLayersBarrel.size() << " nEnd=" << hitPseudoLayersEndcap.size() << endl;
    cout << "   WARNING do not yet know how to deal with this!!" << endl;
    cout << "    will not break the core" << endl;

    // cout << "barrel hits " << endl;
    // for ( std::map <int, int>::iterator itt = hitPseudoLayersBarrel.begin(); itt!=hitPseudoLayersBarrel.end(); itt++) {
    //   cout << "  " << itt->first << " " << itt->second << endl;
    // }
    // cout << "endcap hits " << endl;
    // for ( std::map <int, int>::iterator itt = hitPseudoLayersEndcap.begin(); itt!=hitPseudoLayersEndcap.end(); itt++) {
    //   cout << "  " << itt->first << " " << itt->second << endl;
    // }



  }

  ExtendedCluster2* coreClus = new ExtendedCluster2();
  coreClus->addHits( coreHits );
  coreClus->setReferenceDirection( clusterDir );

  return coreClus;  
}

void ECALGarlicCluster::mergeUnclusteredAndElectronsPhotons_byPDF(std::vector < ExtendedHit2* > & unclusteredHits, 
								  vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > & allElectrons,
								  std::map < CalorimeterHit*, ExtendedCluster2* > & clusters ) {

  //  cout << " hello from mergeUnclusteredAndElectronsPhotons_byPDF " << unclusteredHits.size() << endl;


  std::vector < ExtendedCluster2* > allClusters;
  for ( size_t i=0; i<allElectrons.size(); i++) {
    allClusters.push_back(allElectrons[i].second);
  }
  for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt=clusters.begin(); itt!=clusters.end(); itt++) {
    allClusters.push_back(itt->second);
  }

  std::vector < std::vector < ExtendedHit2* > > hitsToAdd;
  for (size_t ie=0; ie<allClusters.size(); ie++ ) {
    hitsToAdd.push_back(  std::vector < ExtendedHit2* > () );
  }

  const int nsteps=2; // in first step get the closer ones
  for (int iss=0; iss<nsteps; iss++) {
    float maxJoinDist = ((iss+1.0)/nsteps)*
      ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist()*
      ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius();

    //ECALGarlicAlgorithmParameters::Instance().Get_MaxTransMergeDist ();

    for (size_t ih=0; ih<unclusteredHits.size(); ih++) {
      ExtendedHit2* uhit = unclusteredHits[ih];
      float mindist = 9999;
      int imindist=-1;
      for (size_t ie=0; ie<allClusters.size(); ie++ ) {
	ExtendedCluster2* ele = allClusters[ie];
	if ( find ( ele->getHits()->begin(), ele->getHits()->end(), uhit )!=ele->getHits()->end() ) {
	  cout << "WARNING: this hit is already in a cluster! " << ie << endl;
	}
	float dist = ele->getDistToClusterAxis(uhit->getCaloHit()->getPosition(), false);
	float distIP = ele->getDistToClusterAxis(uhit->getCaloHit()->getPosition(), true);
	float distCOG(0);
	for (int i=0; i<3; i++) distCOG+=pow( uhit->getCaloHit()->getPosition()[i] - ele->getCentreOfGravity()[i] , 2 );
	distCOG=sqrt(distCOG);
	float mm = std::min( std::min( dist, distIP ), distCOG );
	if (mm<mindist) {
	  mindist=mm;
	  imindist=ie;
	}
      }
      if ( mindist<maxJoinDist ) {
	hitsToAdd[imindist].push_back( uhit );
      }
    }
  
    for (size_t ie=0; ie<allClusters.size(); ie++ ) {

      if ( hitsToAdd[ie].size()==0 ) continue;

      // does adding the hits make it better?

      bool lowEn = allClusters[ie]->getEnergy()<1.0;
      bool isElec = allClusters[ie]->getAssociatedTrack();

      float oldPDF = isElec ? allClusters[ie]->get_total_electronpdf( lowEn ) : allClusters[ie]->get_total_pdf( lowEn ) ;

      ExtendedCluster2* combo = new ExtendedCluster2();

      combo->addHits( *(allClusters[ie]->getHits()) );

      combo->addHits( hitsToAdd[ie] );

      bool combineHits = false;

      if ( !isElec ) { // use the photon cut selection

	int origCutSel = allClusters[ie]->getPhotonCutSel();
	int comboCutSel = combo->getPhotonCutSel();

	if ( comboCutSel>=origCutSel ) combineHits = true;

      } else { // use PDF for electron

	combo->setAssociatedTrack( allClusters[ie]->getAssociatedTrack() );

	float newPDF = isElec ? combo->get_total_electronpdf( lowEn ) : combo->get_total_pdf( lowEn ) ;

	float pdf_cut_loose(0);
	float pdf_cut_tight(0);

	if ( isElec ) {
	  pdf_cut_loose = lowEn ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_hiE() ;
	  pdf_cut_tight = lowEn ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_hiE() ;
	} else {
	  pdf_cut_loose = lowEn ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_hiE() ;
	  pdf_cut_tight = lowEn ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_hiE() ;
	}

	combineHits = ( newPDF > oldPDF && newPDF > pdf_cut_loose ) || ( newPDF > pdf_cut_tight );

      }

      if (_verbose) 
	cout << "combining nearby hits to cluster of energy " << allClusters[ie]->getEnergy() << " ? " << combineHits << endl;

      if ( combineHits ) {

	if (_verbose) 
	  cout << "   combined cluster energy = " << combo->getEnergy() << endl;

	allClusters[ie]->addHits( hitsToAdd[ie] );

	// then remove from orig list
	if ( hitsToAdd[ie].size()>0 ) {
	  for (size_t j=0; j<hitsToAdd[ie].size(); j++) {
	    if ( find( unclusteredHits.begin(), unclusteredHits.end(), hitsToAdd[ie][j] ) != unclusteredHits.end() ) {
	      unclusteredHits.erase( find( unclusteredHits.begin(), unclusteredHits.end(), hitsToAdd[ie][j] ) );
	      if ( find( unclusteredHits.begin(), unclusteredHits.end(), hitsToAdd[ie][j] ) != unclusteredHits.end() ) cout << "failed to erase?" << endl;
	    } else {
	      cout << "mergeUnclusteredAndElectrons_byPDF: Strange? shouldn;t end up here..." << endl;
	    }
	  }
	}
      }

      hitsToAdd[ie].clear();
      delete combo;

    }

  } // steps


  return;
}


void ECALGarlicCluster::mergeSatellitesAndElectrons_byPDF(std::map < CalorimeterHit*, ExtendedCluster2* > & clusters ,
							  vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > & allElectrons ) {

  //  cout << " hello from mergeSatellitesAndElectrons_byPDF " << clusters.size() << " " << allElectrons.size() << endl;

  if ( clusters.size()>0 && allElectrons.size()>0 ) {

    std::map < CalorimeterHit*, ExtendedCluster2* >::iterator ittt;

    // first order clusters by energy
    std::map < float, ExtendedCluster2* > enClusters;
    for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
      enClusters[-1*(ittt->second->getEnergy())]=ittt->second; // *-1 to avoid using reverse iterator...
    }
    assert( enClusters.size()==clusters.size() );


    //    cout << "Hello1" << endl;

    std::map < float, ExtendedCluster2* >::iterator itt;

    bool changed=true;
    while (changed) {

      // cout << "Hello2" << endl;

      changed=false;
      for (size_t ie=0; ie<allElectrons.size(); ie++) {

	// cout << "Hello3 " << ie << endl;

	for ( itt=enClusters.begin(); itt!=enClusters.end(); itt++ ) {

	  if ( ! mergeCandidate( allElectrons[ie].second, itt->second ) ) continue;

	  // get orig PDFs
	  bool electron_is_lowE = allElectrons[ie].second->getEnergy()<1.0;
	  float pdfE = allElectrons[ie].second->get_total_electronpdf( electron_is_lowE );

	  //	  bool photon_is_lowE = itt->second->getEnergy()<1.0;
	  //	  float pdfG = itt->second->get_total_pdf( photon_is_lowE );

	  float pdf_cut_loose = electron_is_lowE ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_loosecut_hiE();
	  
	  float pdf_cut_tight = electron_is_lowE ? 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_lowE() : 
	    ECALGarlicAlgorithmParameters::Instance().Get_electronPDF_tightcut_hiE();

	  // try to combine them
	  ExtendedCluster2* combo = new ExtendedCluster2();
	  combo->addHits( *(allElectrons[ie].second->getHits()) );
	  combo->addHits( *(itt->second->getHits()) );

	  combo->setAssociatedTrack( allElectrons[ie].second->getAssociatedTrack() );

	  float newPDF = combo->get_total_electronpdf( electron_is_lowE );

	  // here we will need to decide whether to combine or not
	  bool combine(false);
	  if ( ( newPDF > pdfE && newPDF>pdf_cut_loose ) || // improved pdf
	       ( newPDF > pdf_cut_tight )                   // combined gives a good electron pdf, then combine
	       ) combine = true;

	  if ( combine ) {

	    //cout << "combining electron and photon clusters! photonClusEn=" << itt->second->getEnergy() << " " << itt->second->getHits()->size() << 
	    //  " : " << newPDF << " " << pdfE << " " << pdf_cut_loose << " " << pdf_cut_tight << endl;
	    //
	    //cout << "cluster starts: " << itt->second->getStart() << " " << allElectrons[ie].second->getStart() << endl;


	    // add hits
	    allElectrons[ie].second->addHits( *(itt->second->getHits()) );

	    // cout << "blikka" << endl;

	    for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
	      if ( ittt->second == itt->second ) {
		clusters.erase(ittt);
	      }
	    }
	    enClusters.erase(itt);
	    changed=true;

	    delete combo;
	  
	    // cout << "    merged cluster energy = " << allElectrons[ie].second->getEnergy() << endl;

	    break;
	  } else {

	    delete combo;
	    // cout << "not combining : " << newPDF << " " << pdfE << " " << cogDist << endl;

	  }
	}
	if (changed) break;
      }
    }

    // cout << "Hello6" << endl;

  }



  return;

}



void ECALGarlicCluster::mergeSatellites_byPDF( std::map < CalorimeterHit*, ExtendedCluster2* > & clusters , bool use_PDF_or_cuts ) {

  if ( _verbose ) cout << "hello from mergeSatellites_byPDF " << endl;

  std::map < CalorimeterHit*, ExtendedCluster2* >::iterator ittt;

  const float mass_limit = 0.1; // a little smaller than pi0 mass
  const float max_mass_imbalance_for_pi0=100;


  // first order clusters by energy
  std::map < float, ExtendedCluster2* > enClusters;
  for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
    if ( _verbose ) cout << "orig cluster " << ittt->second->getEnergy() << endl;
    enClusters[-1*(ittt->second->getEnergy())]=ittt->second; // *-1 to avoid using reverse iterator...
  }
  assert( enClusters.size()==clusters.size() );


  std::map < float, ExtendedCluster2* >::iterator itt;
  std::map < float, ExtendedCluster2* >::iterator jtt;

  std::vector < ExtendedCluster2* > deletedClusters;

  bool changed = true;
  while (changed) {
    changed = false;

    for ( itt = enClusters.begin(); itt!=enClusters.end(); itt++) {
      if ( itt==enClusters.begin() ) continue;

      ExtendedCluster2* cl_lo = itt->second;
      if (!cl_lo) continue;

      if ( _verbose ) cout << "trying to combine cluster with energy " << cl_lo->getEnergy() << endl;

      std::map < float, ExtendedCluster2* > allScores;
      std::map < float, ExtendedCluster2* > allDists;

      // look at all higher energy clusters
      for ( jtt = enClusters.begin(); jtt!=enClusters.end(); jtt++) {
	if ( jtt==itt ) break;


	ExtendedCluster2* cl_hi = jtt->second;
	if (!cl_hi) continue;

	bool mergeCand = mergeCandidate( cl_hi, cl_lo );

	if ( _verbose ) cout << "clusters: " << cl_lo->getEnergy() << " " << cl_hi->getEnergy() << " merge candidate? " << mergeCand << endl;

	if ( ! mergeCand ) 
	  continue;

	std::pair < float, float > comboMassErr = cl_hi->getCombinedInvariantMass(cl_lo);

	if ( comboMassErr.first+2*comboMassErr.second>mass_limit &&
	     cl_hi->getEnergy()/cl_lo->getEnergy() < max_mass_imbalance_for_pi0 ) {

	  if ( _verbose ) {
	    cout << "refusing to combine clusters whose mass is close to pizero mass! " << 
	      comboMassErr.first << " " << comboMassErr.second << " " <<
	      comboMassErr.first+2*comboMassErr.second << " " << mass_limit  << endl;
	    cout << " energy imbalance = " << cl_hi->getEnergy()/cl_lo->getEnergy() << " limit " << max_mass_imbalance_for_pi0 << endl;
	  }

	  continue;

	}

	float pdf_hi = cl_hi->getEnergy() > 1 ? cl_hi->get_total_pdf(false) : cl_hi->get_total_pdf(true);

	ExtendedCluster2* combo = new ExtendedCluster2();
	combo->addHits( *(cl_hi->getHits()) );
	combo->addHits( *(cl_lo->getHits()) );

	float pdf_combo = cl_hi->getEnergy() > 1 ? combo->get_total_pdf(false) : combo->get_total_pdf(true); // decide on cl_hi energy, to avoid flipping between hi and lo in comparison

	delete combo; combo=NULL;

	bool comboPassesLoose, comboPassesTight;

	float pdf_cut_tight = cl_hi->getEnergy() > 1 ? 
	  ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_hiE() : 
	  ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_lowE() ;

	float pdf_cut_loose = cl_hi->getEnergy() > 1 ? 
	  ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_hiE() : 
	  ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_lowE() ;

	float comboLikeRatio = pdf_combo/pdf_hi;

	if ( use_PDF_or_cuts ) {
	  comboPassesLoose = pdf_combo>pdf_cut_loose;
	  comboPassesTight = pdf_combo>pdf_cut_tight;
	} else {
	  int cutsel = cl_hi->getPhotonCutSel();
	  comboPassesTight = cutsel>=2;
	  comboPassesLoose = cutsel>=1;
	}

	float score = 100*comboPassesTight + 10*comboPassesLoose + 1./comboLikeRatio; // the higher the score, the nicer the combination

	if ( _verbose ) cout << "score: " << cl_hi->getEnergy() << " " << comboPassesTight << " " << comboPassesLoose << " " << comboLikeRatio << " : " << score << endl;

	float dist = cl_hi->getDistToClusterAxis( cl_lo->getCentreOfGravity(), 0 );
	if ( _verbose ) cout << "    distance: " << dist << endl;

	if ( comboLikeRatio<2 ) { // don't make the likelihood significantly worse...
	  allDists[dist] = cl_hi;
	  allScores[score] = cl_hi;
	} else {
	  if ( _verbose ) cout << "rejecting combo which significantly degrades likelihood..." << comboLikeRatio << " " << pdf_combo << " " << pdf_hi << endl;
	}

      }

      if ( allScores.size()>0 ) {
	std::map < float, ExtendedCluster2* >::reverse_iterator kj;
	std::map < float, ExtendedCluster2* >::iterator kjj;

	if ( _verbose ) {
	  cout << "all scores for this cluster: parent energy " << cl_lo->getEnergy() << endl;
	  for ( kj=allScores.rbegin(); kj!=allScores.rend(); kj++) {
	    cout << "   child energy, score     " << kj->second->getEnergy() << " " << kj->first << endl;
	  }
	  
	  cout << "   and the distances..." << endl;
	  for ( kjj=allDists.begin(); kjj!=allDists.end(); kjj++) {
	    cout << "   child energy, distance  " << kjj->second->getEnergy() << " " << kjj->first << endl;
	  }
	}

	ExtendedCluster2* bestCluster(0);
	float bestScore(0);
	float bestDist(0);

	if ( allScores.size()==1 ) {
	  bestCluster = allScores.rbegin()->second;
	  bestScore = allScores.rbegin()->first;
	  bestDist = allDists.begin()->first;
	} else if ( allScores.rbegin()->second == allDists.begin()->second ) { // same chosen by both measures
	  bestCluster = allScores.rbegin()->second;	  
	  bestScore = allScores.rbegin()->first;
	  bestDist = allDists.begin()->first;
	} else {  // some controversy! closest cluster doesn;t give best score!

	  if ( _verbose ) cout << "controversy!" << endl;

	  // look at 2 closest clusters
	  kjj=allDists.begin();
	  float dist1 = kjj->first;
	  kjj++;
	  float dist2 = kjj->first;
	  //	  cout << dist1 << " " << dist2 << endl;
	  // if ( dist2 - dist1 > ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius() ) { // difference in distance is large, choose closest one
	  if ( 1 ) { 
	    // choose by distance
	    bestCluster = allDists.begin()->second;
	    bestDist = allDists.begin()->first;
	    for ( kj=allScores.rbegin(); kj!=allScores.rend(); kj++) {
	      if ( kj->second == allDists.begin()->second ) {
		bestScore=kj->first;
		break;
	      }
	    }
	    //	    cout << "choosing by distance, score = " << bestScore << endl;
	  } else { // not so much difference in distance, choose best score
	    // choose by score
	    bestCluster = allScores.rbegin()->second;
	    bestScore = allScores.rbegin()->first;
	    bestDist = kjj->first;

	    //	    cout << "choosing by score..." << endl;
	  }

	}


	if ( bestCluster ) { // then combine
	
	  if ( _verbose ) cout << "combining! best score = " << bestScore << endl;
	
	  changed=true;
	  bestCluster->addHits( *(cl_lo->getHits()) );
	  for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {    
	    if ( ittt->second == cl_lo ) {
	      clusters.erase(ittt);// delete entry from orig map
	    }
	  }
	  enClusters.erase(itt);  // delete entry from the energy-ordered map
	  delete cl_lo;  // delete the cluster object itself
	}
      }

      if (changed) break;

    } // itt loop


    // for ( itt = enClusters.begin(); itt!=enClusters.end(); itt++) {
    // 
    //   // also don't want penultimate one
    //   jtt=itt;
    //   jtt++;
    //   if (jtt==enClusters.end()) continue;
    // 
    //   ExtendedCluster2* cl_i = itt->second;
    //   if (!cl_i) continue;
    // 
    //   for ( jtt = itt; jtt!=enClusters.end(); jtt++) {
    // 
    // 	if (jtt==itt) continue;
    // 	ExtendedCluster2* cl_j = jtt->second;
    // 	if (!cl_j) continue;
    // 
    // 	bool mergeCand = mergeCandidate( cl_i, cl_j );
    // 	if ( ! mergeCand ) continue;
    // 
    // 	// get their combined invariant mass
    // 	// if it's close to or above the pi0, refuse to combine! (at least if energy imbalance is not enormous)
    // 	std::pair < float, float > comboMassErr = cl_i->getCombinedInvariantMass(cl_j);
    // 	cout << "considering combining 2 nearby photon clusters, energies: " << cl_i->getEnergy() << " " << cl_j->getEnergy() << 
    // 	  " combined mass " << comboMassErr.first << "+/-" << comboMassErr.second << endl;
    // 
    // 	if ( comboMassErr.first+2*comboMassErr.second>mass_limit &&            // not too inconsistent with pi0 mass
    // 	     cl_i->getEnergy()/cl_j->getEnergy() < max_mass_imbalance_for_pi0 ) {  // not completely unbalanced in energy
    // 	  cout << "high combo mass, refusing to combine! "
    // 	       << comboMassErr.first << " +/- " << comboMassErr.second << " , "
    // 	       <<  cl_i->getEnergy()/ cl_j->getEnergy() << endl;
    // 	  continue;
    // 	}
    // 
    // 	cout << "    start1, start2 = " << cl_i->getStart() << " " << cl_j->getStart() << endl;
    // 
    // 	float pdfLo_1 = cl_i->get_total_pdf(true);
    // 	float pdfHi_1 = cl_i->get_total_pdf(false);
    // 
    // 	float pdfLo_2 = cl_j->get_total_pdf(true);
    // 	float pdfHi_2 = cl_j->get_total_pdf(false);
    // 	
    // 	cout << "comparing 2 clusters: " << endl;
    // 	cout << "first PDFs: " <<  cl_i->getEnergy() << " " << cl_i->getHits()->size() << " " << pdfLo_1 << " " << pdfHi_1 << endl;
    // 	cout << "second PDFs: " << cl_j->getEnergy() << " " << cl_j->getHits()->size() << " " << pdfLo_2 << " " << pdfHi_2 << endl;
    // 
    // 	// try to combine them
    // 	ExtendedCluster2* combo = new ExtendedCluster2();
    // 	combo->addHits( *(cl_i->getHits()) );
    // 	combo->addHits( *(cl_j->getHits()) );
    // 
    // 	float pdfLo_C = combo->get_total_pdf(true);
    // 	float pdfHi_C = combo->get_total_pdf(false);
    // 
    // 	cout << "combo PDFs: " << combo->getEnergy() << " : " << pdfLo_C << " " << pdfHi_C << endl;
    // 
    // 	bool combine=false;
    // 
    // 
    // 	float likeRatio = cl_i->getEnergy()>=1 ? pdfHi_C/pdfHi_1 : pdfLo_C/pdfLo_1;
    // 	bool comboPassesLoose = cl_i->getEnergy()>=1 ? 
    // 	  pdfHi_C>ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_hiE() : 
    // 	  pdfLo_C>ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_lowE() ;
    // 	bool comboPassesTight = cl_i->getEnergy()>=1 ? 
    // 	  pdfHi_C>ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_hiE() : 
    // 	  pdfLo_C>ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_lowE() ;
    // 
    // 	cout << "like ratio = " << likeRatio << " " << comboPassesLoose << " " << comboPassesTight << endl;
    // 	
    // 	combine = 
    // 	  ( likeRatio<1.0 ) ||
    // 	  ( likeRatio<1.2 && comboPassesLoose ) || 
    // 	  comboPassesTight ;
    // 
    // 	// // does adding the cluster improve the likelihood of the larger cluster? does the combined cluster pass a loose cut?
    // 	// if ( cl_i->getEnergy()>=1 ) {
    // 	//   if ( 
    // 	//       // pdfHi_C > ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_hiE() && // with small seeds, sometimes showers so narrow they fais cut
    // 	//       pdfHi_C > pdfHi_1 ) combine=true;
    // 	// } else {
    // 	//   if ( 
    // 	//       //pdfLo_C > ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_loosecut_lowE() && 
    // 	//       pdfLo_C > pdfLo_1 ) combine=true;
    // 	// }
    // 	// 
    // 	// // if the combined cluster passes a tighter cut, we choose to combine, even if the likelihood is less than the original
    // 	// if ( cl_i->getEnergy()>=1 ) {
    //     //   if ( pdfHi_C > ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_hiE() ) combine=true;
    // 	// } else {
    // 	//   if ( pdfLo_C > ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_lowE() ) combine=true;
    // 	// }
    // 
    // 	if (combine) {
    // 
    // 	  cout << "Combining!" << endl;
    // 	  // cout << "hello1, I wanna combine 2 clusters! cogdist" << endl;
    // 	  // cout << "nhits before " << cl_i->getHits()->size() << " " << cl_j->getHits()->size() << endl;
    // 	  // cout << pdfHi_C << " " << pdfHi_1 << " " << ECALGarlicAlgorithmParameters::Instance().Get_photonPDF_tightcut_hiE() << endl;
    // 
    // 	  cl_i->addHits( *(cl_j->getHits()) );
    // 
    // 	  // cout << "and after: " << cl_i->getHits()->size() << endl;
    // 
    // 	  // delete entry from orig map
    // 	  for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {    
    // 	    if ( ittt->second == cl_j ) {
    // 	      clusters.erase(ittt);
    // 	    }
    // 	  }
    // 
    // 	  // delete entry from the energy-ordered map
    // 	  enClusters.erase(jtt);
    // 
    // 	  // delete the cluster object itself
    // 	  delete cl_j;
    // 
    // 	  // restart the second loop to compare with combined cluster
    // 	  jtt=itt;
    // 	  changed = true;
    // 	  delete combo;
    // 	  break;
    // 	} else {
    // 	  delete combo;
    // 	}
    // 
    //   } // second cluster loop
    // 
    //   if (changed) break;
    // 
    // } // first cluster loop

    //    cout << "changed? " << changed << endl;

  } // while changed

  // remove the dead clusters from the original list

  // cout << "final n clusters after PDF merge..." << clusters.size() << endl;
  // for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
  //   cout << ittt->second << " " << ittt->second->getEnergy() << " " << ittt->second->getHits()->size() << endl;
  // }
  // cout << "done PDF merge!" << endl;

  return;
}


bool ECALGarlicCluster::mergeCandidate( ExtendedCluster2* primary, ExtendedCluster2* secondary ) {
  // chg+neu: primary should be charged with neutral secondary partner
  // neu+neu: primary should have higher energy

  if ( secondary->getAssociatedTrack() ) {
    if (  primary->getAssociatedTrack() ) {
      cout << "asking ECALGarlicCluster::mergeCandidate to make judgement about two charged clusters...it's not designed for this!" << endl;
      return false;
    } else { // switch to put charged in primary position
      ExtendedCluster2* temp = primary;
      primary=secondary;
      secondary=temp;
    }
  } else if ( !primary->getAssociatedTrack() ) { // 2 neutral clusters
    if ( primary->getEnergy() < secondary->getEnergy() ) { // make sure primary has higher energy
      ExtendedCluster2* temp = primary;
      primary=secondary;
      secondary=temp;
    }
  }

  // should make these parameters!
  const float ratioCut = 0.2;
  const float energy_dist_factor = 120;
  const float distanceMultiplier = 1.;

  float cogDist =  primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 0 ); // ip pointing
  if ( primary->getAssociatedTrack() ) cogDist =  std::min( cogDist, primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 1 )); // track direction
  cogDist =  std::min( cogDist, primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 2 )); // cluster axis

  float energyRatio = secondary->getEnergy()/primary->getEnergy();

  //  cout << "deciding if I should merge: dist: " << cogDist << " enRatio: " << energyRatio << " prinamy energy: " << primary->getEnergy() << endl;
  //  cout << "   " << energyRatio << " " << ratioCut << " " << bool(energyRatio>ratioCut) << endl;
  //  cout << "   " << cogDist << " " << distanceMultiplier*ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist() << " " << 
  // bool(cogDist<distanceMultiplier*ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist()) << endl;
  //  cout << "   " << energyRatio << " " << energy_dist_factor << " " << energy_dist_factor/pow(cogDist,2) << " " << bool(energyRatio<energy_dist_factor/pow(cogDist,2)) << endl;

  return energyRatio<ratioCut && // don't merge large energy clusters
    cogDist<distanceMultiplier*ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist()*ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius() && // nor very far apart ones
    energyRatio<energy_dist_factor/pow(cogDist,2);  // energy fraction-dependent distance cut
}



void ECALGarlicCluster::mergeSatellites( std::map < CalorimeterHit*, ExtendedCluster2* > & clusters ) {

  // cout << "hello from mergeSatellites" << endl;

  // the parameters
  float touchingLayerFraction = ECALGarlicAlgorithmParameters::Instance().GetMergeTouchFraction();
  int   initialSeparationCut = ECALGarlicAlgorithmParameters::Instance().GetMergeInitalLayerSeparation();

  std::map < CalorimeterHit*, ExtendedCluster2* >::iterator itt;
  std::map < CalorimeterHit*, ExtendedCluster2* >::iterator jtt;

  const float touchDist = 1.5*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

  // compare all cluster pairs

  bool changed = true;
  while (changed) {
    changed = false;

    for ( itt = clusters.begin(); itt!=clusters.end(); itt++) {

      // also don't want penultimate one
      jtt=itt;
      jtt++;
      if (jtt==clusters.end()) continue;

      ExtendedCluster2* cl_i = itt->second;
      if (!cl_i) continue;

      std::vector<ExtendedHit2*>* hits_i = cl_i->getHits();

      for ( jtt = clusters.begin(); jtt!=clusters.end(); jtt++) {

	if (jtt==itt) continue;

	ExtendedCluster2* cl_j = jtt->second;
	if (!cl_j) continue;


	float cogDist(0);
	for (int ik=0; ik<3; ik++) 
	  cogDist+=pow( cl_i->getCentreOfGravity()[ik] - cl_j->getCentreOfGravity()[ik], 2);
	cogDist=sqrt(cogDist);

	if (cogDist>
	    ECALGarlicAlgorithmParameters::Instance().Get_MaxMergeDist()*ECALGarlicAlgorithmParameters::Instance().GetMoliereRadius()) continue; // mm

	//	cout << "trying to merge clusters with energies " << cl_i->getEnergy() << " " << cl_j->getEnergy() << endl;

	// find min hit-hit distance in each p-layer
	std::map < int, float > dist_player;

	std::vector<ExtendedHit2*>* hits_j = cl_j->getHits();

	for ( size_t i=0; i<hits_i->size(); i++) {
	  ExtendedHit2* h_i = (*hits_i)[i];
	  int pl_i = h_i->getPseudoLayer();
	  for ( size_t j=0; j<hits_j->size(); j++) {
	    ExtendedHit2* h_j = (*hits_j)[j];
	    int pl_j = h_j->getPseudoLayer();
	    if (pl_i!=pl_j) continue;
	    
	    float dist = 0;
	    for (int kk=0; kk<3; kk++) 
	      dist+=pow( h_i->getCaloHit()->getPosition()[kk] - h_j->getCaloHit()->getPosition()[kk], 2);
	    dist = sqrt( dist );

	    if (dist_player.find(pl_i)==dist_player.end()) dist_player[pl_i] = dist;
	    else                                           dist_player[pl_i] = std::min ( dist_player[pl_i], dist);
	  }
	}

	// count no of close-by-layers
	int ntouch(0);
	int initialSeparation(0);
	bool stillSeparated(true);

	for ( std::map < int, float > ::iterator ipp=dist_player.begin(); ipp!=dist_player.end(); ipp++ ) {
	  if ( ipp->second < touchDist ) {
	    stillSeparated=false;
	    ntouch++;
	  } else {
	    if (stillSeparated) initialSeparation++;
	  }
	}

	//	cout << "mergeSatellites:: pair of clusters: min distances " << endl;
	//	for ( std::map < int, float > ::iterator ipp=dist_player.begin(); ipp!=dist_player.end(); ipp++ ) {
	//	  cout << ipp->first << " " << ipp->second << endl;
	//	}
	//	cout << "ntouching layers = " << ntouch << " / " << dist_player.size() << endl;
	//	cout << "separated initial layers = " << initialSeparation << endl;
	
	if ( float(ntouch)/dist_player.size() > touchingLayerFraction && initialSeparation < initialSeparationCut ) { // merge the clusters

	  // cout << "MERGING clusters! " << endl;

	  cl_i->addHits(*hits_j);
	  delete cl_j;
	  jtt->second = NULL;
	  
	  // restart the second loop to compare with combined cluster
	  jtt=itt;

	  changed = true;
	  break;
	}

      } // second cluster loop

      if (changed) break;

    } // first cluster loop
  }

  //  clusters.erase(clusters.begin(), clusters.end());

  for (itt=clusters.end(); itt!=clusters.begin(); itt--) {
    if (itt==clusters.end()) continue;
    if (!itt->second) clusters.erase(itt);
  }

  // just double check...
  for ( itt = clusters.begin(); itt!=clusters.end(); itt++) {
    if (!itt->second) cout << "ERROR still a deleted cluster in the list!" << endl;
  }

  return;
}



