#include "ECALGarlicCluster.hh"

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>

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
  float cutoffdist = ECALGarlicAlgorithmParameters::Instance().GetSeedDistanceCut();
  int maxpseulayer = ECALGarlicAlgorithmParameters::Instance().GetSeedNLayers();

  const int maxseed=300;

  //  std::vector <CalorimeterHit*> seeds;
  std::map <CalorimeterHit*, bool> seeds;

  std::pair < TH2F*, TH2I* > histos = preClus->GetProjectionHistos( energyCutMip, maxpseulayer );
  TH2F* henergy = histos.first;
  TH2I* hhits = histos.second;

  if (!henergy) {
    return seeds;
  }

  TString blah;

  //  if (henergy && hhits && _fhistos) {

  if (_fhistos) _fhistos->cd();

  if (_hnn!="blah") blah="_"+_hnn+"_";
  else blah="_n_";
  blah+="nhits";
  blah+=preClus->getHits()->size();
  blah+="_";
  blah+=_nSaveHist+1;
    
  henergy = (TH2F*) henergy->Clone("hen"+blah);
  henergy->SetNameTitle("hen"+blah, "hen"+blah);

  hhits = (TH2I*) hhits->Clone("hhit"+blah);
  hhits->SetNameTitle("hhit"+blah, "hhit"+blah);

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

  int iseed=0;
  float seedXvec[maxseed]={0};
  float seedYvec[maxseed]={0};

  int ibadseed=0;
  float badseedXvec[maxseed]={0};
  float badseedYvec[maxseed]={0};

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
      if (iseed<maxseed) {
	seedXvec[iseed]=seed_x;
	seedYvec[iseed]=seed_y;
	iseed++;
      } else {
	cout << "WARNING, too many seeds, only keeping first " << maxseed << endl;
	cout << " increase maxseed in ECALGarlicCluster::getSeeds" << endl;
      }
    } else {
      if (ibadseed<maxseed) {
	badseedXvec[ibadseed]=seed_x;
	badseedYvec[ibadseed]=seed_y;
	ibadseed++;
      }
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
    TGraph* gr = new TGraph(iseed, seedXvec, seedYvec);
    gr->SetNameTitle("seeds"+blah,"seeds"+blah);
    gr->Write();

    TGraph* gr2 = new TGraph(ibadseed, badseedXvec, badseedYvec);
    gr2->SetNameTitle("badseeds"+blah,"badseeds"+blah);
    gr2->Write();
    _nSaveHist++;
  }

  if (henergy) delete henergy;
  if (hhits) delete hhits;

  return seeds;
}

vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > ECALGarlicCluster::getElectrons( ExtendedCluster2* preCluster, vector <ExtendedTrack* > trks ) {

  //  cout << "hello from ECALGarlicCluster::getElectrons "<< trks.size() << " " << preCluster->getNhits() << endl;
  
  vector < std::pair < ExtendedTrack*, ExtendedCluster2* > > electrons;

  float addCut = 1.9*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  int nIterations = ECALGarlicAlgorithmParameters::Instance().GetClusterNIterations();
  const float origin[3]={0,0,0};
  
  for ( size_t it=0; it<trks.size(); it++ ) {

    // position and direction of track at ECAL front face
    const float* position  = trks[it]->getEcalEntryPos();
    const float* direction = trks[it]->getEcalEntryDir();

    // build a "seed"
    CalorimeterHitImpl seed;
    seed.setEnergy(1);
    seed.setPosition(position);

    // make the core
    ExtendedCluster2* electronCore = BuildCore( &seed, preCluster, direction);

    // find hits not too far from core axis
    std::vector < ExtendedHit2* > possibleHits;
    for ( size_t ij=0; ij<preCluster->getHits()->size(); ij++) {
      ExtendedHit2* hh = preCluster->getHits()->at(ij);
      if ( find(possibleHits.begin(), possibleHits.end(), hh)!=possibleHits.end() ) {
	float dist = ECALGarlicGeometryHelpers::GetDistToLine(hh->getCaloHit()->getPosition(), origin, electronCore->getCentreOfGravity());
	if (dist<40.) possibleHits.push_back(hh);
      }
    }

    for (int iiter=0; iiter<nIterations; iiter++) {
      std::vector < ExtendedHit2* > hitsToAdd;
      for ( int ij=possibleHits.size()-1; ij>=0; ij--) {
	
	std::vector < ExtendedHit2* > * corehits = electronCore->getHits();

	for (size_t ik=0; ik<corehits->size(); ik++) {
	  ExtendedHit2* chit = corehits->at(ik);
	  float dist(0);
	  for (int i=0; i<3; i++) 
	    dist += pow ( possibleHits[ij]->getPosition()[i] - chit->getCaloHit()->getPosition()[i], 2);
	  dist = sqrt(dist);
	  if (dist<addCut) {
	    hitsToAdd.push_back(chit);
	    break;
	  }
	}

      }

      if (hitsToAdd.size()>0) {
	electronCore->addHits(hitsToAdd);
      }

    }

    float EonP = electronCore->getEnergy()/trks[it]->getTotalMomentum();

    // make cuts to find electron clusters
    if (electronCore->getNhits()>0 && 
	EonP>0.5 && EonP<1.5 &&
	electronCore->getLongEn()[0]>0.1 && electronCore->getLongEn()[0]+electronCore->getLongEn()[1]>0.5
	) {

      // cout << "found a cluster seeded by track" << endl;
      // cout << " track momentum = " << trks[it]->getTotalMomentum() << " ";
      // cout << " nhits = " << electronCore->getNhits() << " energy = " << electronCore->getEnergy() << " ";
      // cout << " E/p = " << EonP << endl;
      // cout << " long fracs = " << electronCore->getLongEn()[0] << " " << electronCore->getLongEn()[1] ;
      // cout << " " << electronCore->getLongEn()[2] << " " << electronCore->getLongEn()[3] << endl;

      electrons.push_back( std::pair < ExtendedTrack*, ExtendedCluster2* > ( trks[it], electronCore ) );
    }

  } // loop over tracks


  return electrons;
}




std::map < CalorimeterHit*, ExtendedCluster2* > ECALGarlicCluster::getClusters( ExtendedCluster2* preCluster, std::map < CalorimeterHit*, ExtendedCluster2* > cores ) {
  // assign unassigned hits to "best" core

  float addCut = 1.9*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  int nIterations = ECALGarlicAlgorithmParameters::Instance().GetClusterNIterations();

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

  for (int iiter=0; iiter<nIterations; iiter++) {
    std::map < CalorimeterHit*, std::vector < ExtendedHit2* > > preassignedHits;
    // find distance to closest hit in core
    for ( int ihh=unassignedHits.size()-1; ihh>=0; ihh-- ) {
      float mindist=99999;
      CalorimeterHit* closestCore=0;
      for ( std::map < CalorimeterHit*, ExtendedCluster2* >::iterator clit = clusters.begin(); clit!=clusters.end(); clit++ ) {
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
  float distCut = 1.5*ECALGarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

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

  return coreClus;  
}


void ECALGarlicCluster::mergeSatellites( std::map < CalorimeterHit*, ExtendedCluster2* > & clusters ) {

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

	if (cogDist>100.0) continue; // mm

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



