#include "GarlicClusterAlgos.hh"

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

#include "GarlicAlgorithmParameters.hh"
#include "GarlicGeometryParameters.hh"
#include "GarlicGeometryHelpers.hh"

#include "GarlicExtendedCluster.hh"
#include "GarlicExtendedHit.hh"

std::map <CalorimeterHit*, bool> GarlicClusterAlgos::getSeeds(GarlicExtendedCluster* preClus) {

  float energyCutMip = GarlicAlgorithmParameters::Instance().GetSeedHitEnergyCut();
  float seedThresholdMIP = GarlicAlgorithmParameters::Instance().GetSeedEnergyCut();
  float cutoffdist =
    GarlicAlgorithmParameters::Instance().GetSeedDistanceCut()*
    GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  if ( cutoffdist>GarlicAlgorithmParameters::Instance().GetMoliereRadius() ) {
    cout << "WARNING, requested scale for seeding is larger than Moliere radius...setting to Moliere radius: " <<
      " distance (cells)= " << GarlicAlgorithmParameters::Instance().GetSeedDistanceCut() <<
      " cell size = " << GarlicGeometryParameters::Instance().Get_padSizeEcal()[0] <<
      " distance (mm) = " << GarlicAlgorithmParameters::Instance().GetSeedDistanceCut()*
      GarlicGeometryParameters::Instance().Get_padSizeEcal()[0] <<
      " Moliere: " << GarlicAlgorithmParameters::Instance().GetMoliereRadius() << endl;
    cutoffdist=GarlicAlgorithmParameters::Instance().GetMoliereRadius();
  }

  int maxpseulayer = GarlicAlgorithmParameters::Instance().GetSeedNLayers();

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

  // need to clone the histos and rebin them to approx cell size
  float cellsize = GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  int ics = int(cellsize / henergy->GetXaxis()->GetBinWidth(1) );

  henergy->Rebin2D(ics, ics);
  hhits->Rebin2D(ics, ics);

  std::vector < std::pair < int, int > > seedCells;

  seedCells.clear();

  // find highest energy bin
  int maxbin = henergy->GetMaximumBin();
  float maxen = henergy->GetBinContent( maxbin );

  std::vector <float> seedXvec;
  std::vector <float> seedYvec;

  std::vector <float> badseedXvec;
  std::vector <float> badseedYvec;

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

    bool goodSeed = seed_nhits>=GarlicAlgorithmParameters::Instance().GetSeedMinHits() && toten>=seedThresholdMIP ;

    float pos[3];
    preClus->GetGlobalPositionFromLocal(seed_x, seed_y, pos);

    CalorimeterHitImpl* a_seed=new CalorimeterHitImpl();
    a_seed->setEnergy(seed_en);
    a_seed->setPosition(pos);
    seeds[a_seed] = goodSeed;

    if (goodSeed) {
      seedXvec.push_back(seed_x);
      seedYvec.push_back(seed_y);
    } else {
      badseedXvec.push_back(seed_x);
      badseedYvec.push_back(seed_y);
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

vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > GarlicClusterAlgos::getElectrons( GarlicExtendedCluster* preCluster, vector <GarlicExtendedTrack* > trks, bool force_a_cluster ) {

  // try to make a narrow cluster around the track, and check if it looks like an electron
  // force_a_cluster overrides the electron selection (we will always make a cluster attached to each track)
  //    (this feature intended to be used on conversion tracks)

  if (_verbose)
    cout << "hello from GarlicClusterAlgos::getElectrons; ntracks="<< trks.size() << " ; nECAL hits " << preCluster->getNhits() << endl;

  vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > electrons;

  // try a narrow tube
  //  float maxAddDist = 0.75*GarlicAlgorithmParameters::Instance().GetMoliereRadius();
  float maxAddDist = GarlicAlgorithmParameters::Instance().GetElectronTransTubeStepSize()*GarlicAlgorithmParameters::Instance().GetMoliereRadius();

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
    GarlicExtendedCluster* electronCore = BuildCore( &seed, preCluster, direction);
    electronCore->setAssociatedTrack(trks[it]); // ref direction of cluster is track direction at eacl entry

    if (_verbose)
      cout << "core hits: " << electronCore->getHits()->size() << endl;

    if ( electronCore->getHits()->size() ==0 ) {
      delete electronCore;
      electronCore=NULL;
      continue;
    }

    if (_verbose) {
      cout << "electron core selected ? " << electronCore->getElectronCutSel() << " : ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::TRANS ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::LONGT  ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::MISCL ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::POINT ) << " ";
      cout << endl;
    }

    // does the core look like an electron in long. direction? if not, reject
    if ( electronCore->getElectronCutScore(GarlicExtendedCluster::LONGT ) == 0 ) {
      if (_verbose) {
        cout << "rejecting electron core which doesn't look like an electron..." << endl;
      }
      delete electronCore;
      electronCore=NULL;
      continue;
    }


    if (_verbose) {
      cout << "selected electron core..." << endl;
    }

    // find hits not too far from core axis
    float point2[3]; // second point along shower axis
    for (int i=0; i<3; i++)
      point2[i]=electronCore->getCentreOfGravity()[i]-direction[i];

    int nsteps=GarlicAlgorithmParameters::Instance().GetElectronTransNSteps();
    //    const int nsteps=3;

    std::vector < std::vector < GarlicExtendedHit* >  > possibleHits;
    for (int i=0; i<nsteps; i++) {
      std::vector < GarlicExtendedHit* > aa;
      possibleHits.push_back(aa);
    }

    // std::vector < GarlicExtendedHit* > possibleHits[nsteps];
    for ( size_t ij=0; ij<preCluster->getHits()->size(); ij++) {
      GarlicExtendedHit* hh = preCluster->getHits()->at(ij);
      if ( find(electronCore->getHits()->begin(), electronCore->getHits()->end(), hh)==electronCore->getHits()->end() ) {
        // this assumed projective to IP...
        //      float dist = GarlicGeometryHelpers::GetDistToLine(hh->getCaloHit()->getPosition(), origin, electronCore->getCentreOfGravity());
        float dist = GarlicGeometryHelpers::GetDistToLine(hh->getCaloHit()->getPosition(), point2, electronCore->getCentreOfGravity());

        for (int kk=0; kk<nsteps; kk++) {
          if ( dist < (kk+2.)*maxAddDist/nsteps ) {
            possibleHits[kk].push_back(hh);
            break;
          }
        }

      }
    }

    if (_verbose) {
      for (int kk=0; kk<nsteps; kk++) {
        cout << " got " << possibleHits[kk].size() << " possible hits within a distance " << (kk+2)*maxAddDist/nsteps << " from axis " << endl;
      }
    }

    // at each iteration, we should check that we don't screw up electron properties.
    for (int istep=0; istep<nsteps; istep++) {

      GarlicExtendedCluster* newcluster = new GarlicExtendedCluster();
      newcluster->addHits(*(electronCore->getHits()));
      newcluster->addHits(possibleHits[istep]);
      newcluster->setAssociatedTrack(electronCore->getAssociatedTrack());
      newcluster->setReferenceDirection(electronCore->getReferenceDirection());

      float prevScore = electronCore->getElectronCutSel();
      // electronCore->getElectronCutScore(GarlicExtendedCluster::TRANS ) +
      //   electronCore->getElectronCutScore(GarlicExtendedCluster::LONGT ) +
      //   electronCore->getElectronCutScore(GarlicExtendedCluster::MISCL ) +
      //   electronCore->getElectronCutScore(GarlicExtendedCluster::POINT ) ;

      float updatedScore = newcluster->getElectronCutSel();
      //  newcluster->getElectronCutScore(GarlicExtendedCluster::TRANS ) +
      //  newcluster->getElectronCutScore(GarlicExtendedCluster::LONGT ) +
      //  newcluster->getElectronCutScore(GarlicExtendedCluster::MISCL ) +
      //  newcluster->getElectronCutScore(GarlicExtendedCluster::POINT ) ;

      if (_verbose) {
        cout << "previous core selected ? " ;
        cout << electronCore->getElectronCutSel() << " : ";
        cout << electronCore->getElectronCutScore(GarlicExtendedCluster::TRANS ) << " ";
        cout << electronCore->getElectronCutScore(GarlicExtendedCluster::LONGT  ) << " ";
        cout << electronCore->getElectronCutScore(GarlicExtendedCluster::MISCL ) << " ";
        cout << electronCore->getElectronCutScore(GarlicExtendedCluster::POINT ) << " : ";
        cout << prevScore;
        cout << endl;

        cout << "updated core selected ? " ;
        cout << newcluster->getElectronCutSel() << " : ";
        cout << newcluster->getElectronCutScore(GarlicExtendedCluster::TRANS ) << " ";
        cout << newcluster->getElectronCutScore(GarlicExtendedCluster::LONGT ) << " ";
        cout << newcluster->getElectronCutScore(GarlicExtendedCluster::MISCL ) << " ";
        cout << newcluster->getElectronCutScore(GarlicExtendedCluster::POINT ) << " : ";
        cout << updatedScore;
        cout << endl;
      }

      if ( updatedScore >= prevScore ) {
        if (_verbose)
          cout << "UPDATING ELECTRON CORE!" << endl;
        electronCore->addHits(possibleHits[istep]);
        delete newcluster;
        newcluster=NULL;
      } else {
        if (_verbose)
          cout << "NOT UPDATING ELECTRON CORE! stop adding" << endl;
        delete newcluster;
        newcluster=NULL;
        break;
      }

    }


    float EonP = electronCore->getEnergy()/electronCore->getAssociatedTrack()->getTotalMomentum();

    int cutSel = int(electronCore->getElectronCutSel());

    if (_verbose) {
      cout << "final electron" << endl;
      cout << " core hits " << electronCore->getNhits() << " start " << electronCore->getStart() << " cut sel? " << cutSel << " e/p " << EonP << endl;
      cout << " selected ? " << endl;
      cout << electronCore->getElectronCutSel() << " : ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::TRANS ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::LONGT ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::MISCL ) << " ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::POINT ) << " : ";
      cout << electronCore->getElectronCutScore(GarlicExtendedCluster::E_ONP ) << " : ";
      cout << endl;
    }

    // make cuts to find electron clusters
    if ( electronCore->getNhits()>0 && ! electronCore->getIsMipLike() && ( cutSel > 0 || force_a_cluster ) ) {

      if (_verbose) {
        cout << "accepting electron with e/p = " << EonP << " , mip-like ? " << electronCore->getIsMipLike() << endl;
      }

      electrons.push_back( std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > ( trks[it], electronCore ) );
      trks[it]->setElectronSel( 1 );

    } else {

      if (_verbose) 
	cout << "rejecting cluster with energy " << electronCore->getEnergy() << endl;

      // delete the new cluster
      delete electronCore; electronCore=NULL;
    }
  } // loop over tracks

  return electrons;
}




std::map < CalorimeterHit*, GarlicExtendedCluster* > GarlicClusterAlgos::getClusters( GarlicExtendedCluster* preCluster, std::map < CalorimeterHit*, GarlicExtendedCluster* > cores ) {
  // assign unassigned hits to "best" core

  //  float addCut = 1.9*GarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  //  float addCut = 2.4*GarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  float addCut = GarlicAlgorithmParameters::Instance().GetTouchingCellDistance()*GarlicGeometryParameters::Instance().Get_padSizeEcal()[1]; // include diagonal but not 2 cells away
  //  int nIterations = GarlicAlgorithmParameters::Instance().GetClusterNIterations();

  float maxDistCut = GarlicAlgorithmParameters::Instance().GetClusterMaxDist()*GarlicAlgorithmParameters::Instance().GetMoliereRadius();

  // duplicate the cores into the clusters
  std::map < CalorimeterHit*, GarlicExtendedCluster* > clusters;
  for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator clit = cores.begin(); clit!=cores.end(); clit++ ) {
    GarlicExtendedCluster* cluster = new GarlicExtendedCluster();
    (*cluster) = *(clit->second);
    clusters[clit->first] = cluster;
  }

  // first extract the unassigned hits within precluster
  std::vector < GarlicExtendedHit* > unassignedHits = *(preCluster->getHits());
  for ( int ihh=unassignedHits.size()-1; ihh>=0; ihh-- ) {
    for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator clit = clusters.begin(); clit!=clusters.end(); clit++ ) {
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
    std::map < CalorimeterHit*, std::vector < GarlicExtendedHit* > > preassignedHits;
    // find distance to closest hit in core
    for ( int ihh=unassignedHits.size()-1; ihh>=0; ihh-- ) {
      float mindist=99999;
      CalorimeterHit* closestCore=0;
      for ( std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator clit = clusters.begin(); clit!=clusters.end(); clit++ ) {

        // distance from hit to core axis
        float distToCore = clit->second->getDistToClusterAxis( unassignedHits[ihh]->getPosition(), 0 ); // this assumes pointing photon
        if ( distToCore>maxDistCut ) continue; // do not consider this combination

        std::vector < GarlicExtendedHit* >* corehits = clit->second->getHits();
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
          std::vector < GarlicExtendedHit* > temp;
          temp.push_back( unassignedHits[ihh] );
          preassignedHits[closestCore]=temp;
        }
        unassignedHits.erase(unassignedHits.begin()+ihh);
      }
    }
    for ( std::map < CalorimeterHit*, std::vector < GarlicExtendedHit* > >::iterator preass = preassignedHits.begin(); preass!=preassignedHits.end(); preass++) {
      //      cout << "adding preassigned hits to core" << endl;
      clusters[preass->first]->addHits(preass->second);
    }
  } // iterations

  return clusters;
}


std::map < CalorimeterHit*, GarlicExtendedCluster* > GarlicClusterAlgos::getCores( GarlicExtendedCluster* preCluster, std::vector <CalorimeterHit*> seeds ) {
  std::map < CalorimeterHit*, GarlicExtendedCluster* > cores;
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
  std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator itt;
  std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator jtt;

  for (itt=cores.begin(); itt!=cores.end(); itt++) {
    std::vector<GarlicExtendedHit*>* ih = itt->second->getHits();
    if (ih->size()==0) continue;
    for (jtt=itt; jtt!=cores.end(); jtt++) {
      if (jtt==itt) continue;
      if (jtt==cores.end()) continue;
      for (int ii= (int) ih->size()-1; ii>=0; ii--) {
        std::vector<GarlicExtendedHit*>* jh = jtt->second->getHits();
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


GarlicExtendedCluster* GarlicClusterAlgos::BuildCore(CalorimeterHit* mySeed, GarlicExtendedCluster* preCluster, const float *clusterDir) {

  int nl_1         = GarlicAlgorithmParameters::Instance().GetCoreLayersSection1();
  int allowedGap_1 = GarlicAlgorithmParameters::Instance().GetCoreMaxHoleSection1();
  int allowedGap_2 = GarlicAlgorithmParameters::Instance().GetCoreMaxHoleSection2();

  GarlicExtendedHit seedHit(mySeed);

  std::vector <GarlicExtendedHit*> coreHits;

  std::vector <GarlicExtendedHit*> allPreclusterHits = *(preCluster->getHits());

  std::map <int, int> hitPseudoLayersBarrel;
  std::map <int, int> hitPseudoLayersEndcap;

  // define the distance window to add hits to the core
  //  float distCut = 1.5*GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];
  float distCut =
    GarlicAlgorithmParameters::Instance().GetCoreDistanceCut()*
    GarlicGeometryParameters::Instance().Get_padSizeEcal()[0];

  //  cout << "cell size, dist cut = " << GarlicGeometryParameters::Instance().Get_padSizeEcal()[0] << " " << distCut << endl;

  for (size_t ih=0; ih<allPreclusterHits.size(); ih++) {

    GarlicExtendedHit* thisHit = allPreclusterHits[ih];

    // project all hit onto plane defined by seed position and clusterDir
    const float* hitpos = thisHit->getPosition();
    float proj[3];
    GarlicGeometryHelpers::GetGeneralPointPlaneProjection(hitpos, mySeed->getPosition(), mySeed->getPosition(), clusterDir, proj);

    // distance from seed position to projected hit
    float dist(0);
    for (int i=0; i<3; i++) dist+= pow( proj[i] - mySeed->getPosition()[i], 2);
    dist=sqrt(dist);

    // if (dist>1000) {
    //   cout << "GarlicClusterAlgos::BuildCore : probable ERROR in projection measurement! too large distance " << dist << endl;
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

    for (int i=0; i< GarlicGeometryParameters::Instance().Get_nPseudoLayers(); i++) {
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

  GarlicExtendedCluster* coreClus = new GarlicExtendedCluster();
  coreClus->addHits( coreHits );
  coreClus->setReferenceDirection( clusterDir );

  return coreClus;
}



void GarlicClusterAlgos::mergeAllSatellitesAndElectrons(std::map < CalorimeterHit*, GarlicExtendedCluster* > & clusters ,
                                                       vector < std::pair < GarlicExtendedTrack*, GarlicExtendedCluster* > > & allElectrons ) {

  // this one tries to merge photon and electron clusters simulataneously

  if (_verbose) {
    cout << "hello from mergeAllSatellitesAndElectrons!" << endl;
    cout << clusters.size() << " clusters, " << allElectrons.size() << " electrons" << endl;
  }

  // so we can easily see whcih clusters are electrons
  std::vector < GarlicExtendedCluster* > elecCls;

  // order all clusters and electrons by energy
  std::map < float, GarlicExtendedCluster* > enClusters; // energy ordered list of clusters

  std::map < CalorimeterHit*, GarlicExtendedCluster* >::iterator ittt;

  for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
    enClusters[ -(ittt->second->getEnergy()) ] = ittt->second;   // *-1 to avoid using reverse iterator...
  }
  for (size_t i=0; i<allElectrons.size(); i++) {
    enClusters[ -allElectrons[i].second->getEnergy() ] = allElectrons[i].second;
    elecCls.push_back(allElectrons[i].second);
  }


  std::map < float, GarlicExtendedCluster* >::iterator itt, jtt;

  bool changed = true;
  while (changed) {
    changed = false;
    for ( itt = enClusters.begin(); itt!=enClusters.end(); itt++) {
      if ( itt==enClusters.begin() ) continue;
      GarlicExtendedCluster* cl_lo = itt->second;
      if (!cl_lo) continue;

      bool lo_is_electron = find ( elecCls.begin(), elecCls.end(), cl_lo ) != elecCls.end();

      if ( _verbose ) cout << "trying to combine lower energy cluster with energy " << cl_lo->getEnergy() << " electron? " << lo_is_electron << endl;

      std::map < float, GarlicExtendedCluster* > allScores;

      // look at all higher energy clusters
      for ( jtt = enClusters.begin(); jtt!=enClusters.end(); jtt++) {
        if ( jtt==itt ) break;

        GarlicExtendedCluster* cl_hi = jtt->second;
        if (!cl_hi) continue;

        bool hi_is_electron = find ( elecCls.begin(), elecCls.end(), cl_hi ) != elecCls.end();

        if ( _verbose ) cout << "    considering higher energy cluster " << cl_hi->getEnergy() << " electron? " << hi_is_electron << endl;

        if ( 
	    (lo_is_electron && hi_is_electron) || // don't merge 2 electrons
	    ! mergeCandidate( cl_hi, cl_lo ) ) {  // clusters should not be merged
          if ( _verbose ) cout << "NOT A CANDIDATE" << endl;
          continue;
        }

        // original score of higher energy cluster

        float orig_score(0);
        // is this cluster a photon or electron?
        if ( hi_is_electron ) {
          orig_score = 100*cl_hi->getElectronCutSel( ); // + // overall selection
          // cl_hi->getElectronCutSel(GarlicExtendedCluster::TRANS ) +
          // cl_hi->getElectronCutSel(GarlicExtendedCluster::LONGT ) +
          // cl_hi->getElectronCutSel(GarlicExtendedCluster::MISCL ) +
          // cl_hi->getElectronCutSel(GarlicExtendedCluster::POINT ) ;
        } else {
          orig_score = 100*cl_hi->getPhotonCutSel( ); // +
          // cl_hi->getPhotonCutSel(GarlicExtendedCluster::TRANS ) +
          // cl_hi->getPhotonCutSel(GarlicExtendedCluster::LONGT ) +
          // cl_hi->getPhotonCutSel(GarlicExtendedCluster::MISCL ) +
          // cl_hi->getPhotonCutSel(GarlicExtendedCluster::POINT ) ;
        }

        // distance between cluster axes
        int directionassump = hi_is_electron ? 2 : 0 ; // track direction for electron, IP direction for photon
        float dist = cl_hi->getDistToClusterAxis( cl_lo->getCentreOfGravity(), directionassump );

        // now make a combined cluster
        GarlicExtendedCluster* combo = new GarlicExtendedCluster();
        combo->addHits( *(cl_hi->getHits()) );
        combo->addHits( *(cl_lo->getHits()) );
        combo->setAssociatedTrack(cl_hi->getAssociatedTrack());

        float new_score(0);
        // is this cluster a photon or electron?
        if ( hi_is_electron ) {
          new_score = 100*cl_hi->getElectronCutSel( ) ; // overall selection
          //cl_hi->getElectronCutSel(GarlicExtendedCluster::TRANS ) +
          //cl_hi->getElectronCutSel(GarlicExtendedCluster::LONGT  ) +
          //cl_hi->getElectronCutSel(GarlicExtendedCluster::MISCL ) +
          //cl_hi->getElectronCutSel(GarlicExtendedCluster::POINT ) ;
        } else {
          new_score = 100*cl_hi->getPhotonCutSel( ) ;
          //cl_hi->getPhotonCutSel(GarlicExtendedCluster::TRANS ) +
          //cl_hi->getPhotonCutSel(GarlicExtendedCluster::LONGT ) +
          //cl_hi->getPhotonCutSel(GarlicExtendedCluster::MISCL ) +
          //cl_hi->getPhotonCutSel(GarlicExtendedCluster::POINT ) ;
        }


        float improvement = new_score - orig_score;
        // if combinations give same improvement, decide by distance
        // (however, keep it above 100, 200 etc so we can easily see if it's loose, tight)
        improvement += (10.-dist/100.);

        if ( _verbose ) {
          cout << "    orig, new scores = " << orig_score << " " << new_score << " : " << improvement << endl;
          cout << "    distance = " << dist << endl;
        }


        if ( improvement > -50 ) { // doesn't turn tight into loose, loose into rejected
          allScores[improvement] = cl_hi;
        }

      } // second loop over higher energy clusters


      if ( allScores.size()==0 ) continue;

      std::map < float, GarlicExtendedCluster* >::reverse_iterator kj;
      std::map < float, GarlicExtendedCluster* >::iterator kjj;

      if ( _verbose ) {
        cout << "all scores for this cluster: parent energy " << cl_lo->getEnergy() << " " << allScores.size() << endl;
        for ( kj=allScores.rbegin(); kj!=allScores.rend(); kj++) {
          cout << "   child energy, score     " << kj->second->getEnergy() << " " << kj->first << endl;
        }
      }

      GarlicExtendedCluster* bestCluster = allScores.rbegin()->second;
      float bestScore = allScores.rbegin()->first;

      if ( bestCluster ) { // then combine
        changed=true;
        bestCluster->addHits( *(cl_lo->getHits()) );

        if ( _verbose ) cout << "combining! best score = " << bestScore << " new energy = " << bestCluster->getEnergy() << endl;

        for ( ittt = clusters.begin(); ittt!=clusters.end(); ittt++) {
          if ( ittt->second == cl_lo ) {
            clusters.erase(ittt);// delete entry from orig map
          }
        }
        enClusters.erase(itt);  // delete entry from the energy-ordered map
        delete cl_lo;  // delete the cluster object itself
      }

      if (changed) break;

    } // first loop over clusters

  } // while changed



}



bool GarlicClusterAlgos::mergeCandidate( GarlicExtendedCluster* primary, GarlicExtendedCluster* secondary ) {
  // chg+neu: primary should be charged with neutral secondary partner
  // neu+neu: primary should have higher energy

  if ( secondary->getAssociatedTrack() ) {
    if (  primary->getAssociatedTrack() ) {
      cout << "asking GarlicClusterAlgos::mergeCandidate to make judgement about two charged clusters...it's not designed for this!" << endl;
      return false;
    } else { // switch to put charged in primary position
      GarlicExtendedCluster* temp = primary;
      primary=secondary;
      secondary=temp;
    }
  } else if ( !primary->getAssociatedTrack() ) { // 2 neutral clusters
    if ( primary->getEnergy() < secondary->getEnergy() ) { // make sure primary has higher energy
      GarlicExtendedCluster* temp = primary;
      primary=secondary;
      secondary=temp;
    }
  }

  float maxClEn = std::max( primary->getEnergy(), secondary->getEnergy() );

  float ratioCut                   = GarlicAlgorithmParameters::Instance().GetMergeRatioCut              ();
  float energy_dist_factor         = GarlicAlgorithmParameters::Instance().GetMergeEnergyDistFactor      ();
  float absoluteLargestDist        = GarlicAlgorithmParameters::Instance().GetMergeAbsoluteLargestDist   ();
  float mass_limit                 = GarlicAlgorithmParameters::Instance().GetMergePi0MassLimit          ();
  float max_mass_imbalance_for_pi0 = GarlicAlgorithmParameters::Instance().GetMergePi0MaxEnergyImbalance ();

  // energy-dependent max distance between clusters
  float min_max_merge_dist  = GarlicAlgorithmParameters::Instance().GetMergeMaxDistAtLowEn();
  float edap_max_merge_dist = GarlicAlgorithmParameters::Instance().GetMergeMaxDistEnDep();
  float max_merge_dist_rm = std::max ( float(min_max_merge_dist), float(min_max_merge_dist + edap_max_merge_dist*log10( maxClEn ) ));
  float max_merge_dist_mm = max_merge_dist_rm*GarlicAlgorithmParameters::Instance().GetMoliereRadius();

  // first a cut on the simple distance between centres of gravity
  float simpleCCdist(0);
  for (int i=0; i<3; i++)
    simpleCCdist+=pow(primary->getCentreOfGravity()[i] - secondary->getCentreOfGravity()[i],2);
  simpleCCdist=sqrt(simpleCCdist);
  if ( simpleCCdist>absoluteLargestDist ) return false;

  // then the distance from cluster axis (energy-dependent cut)
  float cogDist =  primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 0 ); // ip pointing
  if ( primary->getAssociatedTrack() ) cogDist =  std::min( cogDist, primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 1 )); // track direction
  cogDist =  std::min( cogDist, primary->getDistToClusterAxis(  secondary->getCentreOfGravity(), 2 )); // cluster axis

  float energyRatio = secondary->getEnergy()/primary->getEnergy();

  //cout << " dist " << cogDist << " energyRatio " << energyRatio << endl;

  bool ok_enRatio = energyRatio<ratioCut;                               // don't merge similar-energy clusters
  bool ok_enRatioDist = energyRatio<energy_dist_factor/pow(cogDist,2);  // energy fraction-dependent distance cut
  bool ok_Dist = cogDist<max_merge_dist_mm; // distance

  if ( _verbose ) {
    if ( !ok_Dist ) {
      cout << "failed Dist cut " << cogDist << " > " << max_merge_dist_mm << endl;
    }
  }

  // this is to avoid merging pi0s, based on inv mass
  bool ok_Pi0(true);
  if ( ! primary->getAssociatedTrack() && ! secondary->getAssociatedTrack() ) {
    std::pair < float, float > comboMassErr = primary->getCombinedInvariantMass(secondary);
    if ( comboMassErr.first+2*comboMassErr.second>mass_limit &&
         primary->getEnergy()/secondary->getEnergy() < max_mass_imbalance_for_pi0 ) {
      ok_Pi0=false;
      if ( _verbose ) {
        cout << "refusing to combine clusters whose mass is close to pizero mass! " <<
          comboMassErr.first << " " << comboMassErr.second << " " <<
          comboMassErr.first+2*comboMassErr.second << " " << mass_limit  << endl;
        cout << " energy imbalance = " << primary->getEnergy()/secondary->getEnergy() << " limit " << max_mass_imbalance_for_pi0 << endl;
      }
    }
  }

  if (_verbose)
    cout << "  merge? ok_Pi0, ok_enRatio, ok_enRatioDist, ok_Dist = " << ok_Pi0 <<  " " << ok_enRatio << " " << ok_enRatioDist << " " << ok_Dist << endl;

  /*
  // check cluster start poistions

  // for the moment don't use shower start information...didn't seem too useful
  // maybe revisit later....djeans, jan2015

  float* posP = primary->getClusterStartPosition();
  float* posS = secondary->getClusterStartPosition();

  int pLayP = primary->getClusterStartPlayer();
  int pLayS = secondary->getClusterStartPlayer();

  int dPlS = pLayS - pLayP;

  float magP(0);
  float magS(0);
  float dot(0);
  for (int i=0; i<3; i++) {
    magP+=pow(posP[i],2);
    magS+=pow(posS[i],2);
    dot+=posP[i]*posS[i];
  }
  magP=sqrt(magP);
  magS=sqrt(magS);
  dot/=magP*magS;
  float angle = acos(dot);
  float dT = sin(angle) * std::min( magS, magP );
  float dR = magS - magP;

  float aveLayerThickness =
    0.5*( GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[3]-GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[1] );

  if (_verbose) {
    cout << "merge candidate....start position: ";
    for (int i=0; i<3; i++) cout << primary->getClusterStartPosition()[i] << " " ;
    cout << " : ";
    for (int i=0; i<3; i++) cout << secondary->getClusterStartPosition()[i] << " " ;
    cout << " : " << magP << " " << magS << " " << dot << endl;
    cout << " dT, dR = " << dT << " " << dR << " cc dist = " << simpleCCdist << endl;
  }

  // we expect the smaller cluster to be later than the primary one


  bool ok_showerStart =
    dPlS > 3 ||
    fabs(dR) > 3.*aveLayerThickness ||
    dT < fabs(dR); // 45 degree angle

  if ( _verbose && ok_Pi0 && ok_enRatio && ok_enRatioDist && ok_Dist ) {
    cout << "combining clusters; dT, dR = " << dT << " " << dR << " d_p_layer = " << dPlS <<
      " ( sinpleCC " << simpleCCdist << " enRatio " << energyRatio << " " <<
      primary->getEnergy() << " " << secondary->getEnergy() << " : shower start? " << ok_showerStart << " " << dPlS*aveLayerThickness << " " << dT << endl;
  }
  */


  //  return ok_showerStart &&
  return ok_Pi0 && ok_enRatio && ok_enRatioDist && ok_Dist;

}

