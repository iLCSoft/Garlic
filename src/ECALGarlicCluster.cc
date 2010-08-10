#include "ECALGarlicCluster.hh"

#include <iostream>
using std::cout;
using std::endl;

#include "TVector3.h"
#include "TF1.h"
#include "ECALGarlicClusterHelpers.hh"
#include "ECALGarlicAlgorithmParameters.hh"
#include "ECALGarlicGeometryHelpers.hh"


void ECALGarlicCluster::Cluster(ExtendedCluster &preCluster, vector<vec3> &possibleSeeds, vector<ExtendedCluster* > &clusters,vector<ExtendedTrack*> tracks )
{

  vector<vec3>::iterator seedIt;
  for(seedIt=possibleSeeds.begin();seedIt!=possibleSeeds.end();seedIt++) {
    vec3 seed=*seedIt;

    // 1.) to find nearest hit to seed : REMOVED

    if(_algoParams->GetDebug()>2)
      cout << "The seed is at (" << seed.x << ", " << seed.y << ", " << seed.z << ")" << endl;

    vector<ExtendedHit* > myCluster;

    // 2.) Build central cluster core
    vec3 clusterDirection;
    clusterDirection.x=seed.x;
    clusterDirection.y=seed.y;
    clusterDirection.z=seed.z;
    AddCore(myCluster, &seed, preCluster, &clusterDirection);

    // 3.) the core is built so start clustering neighbours
    if(myCluster.size()==0) // no core found, check next seed
      continue;
    if(_algoParams->GetDebug()>2) cout  << myCluster.size() << " hits in core" << endl;

    BuildClusterFromNeighbours(myCluster, &clusterDirection, preCluster);
    if(_algoParams->GetDebug()>2) cout << "Clusterd " << myCluster.size() << " hits" << endl;

    if(int(myCluster.size())<_algoParams->GetNHitsMin()) {
      _clusterHelper->FreeHits(myCluster);
      if(_algoParams->GetDebug()>1)
	cout << "Not enough hits in cluster!" << endl;
      continue;
    }
    if(_algoParams->GetDebug()>1)
      cout << "Finished clustering of " << myCluster.size() << " hits, now checking cluster criteria" << endl;
    
    ExtendedCluster *aRealCluster = new ExtendedCluster;
    aRealCluster->hitVec=myCluster;
    aRealCluster->seededFrom=seed;
    aRealCluster->PreCluster=&preCluster;
    aRealCluster->location=preCluster.location;
    double En_a=0;
    int nHits_a=aRealCluster->hitVec.size();
    for(int hit_i=0;hit_i<nHits_a;hit_i++) {
      //int layer = (((aRealCluster->hitVec[hit_i])->hit)->getCellID0() >> 24)+1;
      //if(((aRealCluster->hitVec[hit_i])->preShower!=0))
      En_a+=(((aRealCluster->hitVec[hit_i])->hit)->getEnergy());
      ((aRealCluster->hitVec[hit_i])->cluster) = (aRealCluster);
    }
    aRealCluster->rawEn=En_a;
    
    bool isReal = _clusterHelper->CheckClusterCriteria(aRealCluster,tracks );
    if(isReal) {   // build an ExtendedCluster from it
      
      clusters.push_back(aRealCluster);
      sort(clusters.begin(),clusters.end(),ExtendedCluster::higherRawEnergy);
    }
    else {
      _clusterHelper->FreeHits(aRealCluster);
      delete aRealCluster;      
      continue;
    }
  }
  if(_algoParams->GetDebug()>1) {
    if(clusters.size()>0) {
      cout << "Found " << clusters.size() << " clusters from " << possibleSeeds.size() << " seeds" << endl;
      for(seedIt=possibleSeeds.begin();seedIt!=possibleSeeds.end();seedIt++) {
	vec3 seed=*seedIt;
	cout << "at " << seed.x << ", " << seed.y << ", " << seed.z << ")" << endl;
      }
    }
  }
}

void ECALGarlicCluster::AddCore(vector<ExtendedHit* > &aCluster, vec3 *mySeed, ExtendedCluster &preCluster, vec3 *clusterDir) {  
  int _n_ly_10X0=16;
  int _n_ly_hole_cut=20;// This is actually not 10X0 anymore but rather the second stack
  bool no2ndHit=true;
  bool hole_after_cut=false;
  bool continueAfterGapJump=false;
  int gapJumpTo=0;
  vec3 afterJump;
  TVector3 dir(mySeed->x,mySeed->y,mySeed->z);
  double Theta = dir.Theta();

  for(int lay=0;lay<_geomParams->Get_nPseudoLayers();lay++) {

    if(lay>(_algoParams->GetNLayersForSeeding()+1) && no2ndHit) {
      if(_algoParams->GetDebug()>2) cout << "No 2nd hit until layer " << (_algoParams->GetNLayersForSeeding()+1) << ", cancelling search for core" << endl;
      _clusterHelper->FreeHits(aCluster);
      break;
    }
    if(lay>_n_ly_10X0 && aCluster.size()<((unsigned int) _algoParams->GetMinHits10X0())) { // have minimum 4 hits after 10X0
      if(_algoParams->GetDebug()>2) cout << "Not enough hits in 10X0, cancelling search for core" << endl;
      _clusterHelper->FreeHits(aCluster);
      break;
    }

    if(_algoParams->GetDebug()>2)
      cout << "Following cluster from: " << mySeed->x << ", " << mySeed->y << ", " << mySeed->z << "(seed) for " << lay << " layers" << endl;

    bool gapHoleVeto=false;
    gapJumpTo=0;
    vec3 next_search_spot=_clusterHelper->FollowClusterDirection(mySeed,lay,clusterDir,preCluster.location,gapHoleVeto,gapJumpTo);
    if(gapHoleVeto) {
      if (_algoParams->GetDebug()>0) cout << "Couldn't bridge gap with layer " << lay << endl;
      continue;
    }
    if(gapJumpTo>0) {
      continueAfterGapJump=true;
      afterJump=next_search_spot;
      if(_algoParams->GetDebug()>2) cout << "Will continue after gap jump..." << endl;
      break;
    } else {
      bool stop_it=false;
      for(int n_i=0;n_i<10;n_i++) {
	ExtendedHit *nearestHit=NULL;

	_geomHelper->GetNearestHitInPseudoLayer(lay,&next_search_spot,&preCluster,nearestHit);
	if(nearestHit==0) {
	  if(lay>_n_ly_hole_cut && n_i ==0) {
	    if(hole_after_cut==true) {
	      if(_algoParams->GetDebug()>2) 
		cout << "Second hole in core after 10 X0 in layer " << lay << ": core finished" << endl;
	      stop_it=true;
	      break;
	    }
	    hole_after_cut=true;
	  }
	  if(_algoParams->GetDebug()>2) cout << "Braking at n_i = " << n_i << " in layer " << lay << endl;
	  break;  
	}
	vec3 nearestPos;
	nearestPos.x=nearestHit->hit->getPosition()[0];
	nearestPos.y=nearestHit->hit->getPosition()[1];
	nearestPos.z=nearestHit->hit->getPosition()[2];
	double dist=_geomHelper->Get2dProjDistance(&nearestPos,mySeed); 
	if(_algoParams->GetDebug()>2) cout << "Distance is: " << dist << endl;

	double modifier =1;
	if(preCluster.location==2)
	  modifier=fabs(cos(Theta));
	else
	  modifier=sin(Theta);

	if(dist<(1.*(_geomParams->Get_padSizeEcal()[1]))/modifier) { 		    
	  nearestHit->clusterHitVec=&aCluster;
	  aCluster.push_back(nearestHit);
	  if(aCluster.size()>1)
	    no2ndHit=false;
	  if(_algoParams->GetDebug()>2) cout << "Added hit: " << nearestHit->hit->getCellID0() << 
			 " at (" << nearestPos.x << ", " << nearestPos.y << ", " << nearestPos.z << "),in layer " << lay << endl;
	  if(lay>_n_ly_hole_cut) {
	    hole_after_cut=false;
	  }
	} else {
	  if(n_i>0)
	    break;
	  if(lay>_n_ly_hole_cut) { // this was after 10 X0... no after hole cut
	    if(hole_after_cut==true) {
	      if(_algoParams->GetDebug()>2) 
		cout << "Second hole in core after 10 X0 in layer " << lay << ": core finished" << dist << endl;
	      stop_it=true;
	      break;
	    }
	    else {
	      if(lay==_geomParams->Get_nPseudoLayers()) {
		if(_algoParams->GetDebug()>2) 
		  cout << "End of ECAL in layer " << lay << ": core finished" << dist << endl;
		stop_it=true;
		break;
	      }
	      int nextLayerIs=lay+1;
	      hole_after_cut=true;
	      bool gapHoleVeto=false; // TODO: what if there is a hole in front of the cable gap?
	      int gap_jump_to=0;
	      vec3 after_next_search_spot=_clusterHelper->FollowClusterDirection(&next_search_spot,1,clusterDir,preCluster.location,gapHoleVeto,gap_jump_to);
	      ExtendedHit *next_layer_nearestHit=NULL;
	      _geomHelper->GetNearestHitInPseudoLayer(nextLayerIs,&next_search_spot,&preCluster,next_layer_nearestHit); // ...start looking for holes of two layers
	      if(next_layer_nearestHit==0) {
		if(_algoParams->GetDebug()>2) 
		  cout << "No nearest hit in next layer : core finished"  << endl;
		stop_it=true;
		break;
	      }
	      vec3 next_layer_nearestPos;
	      next_layer_nearestPos.x=next_layer_nearestHit->hit->getPosition()[0];
	      next_layer_nearestPos.y=next_layer_nearestHit->hit->getPosition()[1];
	      next_layer_nearestPos.z=next_layer_nearestHit->hit->getPosition()[2];
	      double dist_nl=_geomHelper->Get2dProjDistance(mySeed,&after_next_search_spot);

	      double modifier =1;
	      if(preCluster.location==2)
		modifier=fabs(cos(Theta));
	      else
		modifier=sin(Theta);

	      if(dist_nl>(1.*(_geomParams->Get_padSizeEcal()[0]))/modifier) {
		if(_algoParams->GetDebug()>2) 
		  cout << "Nearest hit in next layer is not close enough : core finished"  << endl;
		stop_it=true;
		break;
	      }
	    }
	  }
	}
      }
      // get out of layer loop ...
      if(stop_it==true)
	break;
    }
  }
  
  if(continueAfterGapJump==true) {
    if(_algoParams->GetDebug()>2) 
      cout << "...at " << afterJump.x << ", " << afterJump.y << ", " << afterJump.z << endl;
    //now continue in the endcap: since the shower seems to be spread a lot wider, we are more generous with the distance cuts
    for(int lay=0;lay<_geomParams->Get_nPseudoLayers();lay++) {
      bool gapHoleVeto=false;
      int gap_jump_to=0;
      vec3 next_search_spot=_clusterHelper->FollowClusterDirection(&afterJump,lay,clusterDir,preCluster.location,gapHoleVeto,gap_jump_to);
      if(gapHoleVeto==true) {
	cout << "Why did you get here? You are already in the endcap!" << endl;
	continue;
      }
      bool stop_it=false;
      for(int n_i=0;n_i<9;n_i++) {
	ExtendedHit *nearestHit=NULL;
	_geomHelper->GetNearestHitInPseudoLayer(lay,&next_search_spot,&preCluster,nearestHit);
	if(nearestHit==0) {
	  if(lay>_n_ly_hole_cut && n_i ==0) {
	    if(hole_after_cut==true) {
	      if(_algoParams->GetDebug()>2) 
		cout << "Second hole in core after 10 X0 in Endcap layer " << lay << ": core finished" << endl;
	      stop_it=true;
	      break;
	    }
	    hole_after_cut=true;
	  }
	  if(_algoParams->GetDebug()>2)
	    cout << "Braking at n_i = " << n_i << " in Endcap layer " << lay << endl;
	  break;  
	}
	vec3 nearestPos;
	nearestPos.x=nearestHit->hit->getPosition()[0];
	nearestPos.y=nearestHit->hit->getPosition()[1];
	nearestPos.z=nearestHit->hit->getPosition()[2];
	double dist=_geomHelper->Get2dProjDistance(&nearestPos,&afterJump); 
	if(_algoParams->GetDebug()>2)
	  cout << "Distance is: " << dist << endl;
	if(dist<sqrt(8.)/fabs(cos(Theta))*(_geomParams->Get_padSizeEcal()[0])) { 		    
	  nearestHit->clusterHitVec=&aCluster;
	  aCluster.push_back(nearestHit);
	  if(_algoParams->GetDebug()>2)
	    cout << "Added hit: " << nearestHit->hit->getCellID0() << " at (" << nearestPos.x << ", " << nearestPos.y << ", " << nearestPos.z << "),in Endcap layer " << lay << endl;
	  if(lay>_n_ly_hole_cut) {
	    hole_after_cut=false;
	  }
	}
	else {
	  if(n_i>0)
	    break;
	  if(lay>_n_ly_hole_cut) { // This was after 10 X0... no after hole cut
	    if(hole_after_cut==true) {
	      if(_algoParams->GetDebug()>2) 
		cout << "Second hole in core after 10 X0 in Ecdcap layer " << lay << ": core finished" << dist << endl;
	      stop_it=true;
	      break;
	    }
	    else {
	      if(lay==_geomParams->Get_nPseudoLayers()) {
		if(_algoParams->GetDebug()>2) 
		  cout << "End of ECAL in Endcap layer " << lay << ": core finished" << dist << endl;
		stop_it=true;
		break;
	      }
	      int nextLayerIs=lay+1;
	      hole_after_cut=true;
	      bool gapHoleVeto=false; // TODO: what if there is a hole in front of the cable gap?
	      int gap_jump_to=0;
	      vec3 after_next_search_spot=_clusterHelper->FollowClusterDirection(&next_search_spot,1,clusterDir,preCluster.location,gapHoleVeto,gap_jump_to);
	      ExtendedHit *next_layer_nearestHit=NULL;
	      _geomHelper->GetNearestHitInPseudoLayer(nextLayerIs,&next_search_spot,&preCluster,next_layer_nearestHit); // ...start looking for holes of two layers
	      if(next_layer_nearestHit==0) {
		if(_algoParams->GetDebug()>2) 
		  cout << "No nearest hit in next Endcap layer : core finished"  << endl;
		stop_it=true;
		break;
	      }
	      vec3 next_layer_nearestPos;
	      next_layer_nearestPos.x=next_layer_nearestHit->hit->getPosition()[0];
	      next_layer_nearestPos.y=next_layer_nearestHit->hit->getPosition()[1];
	      next_layer_nearestPos.z=next_layer_nearestHit->hit->getPosition()[2];
	      double dist_nl=_geomHelper->Get2dProjDistance(&afterJump,&after_next_search_spot);
	      if(dist_nl>sqrt(5.)/fabs(cos(Theta))*(_geomParams->Get_padSizeEcal()[lay])) {
		if(_algoParams->GetDebug()>2) 
		  cout << "Nearest hit in next Endcap layer is not close enough : core finished"  << endl;
		stop_it=true;
		break;
	      }
	    }
	  }
	}
      }
      // get out of layer loop ...
      if(stop_it==true)
	break;
    }
  }
}

void ECALGarlicCluster::BuildClusterFromNeighbours(vector<ExtendedHit* > &myCluster, vec3 *clusterDir, ExtendedCluster &preCluster) 
{
  TVector3 dir(clusterDir->x,clusterDir->y,clusterDir->z);
  double Theta = dir.Theta();
  for(int iteration=0;iteration<_algoParams->GetNIterations(); iteration++) {
    //cout << " Iteration: " << iteration << endl;
    // create a copy from the cluster to restrain added hits in 1 iteration
    vector<ExtendedHit* > myCluster_copy=myCluster;
    if(_algoParams->GetDebug()>2) cout << "Cluster copied... now working on copy!" << endl;
    for(int layer_i=(myCluster_copy[0]->pseudoLayer);layer_i<_geomParams->Get_nPseudoLayers();layer_i++) {
      int NHitsInCopy=myCluster_copy.size();
      for( int hit_i = 0; hit_i<NHitsInCopy; hit_i++ ) {
	ExtendedHit *a_copy_hit=dynamic_cast<ExtendedHit *>(myCluster_copy[hit_i]);
	if(a_copy_hit) {
	  if(a_copy_hit->pseudoLayer==layer_i) {
	    //insert the neighbours
	    vector<ExtendedHit* > freeNeighbours;
	    _clusterHelper->GetFreeNeighbours(a_copy_hit, freeNeighbours, &preCluster, iteration, Theta);
	    int NNeighbours=freeNeighbours.size();
	    if(NNeighbours>0) {
	      for(int neigh_i=0;neigh_i<NNeighbours;neigh_i++) {
		ExtendedHit *a_neighbour=dynamic_cast<ExtendedHit *>(freeNeighbours[neigh_i]);
		if(a_neighbour) {
		  if(_algoParams->GetDebug()>2) cout << "Added neighbour hit "  << a_neighbour->hit->getCellID0() << endl;
 		  a_neighbour->clusterHitVec=&myCluster;
		  myCluster.push_back(a_neighbour);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}



void ECALGarlicCluster::MergeSatellites(vector<ExtendedCluster* > *clusters)
{
  // clusters are ordered by energy. 
  // a cluster is considered a satellite if : 1.) invariant mass of the two clusters combinded is lower than pi0 mass, 2.) (minimum) distance from possible main cluster is smaller than user value, 3.) Fits in a given CosTh:EnRatio value
  int NClusters=clusters->size();
  if(NClusters>1) {
    //map<ExtendedCluster* ,vector<ExtendedCluster* >* > main_clusters;
    //vector<ExtendedCluster* > *satellites = new vector<ExtendedCluster* >();
    vector<ExtendedCluster* > main_clusters;
    vector<ExtendedCluster* > satellites;
    vector<ExtendedCluster* >::iterator sat_it;
    vector<ExtendedCluster* >::iterator main_it;
    vector<ExtendedCluster* >::iterator clus_it;

    // initialise functions
    TF1 *para0 = new TF1("para0","[0]+[1]*exp([2]*x)",0,400);
    para0->SetParameter(0,9.99855e-01);
    para0->SetParameter(1,-3.67450e-03);
    para0->SetParameter(2,-1.07102e-01);

    TF1 *para1 = new TF1("para1","[0]+[1]*exp([2]*x)",0,400);
    para1->SetParameter(0,-1.95419);
    para1->SetParameter(1,1.29176e-01);
    para1->SetParameter(2,-6.81975e-02);

    for(int cl_i=0;cl_i<NClusters;cl_i++) {
      ExtendedCluster *a_ext_cluster = dynamic_cast<ExtendedCluster*>((*clusters)[cl_i]);
      ClusterParameters *clus_par = a_ext_cluster->Parameters;
      double clus_en = clus_par->E_GeV;
      if(clus_en>5)
	main_clusters.push_back((*clusters)[cl_i]);
      if(clus_en<_algoParams->GetMaxSatEn() && clus_en>_algoParams->GetMinSatEn()) {
	satellites.push_back((*clusters)[cl_i]);
      }	
    }
    if(main_clusters.size()>0 && satellites.size()>0) {
      for(sat_it=satellites.begin();sat_it!=satellites.end();sat_it++) {
	bool is_satellite=false;
	// find nearest main cluster
	double min_dist = 99999;
	double max_en = 0;
	double distances[3];
	ExtendedCluster *main=0;
	for(main_it=main_clusters.begin();main_it!=main_clusters.end();main_it++) {
	  _geomHelper->GetDistancesBetweenClusters(*sat_it, *main_it, distances);
	  double clus_en = ((*main_it)->Parameters)->E_GeV;
	  if(distances[0]<min_dist) {
	    min_dist = distances[0];
	    max_en = clus_en;
	    main = (*main_it);
	  }
	  if(distances[0]==min_dist)
	    if(clus_en>max_en)
	      {
		min_dist = distances[0];
		max_en = clus_en;
		main = (*main_it);
	      }
	}
	if(main) {
	  double main_en = ((main)->Parameters)->E_GeV;
	  double max_sat_dist = 4*_geomParams->Get_padSizeEcal()[0];
	  if(main_en>15)
	    max_sat_dist =(6*_geomParams->Get_padSizeEcal()[0]) ;
	  if(min_dist< max_sat_dist) {
	    // set limits according to main energy
	    double cosTh_lim = 0.86;
	    double enRatio_lim = 0.08;
	    if(main_en>8) {
	      enRatio_lim = 0.06;
	      cosTh_lim = 0.9;
	    }
	    if(main_en>10) {
	      enRatio_lim = 0.05;
	      cosTh_lim = 0.9;
	    }
	    if(main_en>15) {
	      enRatio_lim = 0.035;
	      cosTh_lim = 0.92;
	    }
	    if(main_en>35) {
	      enRatio_lim = 0.03;
	      cosTh_lim = 0.93;
	    }
	    
	    double sat_en = ((*sat_it)->Parameters)->E_GeV;
	    // initialise functions
	    double par_0 = para0->Eval(main_en);
	    double par_1 = para1->Eval(main_en);
	    TF1 *lin = new TF1("lin","pol1",0,1);
	    lin->SetParameter(0,par_0);
	    lin->SetParameter(1,par_1);
	    
	    double enRatio = sat_en/main_en;
	    double cosTh = fabs(main_en-sat_en)/(main_en+sat_en);
	    //cout << "EnRatio: "  << enRatio << " , CosTh: " << cosTh << endl;
	    if(cosTh>cosTh_lim && enRatio<enRatio_lim) {
	      double ref_val = lin->Eval(enRatio);
	      //cout << " RefVal: " << ref_val << endl;
	      if(fabs(cosTh-ref_val) < 0.002) {
		is_satellite = true;
		//cout << "Satellite!" << endl;
	      }
	    }
	    if(is_satellite) {
	      //transefr hits
	      vector<ExtendedHit*> *main_hitVec = &(main->hitVec);
	      vector<ExtendedHit*> *sat_hitVec = &((*sat_it)->hitVec);
	      int NSatHits = sat_hitVec->size();
	      for(int h_i=0;h_i<NSatHits;h_i++) {
		ExtendedHit *myHit = dynamic_cast<ExtendedHit* > ((*sat_hitVec)[h_i]);
		myHit->cluster = main;
		myHit->clusterHitVec = main_hitVec;
		main_hitVec->push_back(myHit);
		//cout << "Transfered hit" << endl;
	      }
	      //delete satellite from cluster vector
	      for(clus_it=clusters->begin();clus_it!=clusters->end();clus_it++) {
		if(*clus_it==*sat_it) {
		  clusters->erase(clus_it);
		  //cout << "Erased satellite" << endl;
		  break;
		}
	      }
	    }
	    delete lin;
	  }
	}
      }
    }
    delete para0;
    delete para1;
  }
}

