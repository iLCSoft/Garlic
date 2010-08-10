/*********************************************************************
 * ECALPreClustering
 *
 * A Marlin processor for simple grouping of ECAL hits to obtain
 * clusters from which a ROI for the GARLIC processor is deduced.
 *
 * 26 November 2007  Marcel Reinhard (LLR)
 *********************************************************************/
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <marlin/Global.h>
#include <marlin/ProcessorMgr.h>

#include <ECALPreClustering.hh>

using namespace marlin;
using namespace lcio;
using namespace std;

ECALPreClustering anECALPreClustering;



ECALPreClustering::ECALPreClustering()
  : Processor("ECALPreClustering")
{
  // processor description
  _description = "Preclusters ECAL hits to define ROIs for Garlic";


  string ecalEndcapHitCollName = "ECALEndcap";
  registerProcessorParameter("EcalEndcapHitCollection",
			     "CalorimeterHit collection name of ECAL Endcap",
			     _ecalEndcapHitCollectionName,
			     ecalEndcapHitCollName);

  string ecalBarrelHitCollName = "ECALBarrel";
  registerProcessorParameter("EcalBarrelHitCollection",
			     "CalorimeterHit collection name of ECAL Barrel",
			     _ecalBarrelHitCollectionName,
			     ecalBarrelHitCollName);

  string ecalOtherHitCollName = "ECALOther";
  registerProcessorParameter("EcalOtherHitCollection",
			     "CalorimeterHit collection name of ECAL Ring",
			     _ecalOtherHitCollectionName,
			     ecalOtherHitCollName);

  string ecalPreClusterCollName = "ECAL_PreClusters";
  registerProcessorParameter("EcalPreClusterCollection",
			     "Collecion name of the ECAL preClusrters",
			     _ecalPreClusterCollectionName,
			     ecalPreClusterCollName);


  registerProcessorParameter("MinimumHits",
			     "Minimum number of hits to accept a PreCluster",
			     _minHits,
			     int(5));

  registerProcessorParameter("DistanceCut",
			     "Maximum distance of two hits grouped to the same PreCluster",
			     _distanceCut,
			     float(40.0));
  
  registerProcessorParameter("Debug",
			     "Debugging level 0-3",
			     _debug,
			     int(0));
  
}



void ECALPreClustering::init()
{
  // always

  printParameters();
  _nEvents=0;
  _nPreClusters=0;
}



void ECALPreClustering::processRunHeader(LCRunHeader * run)
{
}


void ECALPreClustering::processEvent(LCEvent * evt)
{
  _nEvents++;
  _nPreClustersInEvent=0;
  if(_debug>1) {
    cout << endl;
    cout << "Preclustering Event " << evt->getEventNumber() << endl;
  }
  
  vector<ExtendedHit* > hitVec;
  multimap<int,ClusterImpl* > preClusVec;
  
  if(_debug>2) 
    cout << "Preparing Hits..." << endl;
  // prepare hits for easier clustering
  PrepareHits(evt,hitVec);

  if(_debug>2) 
    cout << "Preclustering Hits..." << endl;
  // begin PreClustering:
  PreCluster(hitVec, preClusVec);

  if(_debug>2) 
    cout << "Writing Collections..." << endl;
  // write PreClusters to as LC collection
  WritePreClusters(evt, preClusVec);
}



void ECALPreClustering::check(LCEvent * evt)
{

}



void ECALPreClustering::end()
{
  cout << endl;
  cout << "ECALPreClustering Report: " << endl;
  cout << _nEvents << " events processed" << endl;
  cout << "Found " << _nPreClusters << " PreCluster" << endl;
  cout << "= " << _nPreClusters/_nEvents  << " /event" << endl;
}



float ECALPreClustering::GetDistance(float *a_pos,float *s_pos)
{
  float dist=sqrt((s_pos[0]-a_pos[0])*(s_pos[0]-a_pos[0])+(s_pos[1]-a_pos[1])*(s_pos[1]-a_pos[1])+(s_pos[2]-a_pos[2])*(s_pos[2]-a_pos[2]));
  return dist;
}



void ECALPreClustering::PrepareHits(const LCEvent *evt, vector<ExtendedHit* > &hitVec) 
{
  LCCollection *hitColl = 0;
  try {
    hitColl = evt->getCollection(_ecalBarrelHitCollectionName);
  }
  catch(DataNotAvailableException &exc) {
    if (_debug>0)
      std::cerr << "In ECALPreClustering::processEvent "
		<< exc.what() << std::endl;
    //   return;
  }
  int NHits = 0;
  if(hitColl) {
    NHits = hitColl->getNumberOfElements();
    if(_debug>2) 
      cout << "Barrel Hit Collection has " << NHits << " hits" << endl;
    for (int hit_i = 0; hit_i < NHits; hit_i++) {
      CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit*>( hitColl->getElementAt(hit_i) );
      ExtendedHit *a_ext_hit = new ExtendedHit();
      a_ext_hit->hit = a_hit;
      a_ext_hit->isClustered = 0;
      hitVec.push_back(a_ext_hit);
      //}
    }
  }
  else
    if (_debug>0) cout << "No Barrel Hits" << endl;

  hitColl = 0;
  try {
    hitColl = evt->getCollection(_ecalEndcapHitCollectionName);
  }
  catch(DataNotAvailableException &exc) {
    if (_debug>0)
      std::cerr << "In ECALPreClustering::processEvent "
		<< exc.what() << std::endl;
    //return;
  }

  if(hitColl) {
    CellIDDecoder<CalorimeterHit> decoder(hitColl); 
    NHits = hitColl->getNumberOfElements();
    if(_debug>2) 
      cout << "Endcap Hit Collection has " << NHits << " hits" << endl;
    for (int hit_i = 0; hit_i < NHits; hit_i++) {
      CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit*>( hitColl->getElementAt(hit_i) );
      /*
      int cell_i = decoder(a_hit)["I"];
      int cell_j = decoder(a_hit)["J"];
      int module = decoder(a_hit)["M"];
      int stave = decoder(a_hit)["S-1"];
      int layer = decoder(a_hit)["K-1"];
      cout << "Endcap hit:" << cell_i << " " << cell_j << " " << stave << " " << module << endl;
      */
      ExtendedHit *a_ext_hit = new ExtendedHit();
      a_ext_hit->hit = a_hit;
      a_ext_hit->isClustered = 0;
      hitVec.push_back(a_ext_hit);
      //}
    }
  }
  else
    if (_debug>0)
      cout << "No Endcap Hits" << endl;


  hitColl = 0;
  try {
    hitColl = evt->getCollection(_ecalOtherHitCollectionName);
  }
  catch(DataNotAvailableException &exc) {
    if (_debug>0)
      std::cerr << "In ECALPreClustering::processEvent "
		<< exc.what() << std::endl;
    //return;
  }

  if(hitColl) {
    NHits = hitColl->getNumberOfElements();
    CellIDDecoder<CalorimeterHit> decoder(hitColl); 
    if(_debug>2) 
      cout << "Ring Hit Collection has " << NHits << " hits" << endl;
    for (int hit_i = 0; hit_i < NHits; hit_i++) {
      CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit*>( hitColl->getElementAt(hit_i) );
      /*  
    int cell_i = decoder(a_hit)["I"];
      int cell_j = decoder(a_hit)["J"];
      int module = decoder(a_hit)["M"];
      int stave = decoder(a_hit)["S-1"];
      int layer = decoder(a_hit)["K-1"];
      cout << "Ring hit:" << cell_i << " " << cell_j << " " << stave << " " << module << endl;
*/
      ExtendedHit *a_ext_hit = new ExtendedHit();
      a_ext_hit->hit = a_hit;
      a_ext_hit->isClustered = 0;
      hitVec.push_back(a_ext_hit);
      //}
    }
  }
  else
    if (_debug>0) cout << "No Ring Hits" << endl;
  
  if(_debug>2) 
    cout << "New extended Hit vector has " << hitVec.size() << " hits" << endl;
  // sort hits by energy deposit
  if(hitVec.size()>0)
    sort(hitVec.begin(),hitVec.end(),ExtendedHit::higherEnergy);
  else return;
}



void ECALPreClustering::PreCluster( vector<ExtendedHit* > &hitVec, multimap<int,ClusterImpl* > &preClusVec)
{
  // 1. regard every non-clustered hit as potentiel source of a cluster
  
  vector<ExtendedHit* >::iterator source_it;
  for(source_it=hitVec.begin();source_it!=hitVec.end();source_it++) {
    
    if(!((*source_it)->isClustered)) {   // if non-clustered start clustering
      
      ClusterImpl *a_cluster = new ClusterImpl();  // open new cluster 
      if(_debug>2) 
	cout << "Opened a new cluster with hit " << (*source_it)->hit->getCellID0() << " as source" << endl;
      a_cluster->addHit((*source_it)->hit,1);
      (*source_it)->isClustered = true;
      
      ClusterNearHits(*source_it, hitVec, *a_cluster);
      if(_debug>2) 
	cout << "This Cluster has now " << (a_cluster->getCalorimeterHits()).size() << " hits" << endl;
      int nHits = (a_cluster->getCalorimeterHits()).size();
      preClusVec.insert(make_pair(nHits,a_cluster));
    }
    else   // else do nothing
      ;
  }

}



void ECALPreClustering::ClusterNearHits(ExtendedHit *s_hit, vector<ExtendedHit* > &hits, ClusterImpl &cluster)
{
  const float *spos=s_hit->hit->getPosition();
  float s_pos[3]={spos[0],spos[1],spos[2]};
  vector<ExtendedHit* >::iterator hit_it;
  for(hit_it=hits.begin(); hit_it!=hits.end(); hit_it++) {
    if(!((*hit_it)->isClustered)) {
      const float *apos=(*hit_it)->hit->getPosition();
      float a_pos[3]={apos[0],apos[1],apos[2]};
      float dist=GetDistance(s_pos,a_pos);
      if(dist < _distanceCut) {
	if(_debug>2) 
	  cout << "Now clustering hit " << (*hit_it)->hit->getCellID0() << " with distance " << dist << " from hit " << s_hit->hit->getCellID0() << endl;
	(*hit_it)->isClustered=true;
	cluster.addHit((*hit_it)->hit,1);
	ClusterNearHits((*hit_it), hits, cluster);
      }
    }
  }
}



void ECALPreClustering::WritePreClusters(LCEvent *evt, multimap<int,ClusterImpl* > &preClusVec)
{
  LCCollection *hitColl = 0;
  try {
    hitColl = evt->getCollection(_ecalBarrelHitCollectionName);
  }
  catch(DataNotAvailableException &exc) {
    //std::cerr << "In ECALPreClustering::WritePreClusters "
    //      << exc.what() << std::endl;
    //return;
  }
  if(hitColl==0) {
    //cout << "Taking Endcap collection" << endl;
    try {
      hitColl = evt->getCollection(_ecalEndcapHitCollectionName);
    }
    catch(DataNotAvailableException &exc) {
      //std::cerr << "In ECALPreClustering::WritePreClusters "
      //	<< exc.what() << std::endl;
      return;
    }
  }
  LCCollectionVec *preClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  preClusterColl->setFlag( preClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  preClusterColl->parameters().setValue("CellIDEncoding", hitColl->getParameters().getStringVal("CellIDEncoding"));

  _nPreClusters += preClusVec.size();
  _nPreClustersInEvent = preClusVec.size();
  multimap<int,ClusterImpl*>::reverse_iterator preClusIt;
  for (preClusIt=preClusVec.rbegin();preClusIt!=preClusVec.rend() ; preClusIt++) {
    if(preClusIt->first>=_minHits) {
      if(_debug>2) 
	cout << "Adding PreCluster with " << (preClusIt->first) << " hits" << endl;
      preClusterColl->addElement(preClusIt->second);
      _nPreClusters ++;
      _nPreClustersInEvent++;
    }
  }
  if(_debug>2) 
    cout << preClusterColl->getNumberOfElements() << " PreClusters in event " << evt->getEventNumber() << endl;

  evt->addCollection(preClusterColl, _ecalPreClusterCollectionName);
  if(_debug>2) 
    cout << "PreCluster Collection written!" << endl;
}






