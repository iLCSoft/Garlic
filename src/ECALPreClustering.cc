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

#include <ECALGarlicExtendedHit.hh>

using namespace marlin;
using namespace lcio;
using namespace std;

ECALPreClustering anECALPreClustering;

ECALPreClustering::ECALPreClustering() : Processor("ECALPreClustering_v30") {

  // processor description
  _description = "Preclusters ECAL hits to define ROIs for Garlic";

  registerProcessorParameter("Debug",
			     "Debugging level 0-3",
			     _debug,
			     int(0));
  
  StringVec ecalHitCollNames;
  ecalHitCollNames.push_back("ECALEndcap");
  ecalHitCollNames.push_back("ECALBarrel");
  ecalHitCollNames.push_back("ECALOther");
  registerProcessorParameter( "EcalHitCollections",
			      "names of ECAL hit collections",
			      _ecalHitCollNames,
			      ecalHitCollNames);

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
  
}

ECALPreClustering::~ECALPreClustering() {
}


void ECALPreClustering::init()
{
  printParameters();
  _nEvents=0;
  _nPreClusters=0;
  _decodeString="";
  return;
}



void ECALPreClustering::processRunHeader(LCRunHeader * run)
{
  return;
}


void ECALPreClustering::processEvent(LCEvent * evt)
{
  _nEvents++;
  _nPreClustersInEvent=0;
  if(_debug>1 || evt->getEventNumber()%200==0) {
    cout << endl << "Preclustering Event " << evt->getEventNumber() << endl;
  }

  vector<ExtendedHit* > hitVec;
  multimap<int,ClusterImpl* > preClusVec;
  
  if(_debug>2) cout << "Preparing Hits..." << endl;
  // prepare hits for easier clustering
  PrepareHits(evt,hitVec);

  if(_debug>2) cout << "Preclustering Hits..." << endl;
  // begin PreClustering:
  PreCluster(hitVec, preClusVec);

  if(_debug>2) cout << "Writing Collections..." << endl;
  // write PreClusters to as LC collection
  WritePreClusters(evt, preClusVec);

  // clean up ExtendedHits;
  hitVec.erase(hitVec.begin(), hitVec.end());

  return;
}



void ECALPreClustering::check(LCEvent * evt)
{
  return;
}



void ECALPreClustering::end()
{
  cout << endl;
  cout << "ECALPreClustering Report: " << endl;
  cout << _nEvents << " events processed" << endl;
  cout << "Found " << _nPreClusters << " PreCluster" << endl;
  if (_nEvents>0)
    cout << "= " << _nPreClusters/_nEvents  << " /event" << endl;
  return;
}



float ECALPreClustering::GetDistance(float *a_pos,float *s_pos)
{
  //  float dist=sqrt((s_pos[0]-a_pos[0])*(s_pos[0]-a_pos[0])+(s_pos[1]-a_pos[1])*(s_pos[1]-a_pos[1])+(s_pos[2]-a_pos[2])*(s_pos[2]-a_pos[2]));

  float dist(0);
  for (int i=0; i<3; i++) dist+=pow(s_pos[i]-a_pos[i], 2);
  dist=sqrt(dist);
  return dist;
}



void ECALPreClustering::PrepareHits(const LCEvent *evt, vector<ExtendedHit* > &hitVec) {

  /*

  for (size_t ic=0; ic<_ecalHitCollNames.size(); ic++) {
    string colname = _ecalHitCollNames[ic];
    LCCollection *hitColl = 0;
    try {
      hitColl = evt->getCollection(colname);
      if (hitColl) {

        if (_decodeString=="")
	  _decodeString=hitColl->getParameters().getStringVal("CellIDEncoding");

	int NHits = hitColl->getNumberOfElements();
	if(_debug>1) cout << colname << " hit Collection has " << NHits << " hits" << endl;
	for (int hit_i = 0; hit_i < NHits; hit_i++) {
	  CalorimeterHit *a_hit = dynamic_cast<CalorimeterHit*>( hitColl->getElementAt(hit_i) );

	  // check for infinite position
	  for (int ikl=0; ikl<3; ikl++)
	    assert (fabs(a_hit->getPosition()[ikl])<99999999);

	  ExtendedHit *a_ext_hit = new ExtendedHit();
	  a_ext_hit->hit = a_hit;
	  a_ext_hit->isClustered = 0;
	  hitVec.push_back(a_ext_hit);
	}
      }
      else
	if (_debug>0) cout << "No hits found in collection " << colname << endl;

    }
    catch(DataNotAvailableException &exc) {
      if (_debug>0)
	std::cerr << "In ECALPreClustering::processEvent "
		  << exc.what() << std::endl;
    }

  }

  
  if(_debug>2) cout << "New extended Hit vector has " << hitVec.size() << " hits" << endl;

  // daniel decides this is completely superfluous
  //   // sort hits by energy deposit
  if(hitVec.size()>0)
    sort(hitVec.begin(),hitVec.end(),ExtendedHit::higherEnergy);

  */

  return;
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
  }
  return;
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

  LCCollectionVec *preClusterColl = new LCCollectionVec(LCIO::CLUSTER);
  preClusterColl->setFlag( preClusterColl->getFlag() | 1 << LCIO::CLBIT_HITS );
  preClusterColl->parameters().setValue("CellIDEncoding", _decodeString);

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

  if (preClusterColl->getNumberOfElements() > 0) {
    evt->addCollection(preClusterColl, _ecalPreClusterCollectionName);
    if(_debug>2) cout << "PreCluster Collection written!" << endl;
  }
  return;
}






