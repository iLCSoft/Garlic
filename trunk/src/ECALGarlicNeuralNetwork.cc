#include "ECALGarlicNeuralNetwork.hh"

#include <iostream>
#include <vector>
#include <assert.h>

#include "ECALGarlicAlgorithmParameters.hh"


//#include "weights/MLP_Garlic_en0_02_notrkAllbck.C"
#include "weights/MLP_Garlic_en0_02_notrkAllbck.h"
#include "weights/MLP_Garlic_en02_05_notrkAllbck.h"
#include "weights/MLP_Garlic_en05_1_notrkAllbck.h"
#include "weights/MLP_Garlic_en1_3_notrkAllbck.h"
#include "weights/MLP_Garlic_en3_10_notrkAllbck.h"
#include "weights/MLP_Garlic_en10_100_notrkAllbck.h"

#include "weights/MLP_Garlic_en0_02_trkPibck.h"
#include "weights/MLP_Garlic_en02_05_trkPibck.h"
#include "weights/MLP_Garlic_en05_1_trkPibck.h"
#include "weights/MLP_Garlic_en1_3_trkPibck.h"
#include "weights/MLP_Garlic_en3_10_trkPibck.h"
#include "weights/MLP_Garlic_en10_100_trkPibck.h"


using std::vector;
using std::cout;
using std::endl;

void ECALGarlicNeuralNetwork::setup() {

  std::cout << "hello from ECALGarlicNeuralNetwork::setup" << std::endl;

  cutVals = *(ECALGarlicAlgorithmParameters::Instance().GetMLPCuts());
  assert (cutVals.size() == nNN );

  vector<std::string> varnames;

  varnames.push_back("dirErr");
  varnames.push_back("depth");
  varnames.push_back("enFracRel_x0_10");
  varnames.push_back("hitsMeanEn");
  varnames.push_back("hitsRMSEn/hitsMeanEn");
  varnames.push_back("fracDim4");
  varnames.push_back("transRMSmin");
  varnames.push_back("distToTrackPOS");
  varnames.push_back("angleToTrackPOS");

  MLPResponses_trk[0] = new ReadMLP_Garlic_en0_02_trkPibck(varnames);
  MLPResponses_trk[1] = new ReadMLP_Garlic_en02_05_trkPibck(varnames);
  MLPResponses_trk[2] = new ReadMLP_Garlic_en05_1_trkPibck(varnames);
  MLPResponses_trk[3] = new ReadMLP_Garlic_en1_3_trkPibck(varnames);
  MLPResponses_trk[4] = new ReadMLP_Garlic_en3_10_trkPibck(varnames);
  MLPResponses_trk[5] = new ReadMLP_Garlic_en10_100_trkPibck(varnames);

  varnames.clear();

  varnames.push_back("dirErr");
  varnames.push_back("depth");
  varnames.push_back("enFracRel_x0_10");
  varnames.push_back("hitsMeanEn");
  varnames.push_back("hitsRMSEn/hitsMeanEn");
  varnames.push_back("fracDim4");
  varnames.push_back("transRMSmin");

  MLPResponses_noTrk[0] = new ReadMLP_Garlic_en0_02_notrkAllbck(varnames);
  MLPResponses_noTrk[1] = new ReadMLP_Garlic_en02_05_notrkAllbck(varnames);
  MLPResponses_noTrk[2] = new ReadMLP_Garlic_en05_1_notrkAllbck(varnames);
  MLPResponses_noTrk[3] = new ReadMLP_Garlic_en1_3_notrkAllbck(varnames);
  MLPResponses_noTrk[4] = new ReadMLP_Garlic_en3_10_notrkAllbck(varnames);
  MLPResponses_noTrk[5] = new ReadMLP_Garlic_en10_100_notrkAllbck(varnames);

  return;
}

ECALGarlicNeuralNetwork::~ECALGarlicNeuralNetwork() {
  for (int i=0; i<nNN; i++) {
    if (MLPResponses_trk[i]) {
      delete MLPResponses_trk[i];
      MLPResponses_trk[i]=NULL;
    }
    if (MLPResponses_noTrk[i]) {
      delete MLPResponses_noTrk[i];
      MLPResponses_noTrk[i]=NULL;
    }
  }
}



std::pair <float, bool> ECALGarlicNeuralNetwork::getNNoutput(ExtendedCluster2* clus) {

  vector<double> inputVec;

  int ien(-1);

  bool hasTrack = clus->getTrackDist_proj()<100;

  inputVec.push_back(clus->getClusterPointAngle());
  inputVec.push_back(clus->getMeanDepth());
  inputVec.push_back(clus->getRelLongEn()[1]);
  inputVec.push_back(clus->getHitMeanEn());
  inputVec.push_back(clus->getHitRMSEn()/clus->getHitMeanEn());
  inputVec.push_back(clus->getFractalDimension()[1]);
  inputVec.push_back(clus->getTransverseRMS().first);
  if (hasTrack) {
    inputVec.push_back(clus->getTrackDist_proj());
    inputVec.push_back(clus->getTrackAng_proj());
  }

  if (clus->getEnergy()<0.2) {
    ien=0;
  } else if (clus->getEnergy()<0.5) {
    ien=1;
  } else if (clus->getEnergy()<1.0) {
    ien=2;
  } else if (clus->getEnergy()<3.0) {
    ien=3;
  } else if (clus->getEnergy()<10.0) {
    ien=4;
  } else {
    ien=5;
  }


  float MLPout = hasTrack ? MLPResponses_trk[ien]->GetMvaValue(inputVec) : MLPResponses_noTrk[ien]->GetMvaValue(inputVec);
  bool MLPcut = MLPout>cutVals[ien];

  return std::pair <float, bool> (MLPout, MLPcut);

}

