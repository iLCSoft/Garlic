#ifndef ECAL_Garlic_Cluster_Energy_Corrector
#define ECAL_Garlic_Cluster_Energy_Corrector

#include "ECALGarlicExtendedCluster.hh"

class ECALGarlicClusterEnergyCorrector {

public:
  
  ECALGarlicClusterEnergyCorrector(){}
  ~ECALGarlicClusterEnergyCorrector(){}
  
  float getCorrectedEnergy(ExtendedCluster2* cl);

private:



};


#endif
