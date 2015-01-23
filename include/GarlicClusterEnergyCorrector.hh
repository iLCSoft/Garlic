#ifndef Garlic_Cluster_Energy_Corrector
#define Garlic_Cluster_Energy_Corrector

#include "GarlicExtendedCluster.hh"

class GarlicClusterEnergyCorrector {

public:
  
  GarlicClusterEnergyCorrector(){}
  ~GarlicClusterEnergyCorrector(){}
  
  float getCorrectedEnergy(GarlicExtendedCluster* cl);

private:



};


#endif
