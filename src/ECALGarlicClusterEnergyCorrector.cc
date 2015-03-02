#include "ECALGarlicClusterEnergyCorrector.hh"

float ECALGarlicClusterEnergyCorrector::getCorrectedEnergy(ExtendedCluster2* cl) {

  // this is just a placeholder for now
  // eventually include e.g. theta, phi, energy-dependent corrections.

  return cl->getEnergy();

}
