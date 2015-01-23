#include "GarlicClusterEnergyCorrector.hh"

float GarlicClusterEnergyCorrector::getCorrectedEnergy(GarlicExtendedCluster* cl) {

  // this is just a placeholder for now
  // eventually include e.g. theta, phi, energy-dependent corrections.

  return cl->getEnergy();

}
