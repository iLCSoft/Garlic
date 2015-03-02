#ifndef ECALGarlicNeuralNetwork_HH_
#define ECALGarlicNeuralNetwork_HH_

class IClassifierReader;

#include "ECALGarlicGeometryHelpers.hh"
#include "ECALGarlicExtendedCluster.hh"

#include <vector>

class ECALGarlicNeuralNetwork : public ECALGarlicGeometryHelpers {

public:

  static ECALGarlicNeuralNetwork& Instance() {
    static ECALGarlicNeuralNetwork nn;
    return nn;
  }

  std::pair <float, bool> getNNoutput(ExtendedCluster2* clus);

private:

  ECALGarlicNeuralNetwork() {
    for (int i=0; i<nNN; i++) {
      MLPResponses_trk[i]=NULL;
      MLPResponses_noTrk[i]=NULL;
    }
    setup();
  }

  ECALGarlicNeuralNetwork (const ECALGarlicNeuralNetwork&);
  ECALGarlicNeuralNetwork& operator= (const ECALGarlicNeuralNetwork&);

  ~ECALGarlicNeuralNetwork();

  void setup();

  enum {nNN=6};

  IClassifierReader* MLPResponses_trk[nNN];
  IClassifierReader* MLPResponses_noTrk[nNN];

  std::vector <float> cutVals;


};

#endif
