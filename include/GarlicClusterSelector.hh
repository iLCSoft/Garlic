#ifndef CLUSTERSELECTOR_hh
#define CLUSTERSELECTOR_hh

#include "GarlicExtendedCluster.hh"
#include <map>
#include <vector>
#include <utility>

class GarlicClusterSelector {

public:

  GarlicClusterSelector(std::string cutFileName, int pdg=0, bool verb=false) {
    verbose=verb;
    _pdg=pdg;
    readInCutValues(cutFileName);
  }
  ~GarlicClusterSelector() {}

  //---------------

  void readInCutValues(std::string cutFileName);
  
  std::map < std::string, float > cluster_select(GarlicExtendedCluster* ecl);

  void setPdg(int pdg) {_pdg=pdg;};
  int getPdg() {return _pdg;}

  enum { REJECT, VERYLOOSE, LOOSE, TIGHT };

private:

  enum { N_PERC=4 };

  bool verbose;

  std::map < std::pair < int, std::string > , std::map < float, std::vector < float > > > _varCuts;
  std::map < std::pair < int, std::string > , std::vector < std::string > > _varClasses;

  std::vector < int > cluster_classSelection(GarlicExtendedCluster* ecl, std::string className);

  int cluster_select_variable(GarlicExtendedCluster* ecl, std::string variable);

  int _pdg;

  void dump_cutvalues();

};


#endif
