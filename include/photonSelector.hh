#ifndef PHOTONSELECTOR_hh
#define PHOTONSELECTOR_hh

#include "ECALGarlicExtendedCluster.hh"

class photonSelector {

public:

  photonSelector(bool verb=false) {verbose=verb;}
  ~photonSelector() {}

  //---------------
  
  int photon_select(ExtendedCluster2* ecl);

  int photon_longProfile(ExtendedCluster2* ecl);
  int photon_transProfile(ExtendedCluster2* ecl);
  int photon_hitEnergies(ExtendedCluster2* ecl);
  int photon_pointing(ExtendedCluster2* ecl);

  //----------------

  int photon_select_pointAng(ExtendedCluster2* ecl);
  int photon_select_start(ExtendedCluster2* ecl);
  int photon_select_reldepth(ExtendedCluster2* ecl);
  int photon_select_relrelLongE0(ExtendedCluster2* ecl);
  int photon_select_relrelLongE1(ExtendedCluster2* ecl);
  int photon_select_relrelLongE2(ExtendedCluster2* ecl);
  int photon_select_hitEnDistr(ExtendedCluster2* ecl);
  int photon_select_fracDim(ExtendedCluster2* ecl);
  int photon_select_transRmsA(ExtendedCluster2* ecl);
  int photon_select_transRmsB(ExtendedCluster2* ecl);
  int photon_select_eccen(ExtendedCluster2* ecl);
  int photon_select_molA(ExtendedCluster2* ecl);
  int photon_select_molB(ExtendedCluster2* ecl);
  int photon_select_eccenMol(ExtendedCluster2* ecl);
  int photon_select_earlyMolB(ExtendedCluster2* ecl);
  int photon_select_fracPLay(ExtendedCluster2* ecl);

private:

  bool verbose;

};


#endif
