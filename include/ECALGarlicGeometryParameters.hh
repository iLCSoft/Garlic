#ifndef ECALGARLICGEOMETRYPARAMS_HH_
#define ECALGARLICGEOMETRYPARAMS_HH_

#define MAX_NUMBER_OF_LAYERS 200

#include <cmath>
#include <UTIL/CellIDDecoder.h>

#include "ECALGarlicExtendedObjects.hh"
using namespace ECALGarlicExtendedObjects;

class ECALGarlicGeometryParameters {

public:

  static ECALGarlicGeometryParameters& Instance() {
    static ECALGarlicGeometryParameters geompars;
    return geompars;
  }

  void Set_rOfBarrel (float f) {_rOfBarrel=f;}
  void Set_zOfBarrel (float f) {_zOfBarrel=f;}
  void Set_rOfEndcap (float f) {_rOfEndcap=f;}
  void Set_zOfEndcap (float f) {_zOfEndcap=f;}
  void Set_activeThickness (float f) {_activeThickness=f;}
  void Set_passiveThickness (float f) {_passiveThickness=f;}
  //  void Set_firstEndcapLayerOffset (float f) {_firstEndcapLayerOffset=f;}
  //  void Set_firstBarrelLayerOffset (float f) {_firstBarrelLayerOffset=f;}
  void Set_bField (float f) {_bField=f;}
  void Set_rInnerEcalEndcap (float f) {_rInnerEcalEndcap=f;}
  void Set_cosOfBarrel (float f) {_cosOfBarrel=f;}
  void Set_symmetry (int   f) {_symmetry=f;}
  void Set_nPseudoLayers (int   f) {_nPseudoLayers=f;}
  void Set_nCellsPerWafer (int   f) {_nCellsPerWafer=f;}

  void Set_nBarrelEcalLayers (int it) {_nBarrelEcalLayers=it;}
  void Set_nEndcapEcalLayers (int it) {_nEndcapEcalLayers=it;}

  int _nBarrelLayers;
  void Set_positionBarrelLayer (float* f) {
    for (uint i=0; i<MAX_NUMBER_OF_LAYERS; i++)
      _positionBarrelLayer[i]=f[i];
  }
  void Set_absThicknessBarrelLayer (float* f) {
    for (uint i=0; i<MAX_NUMBER_OF_LAYERS; i++)
      _absThicknessBarrelLayer[i]=f[i];
  }
  void Set_absThicknessEndcapLayer (float* f) {
    for (uint i=0; i<MAX_NUMBER_OF_LAYERS; i++)
      _absThicknessEndcapLayer[i]=f[i];
  }
  void Set_padSizeEcal (float* f) {
    for (uint i=0; i<MAX_NUMBER_OF_LAYERS; i++)
      _padSizeEcal[i]=f[i];
  }
  void Set_positionEndcapLayer (float* f) {
    for (uint i=0; i<MAX_NUMBER_OF_LAYERS; i++)
      _positionEndcapLayer[i]=f[i];
  }

  void Set_barrelStaveDir (std::vector< std::vector < float > > gg) {_barrelStaveDir=gg;}
  void Set_defaultDecoder (CellIDDecoder<CalorimeterHit>* dec) {_defaultDecoder=dec;}

  void Set_absorberRadiationLength (float x0) {_absorberRadiationLength=x0;}

  float Get_rOfBarrel () {return _rOfBarrel;}
  float Get_zOfBarrel () {return _zOfBarrel;}
  float Get_rOfEndcap () {return _rOfEndcap;}
  float Get_zOfEndcap () {return _zOfEndcap;}
  float Get_activeThickness () {return _activeThickness;}
  float Get_passiveThickness () {return _passiveThickness;}
  //  float Get_firstEndcapLayerOffset () {return _firstEndcapLayerOffset;}
  //  float Get_firstBarrelLayerOffset () {return _firstBarrelLayerOffset;}
  float Get_bField () {return _bField;}
  float Get_rInnerEcalEndcap () {return _rInnerEcalEndcap;}
  float Get_cosOfBarrel () {return _cosOfBarrel;}
  int   Get_symmetry () {return _symmetry;}
  int   Get_nPseudoLayers () {return _nPseudoLayers;}
  int   Get_nCellsPerWafer () {return _nCellsPerWafer;}
  float Get_absorberRadiationLength () {return _absorberRadiationLength;}

  int Get_nBarrelEcalLayers () {return _nBarrelEcalLayers;}
  int Get_nEndcapEcalLayers () {return _nEndcapEcalLayers;}


  float * Get_positionBarrelLayer () {return _positionBarrelLayer;}
  float * Get_absThicknessBarrelLayer () {return _absThicknessBarrelLayer;}
  float * Get_absThicknessEndcapLayer () {return _absThicknessEndcapLayer;}
  float * Get_padSizeEcal () {return _padSizeEcal;}
  float * Get_positionEndcapLayer () {return _positionEndcapLayer;}

  std::vector< std::vector < float > > Get_barrelStaveDir () {return _barrelStaveDir;}
  CellIDDecoder<CalorimeterHit>* Get_defaultDecoder () {return _defaultDecoder;}

private:


  ECALGarlicGeometryParameters() {init();}

  ECALGarlicGeometryParameters (const ECALGarlicGeometryParameters&);
  ECALGarlicGeometryParameters& operator= (const ECALGarlicGeometryParameters&);

  double _rOfBarrel;
  double _zOfBarrel;
  double _rOfEndcap;
  double _zOfEndcap;

  float _activeThickness;
  float _passiveThickness;
  //  float _firstEndcapLayerOffset;
  //  float _firstBarrelLayerOffset;
  float _bField;
  float _positionBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _absThicknessBarrelLayer[MAX_NUMBER_OF_LAYERS];
  float _absThicknessEndcapLayer[MAX_NUMBER_OF_LAYERS];
  float _padSizeEcal[MAX_NUMBER_OF_LAYERS];
  float _positionEndcapLayer[MAX_NUMBER_OF_LAYERS];
  float _rInnerEcalEndcap;
  float _cosOfBarrel;

  int _symmetry;
  int _nPseudoLayers;
  int _nCellsPerWafer;

  int _nBarrelEcalLayers;
  int _nEndcapEcalLayers;

  float _absorberRadiationLength;

  static const float pi;

  std::vector< std::vector < float > > _barrelStaveDir;

  CellIDDecoder<CalorimeterHit>* _defaultDecoder;

  void init() {
    _rOfBarrel=-999;
    _zOfBarrel=-999;
    _rOfEndcap=-999;
    _zOfEndcap=-999;
    _nBarrelEcalLayers=-999;
    _nEndcapEcalLayers=-999;
    _activeThickness=-999;
    _passiveThickness=-999;
    //    _firstEndcapLayerOffset=-999;
    //    _firstBarrelLayerOffset=-999;
    _bField=-999;
    _rInnerEcalEndcap=-999;
    _cosOfBarrel=-999;
    _symmetry=-999;
    _nPseudoLayers=-999;
    _nCellsPerWafer=-999;
    _absorberRadiationLength=-999;
  }

  
};


#endif
