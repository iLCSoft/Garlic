#include "GarlicExtendedHit.hh"
#include "GarlicAlgorithmParameters.hh"

#include <iostream>
using std::cout;
using std::endl;
#include <assert.h>

void PointInCalo::calculateX0depth() {

  if (_pseudoLayer<0) {
    assignPseudoLayer();
  }

  if (getZone() == CALHITZONE_BARREL) {
    _x0depth = getBarrelCorrAbsThickness();

    if (_x0depth>150)
      cout << "strangely deep barrel hit! " << _x0depth << " zone " << getZone() << " " << _pos[0] << " " << _pos[1] << " " << _pos[2] <<
        " : " << _pseudoLayer << endl;

  } else if (getZone() == CALHITZONE_ENDCAP || getZone() == CALHITZONE_RING) {
    // normal to front face
    float norm[3] = {0,0,1};
    // angle between norm and line of flight
    float angle = GarlicGeometryHelpers::angle(norm, _pos);
    _x0depth = getAbsThickness(_pseudoLayer)/fabs(cos(angle));

    if (_x0depth>100) cout << "strangely deep end/ring hit! " <<
      _x0depth << " zone " << getZone() << " " << _pos[0] << " " << _pos[1] << " " << _pos[2] << " " << angle << " " << cos(angle) << " " << _pseudoLayer << endl;


    if (getZone() == CALHITZONE_ENDCAP) {
      // also check for crossing barrel in overlap region
      // project onto end face of barrel
      float hitProj[3];
      float scale = GarlicGeometryParameters::Instance().Get_zOfBarrel()/fabs(_pos[2]);
      GarlicGeometryHelpers::scale(_pos, scale, hitProj);
      // which barrel pseudolayer does this correspond to?
      PointInCalo projection(hitProj);
      int pzone =  projection.getZone();
      if (pzone == CALHITZONE_BARREL) {
        int extrabar = projection.getBarrelPseudoLayer();
        if (extrabar>=0) {
          float extraX0 = projection.getBarrelCorrAbsThickness();
          if (extraX0>0) _x0depth+=extraX0;
          else if (extrabar>1) cout << "got strange extra x0 from barrel overlap..." << extraX0 << endl;
        }
      }
    }
  }

  // convert from mm of absorber to x0

  _x0depth/=GarlicGeometryParameters::Instance().Get_absorberRadiationLength();

  return;
}

float PointInCalo::getAbsThickness(int layer) {
  if (layer<0 || layer>=GarlicGeometryParameters::Instance().Get_nBarrelEcalLayers()) return -999;
  float absThick(0);
  for (int il=0; il<layer; il++) {
    absThick+=GarlicGeometryParameters::Instance().Get_absThicknessBarrelLayer()[il];
  }
  return absThick;
}

float PointInCalo::getBarrelCorrAbsThickness(int layer, const float* pos) {

  if (getZone()!=CALHITZONE_BARREL) {
    cout << "error from PointInCalo::getBarrelCorrAbsThickness, trying to get barrel info about point not in barrel! " << getZone() << endl;
    cout << "position: " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    return -999;
  }
  int pstave = getPseudoStave();
  float pstavePhi = getPseudoStavePhi(pstave);
  float normal[3] = {cos(pstavePhi), sin(pstavePhi), 0};
  float cosang = cosangle(normal, pos);

  float thickness = getAbsThickness(layer)/cosang;

  if (thickness>200) cout << "PointInCalo::getBarrelCorrAbsThickness very large barrel thickness " << thickness << " stave " << pstave << " stavePhi " << pstavePhi <<
    " normal:" << normal[0] << " " << normal[1] << " " << normal[2] <<
    " position: " << pos[0] << " " << pos[1] << " " << pos[2] <<
    " cos angle: " << cosang << " layer " << layer << " absorberthickness " << getAbsThickness(layer) << endl;

  //if (thickness<=0)
  if (layer>1 && thickness<=0)
    cout << "PointInCalo::getBarrelCorrAbsThickness very small barrel thickness! " << thickness << " stave " << pstave << " stavePhi " << pstavePhi <<
      " normal:" << normal[0] << " " << normal[1] << " " << normal[2] <<
      " position: " << pos[0] << " " << pos[1] << " " << pos[2] <<
      " cos angle: " << cosang << " layer " << layer << " absorberthickness " << getAbsThickness(layer) << endl;


  return thickness;
}


void PointInCalo::assignZone(bool check) {
  int zone = CALHITZONE_UNKNOWN;
  double prodRadiusMax=0.;

  float bar_end_boundary = ( GarlicGeometryParameters::Instance().Get_zOfBarrel() + GarlicGeometryParameters::Instance().Get_zOfEndcap() )/2.;

  // check if position is inside ECAL, and in which zone
  if (fabs(_pos[2])<=bar_end_boundary) {

    for(int istave = 0; istave < GarlicGeometryParameters::Instance().Get_symmetry(); ++istave){
      double prodRadius =
        _pos[0]*GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][0]+
        _pos[1]*GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][1];
      prodRadius/=sqrt( pow(GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][0],2) +
                        pow(GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][1],2) );

      if(prodRadius>prodRadiusMax) prodRadiusMax=prodRadius;
    }

    //    float rmin = GarlicGeometryParameters::Instance().Get_rOfBarrel();
    float rmin = GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [0];

    //float rmax = GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [GarlicGeometryParameters::Instance().Get_nBarrelEcalLayers()-1];
    //float thick = rmax - GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [GarlicGeometryParameters::Instance().Get_nBarrelEcalLayers()-2];
    //rmax += thick; // be generous to deal with gear file definitions...

    float rmax = GarlicGeometryParameters::Instance().Get_rMaxOfBarrel();

    // thick =
    //   GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [1] -
    //   GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [0];
    //
    // rmin -= thick; // be generous to deal with gear file definitions...

    if ( prodRadiusMax>=rmin && prodRadiusMax<=rmax ) {
      zone = CALHITZONE_BARREL;
    } else {

      if (check) {

        if (prodRadiusMax>rmax) cout << "BARREL point beyond ECAL... point radius = " << prodRadiusMax << " " << rmax << endl;
        if (prodRadiusMax<rmin) cout << "BARREL point before ECAL... point radius = " << prodRadiusMax << " " << rmin << endl;
        cout << " Z pos, boundary = " << _pos[2] << " " << bar_end_boundary << endl;

        // cout << _pos[0] << " " << _pos[1] << " " << _pos[2] << endl;
        // cout << "z end barrel = " << GarlicGeometryParameters::Instance().Get_zOfBarrel() << endl;
        // cout << "z begin endcap = " << GarlicGeometryParameters::Instance().Get_zOfEndcap() << endl;
        //
        // cout << GarlicGeometryParameters::Instance().Get_nBarrelEcalLayers() << " " << rmax << endl;
        //
        // for (int i=0; i<GarlicGeometryParameters::Instance().Get_nBarrelEcalLayers(); i++) {
        //      cout << i << " " << GarlicGeometryParameters::Instance().Get_positionBarrelLayer() [i] << endl;
        // }
        // cout << "--------" << endl;
        assert(0);
      }

    }

  } else { // looks like endcap

    float endcap_start = GarlicGeometryParameters::Instance().Get_zOfEndcap();
    float endcap_end = GarlicGeometryParameters::Instance().Get_positionEndcapLayer() [GarlicGeometryParameters::Instance().Get_nEndcapEcalLayers()-1];

    endcap_start-=20;
    endcap_end+=20; // just to be generous...well actually what's taken from gear/XML file seems a not precisely what we expect...

    if (fabs(_pos[2])>=endcap_start &&
        fabs(_pos[2])<=endcap_end ) {
      if ( fabs(_pos[0])<GarlicGeometryParameters::Instance().Get_rInnerEcalEndcap() &&
           fabs(_pos[1])<GarlicGeometryParameters::Instance().Get_rInnerEcalEndcap() ) {
        zone = CALHITZONE_RING;
      } else
        zone = CALHITZONE_ENDCAP;
    } else {

      if (check) {
        cout << endl << "WARNING, hit outside ENDCAP/RING zone..." << endl <<
          " hit position " << _pos[0] << " " << _pos[1] << " " << _pos[2] << endl <<
          " barrel z = " << GarlicGeometryParameters::Instance().Get_zOfBarrel() << endl <<
          " endcap z = " << GarlicGeometryParameters::Instance().Get_zOfEndcap() << endl <<
          " my defn: begin, end of endcap " << endcap_start << " " << endcap_end << endl;
        assert(0);
      }

    }
  }

  if (zone==CALHITZONE_UNKNOWN) {

    if (check) {
      cout << "ERROR could not determine hits detector zone!!" << endl;
      cout << " position " << _pos[0] << " " << _pos[1] << " " << _pos[2] << endl;
      cout << " z of barrel = " << GarlicGeometryParameters::Instance().Get_zOfBarrel() <<
        " z of endcap = " << GarlicGeometryParameters::Instance().Get_zOfEndcap() << " my boundary " << bar_end_boundary << endl;
      assert(0);
    }

  }

  _zone = zone;
  return;
}


void PointInCalo::assignPseudoStave() {
  _pseudoStave = GarlicGeometryHelpers::getPseudoStave(_pos);
  return;
}

int PointInCalo::getBarrelPseudoLayer(const float* pos) {
  if (getZone()!=CALHITZONE_BARREL) {
    cout << "error from PointInCalo::getBarrelPseudoLayer, trying to get barrel info about point not in barrel! " << getZone() << endl;
    cout << "position " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
    return -999;
  }
  float x = pos[0];
  float y = pos[1];
  int bestLayer(-999);
  double barrelLayerDist = 9999.;
  double prodRadiusMax=0.;
  for(int istave = 0; istave < GarlicGeometryParameters::Instance().Get_symmetry(); ++istave){
    double prodRadius =
      x*GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][0]+
      y*GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][1];
    prodRadius/=sqrt( pow(GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][0],2) +
                      pow(GarlicGeometryParameters::Instance().Get_barrelStaveDir()[istave][1],2) );
    if(prodRadius>prodRadiusMax) prodRadiusMax=prodRadius;
  }


  float laythick = GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[1] - GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[0];

  if (prodRadiusMax>=GarlicGeometryParameters::Instance().Get_rOfBarrel() - laythick ) { // hit inside ECAL
    for(int ilayer=0; ilayer<=GarlicGeometryParameters::Instance().Get_nPseudoLayers(); ++ilayer){
      if(fabs(prodRadiusMax-GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[ilayer])<barrelLayerDist){
        barrelLayerDist = fabs(prodRadiusMax-GarlicGeometryParameters::Instance().Get_positionBarrelLayer()[ilayer]);
        bestLayer = ilayer;
      }
    }
  } else cout << "hit in front of ECAL..." << endl;
  if (bestLayer<0) cout << "error strange pseudolayer in barrel! " << bestLayer << endl;


  return bestLayer;
}

int PointInCalo::getEndcapPseudoLayer(const float* pos) {
  if (getZone()!=CALHITZONE_ENDCAP && getZone()!=CALHITZONE_RING) {
    cout << "error, trying to get endcap/ring info about point not in endcap/ring! " << getZone() << endl;
    return -999;
  }

  float z = fabs(pos[2]);
  int bestLayer(-999);
  double endcapLayerDist = 9999.;

  //  cout << " -- " << z << endl;

  for(int ilayer = 0; ilayer <= GarlicGeometryParameters::Instance().Get_nPseudoLayers(); ++ilayer){
    float dist = fabs( z-GarlicGeometryParameters::Instance().Get_positionEndcapLayer()[ilayer] );

    //    cout << "endcap layer " << ilayer << "  z position " << GarlicGeometryParameters::Instance().Get_positionEndcapLayer()[ilayer] << " " << dist << endl;

    if (dist<endcapLayerDist) {
      endcapLayerDist = dist;
      bestLayer = ilayer;
    }
  }
  return bestLayer;
}


void PointInCalo::assignPseudoLayer() {
  int bestLayer(-999);
  if (getZone()==CALHITZONE_BARREL) {
    bestLayer = getBarrelPseudoLayer();
  } else if (getZone()==CALHITZONE_ENDCAP || getZone()==CALHITZONE_RING) {
    bestLayer = getEndcapPseudoLayer();
  } else {
    //    cout << " strange zone? " << getZone() << endl;
  }

  if (bestLayer<0 || bestLayer>GarlicGeometryParameters::Instance().Get_nPseudoLayers()) {
    //    cout << "ERROR problem with pseudolayer assignment " << bestLayer << " " << getZone() << endl;
  } else {
    _pseudoLayer = bestLayer;
  }

  return;
}


void PointInCalo::movePseudoLayerAlongDirection(int dps, const float* dir, float* newpos) {
  // from the point, progress in a given direction for a certain number of pseudolayers

  float newpoint[3]={-99999,-99999,-99999};

  // get the new pseudolayer (present one plus any change)
  int newPLayer = getPseudoLayer()+dps;

  if (getZone()==CALHITZONE_BARREL) { // original point in barrel

    // second point on the direction line
    float point2[3]; sum(_pos, dir, point2);

    // intersection of line with requested p-layer
    getLineBarrelPseudoLayerIntersection(_pos, point2, newPLayer, newpoint);

    // check that this point is still in the barrel
    if ( fabs(newpoint[2])>GarlicGeometryParameters::Instance().Get_zOfBarrel() ) { // new position not inside barrel

      float displacement[3];
      diff(newpoint, _pos, displacement);

      // find where crosses end of barrel
      float sc = ( GarlicGeometryParameters::Instance().Get_zOfBarrel() - fabs(_pos[2]) )/ fabs(displacement[2]);

      float dispScale[3];
      scale(displacement, sc, dispScale);

      float barrelExit[3];
      sum(_pos, dispScale, barrelExit);
      PointInCalo exitPoint(barrelExit);
      int exitPLayer = exitPoint.getPseudoLayer(); // this is barrel pseudolayer at exit

      // how many layers to advance into endcap?
      int endcapLayer = newPLayer - exitPLayer;

      // z position of this endcap layer
      if (endcapLayer>=0 && endcapLayer<GarlicGeometryParameters::Instance().Get_nEndcapEcalLayers()) {
        float newz = GarlicGeometryParameters::Instance().Get_positionEndcapLayer() [endcapLayer];
        sc = (newz - fabs(_pos[2]))/displacement[2];
        scale(displacement, sc, dispScale);
        sum(_pos, dispScale, newpoint);
      }
    }

  } else if (getZone()==CALHITZONE_ENDCAP || getZone()==CALHITZONE_OVERLAP) { // start point in endcap
    int newlayer = getPseudoLayer()+dps;
    if (newlayer>=0 && newlayer<GarlicGeometryParameters::Instance().Get_nEndcapEcalLayers()) {
      float newzpos = GarlicGeometryParameters::Instance().Get_positionEndcapLayer() [newlayer];
      float sc = (newzpos - fabs(_pos[2]))/dir[2];
      float displ[3];
      scale(dir, sc, displ);
      float newendpos[3];
      sum(_pos, displ, newendpos);
      // check this new point is still in ECAL
      // first double-check the z coordinate
      if ( fabs(newendpos[2]) <  GarlicGeometryParameters::Instance().Get_positionEndcapLayer() [0] ||
           fabs(newendpos[2]) >= GarlicGeometryParameters::Instance().Get_positionEndcapLayer() [GarlicGeometryParameters::Instance().Get_nEndcapEcalLayers()-1] ) {
        cout << "error: this should never happen, point outside endcap ecal in z" << endl;
      } else { // check also the radial coordinate
        if (fabs(newendpos[0]) > GarlicGeometryParameters::Instance().Get_rOfEndcap() ) {
          cout << "warning: endcap point has too large a radius" << endl;
        } else { // it seems OK
          for (int i=0; i<3; i++)
            newpoint[i]=newendpos[i];
        }
      }
    } else if (newlayer<0 && (dir[2]*_pos[2])<0 ) { // layer in front of endcap, direction points towards centre: check if new point could be in barrel
      // get point where line hits end of barrel
      float sc = ( fabs(_pos[2])-(GarlicGeometryParameters::Instance().Get_zOfBarrel()-0.1) ) / ( _pos[2]-dir[2] );

      float barrelHitPoint[3];
      float displ[3];
      scale(dir, sc, displ);
      sum(_pos, displ, barrelHitPoint);

      PointInCalo barrelEntry(barrelHitPoint);

      float barrelpos[3];
      barrelEntry.movePseudoLayerAlongDirection(-newlayer, dir, barrelpos);

      if (barrelpos[0]>-99999)
        for (int i=0; i<3; i++)
          newpoint[i]=barrelpos[i];
    }
  } else {
    cout << "error, cannot move along direction for point in unknown zone" << endl;
  }

  for (int i=0; i<3; i++)
    newpos[i]=newpoint[i];

  return;
}

void GarlicExtendedHit::getRegion() {
  if (_layer>=0) return;
  _layer = (*GarlicGeometryParameters::Instance().Get_defaultDecoder())(_hit)["K-1"];
  _stave = (*GarlicGeometryParameters::Instance().Get_defaultDecoder())(_hit)["S-1"];
  _module = (*GarlicGeometryParameters::Instance().Get_defaultDecoder())(_hit)["M"];
  return;
}

float GarlicExtendedHit::getEnergyMIP() {

  if ( _energyMIP>0 ) return _energyMIP; // already calculated

  float origEnergy = getCaloHit() ? getCaloHit()->getEnergy() : -999; // this is the E-M energy
  if ( origEnergy<=0 ) return 0;

  // get the reference absorber thickness
  if ( _thick0<0 ) {
    int i(0);
    while ( _thick0<=0 && i<5 ) { // in case of preshower layers and so on...
      _thick0 = GarlicGeometryParameters::Instance().Get_absThicknessBarrelLayer()[i++];
    }
    if ( _thick0<=0 ) cout << "ERROR getting thickness of first ECAL absorber layer!" << endl;
  }

  // get the absorber thickness for this layer
  int layer = getLayer();
  float thickness(0);
  if ( getZone()==CALHITZONE_BARREL )
    thickness = GarlicGeometryParameters::Instance().Get_absThicknessBarrelLayer()[layer];
  else if ( getZone()==CALHITZONE_ENDCAP || getZone()==CALHITZONE_RING )
    thickness = GarlicGeometryParameters::Instance().Get_absThicknessEndcapLayer()[layer];
  else
    cout << "WEIRD, trying to get MIP energy of hit in strange detector region: " << getZone() << endl;

  float mipEnergyConv = GarlicAlgorithmParameters::Instance().GetEnergyMIPconversion();
  if ( mipEnergyConv <= 0 ) cout << "WARNING, EnergyMIPconversion in GarlicAlgorithmParameters does not seem to have been set..." << endl;

  float calib = (thickness/_thick0)*mipEnergyConv; // at normal incidence
  _energyMIP=calib*origEnergy;

  // angular correcion (path length in silicon)
  float mod(0);
  for (int i=0; i<3; i++) mod+=pow( getPosition()[i], 2);
  mod=sqrt(mod);
  float costh = fabs(getPosition()[2]/mod);
  float th = acos(costh);
  if ( getZone()==CLUS2_LOCATION_BARREL )
    _energyMIP*=sin(th);
  else if ( getZone()==CLUS2_LOCATION_ENDCAP )
    _energyMIP*=costh;



  //  cout << "thick: " << thickness << " thick0: " << _thick0 << " calib: " << calib << " orig: " << origEnergy << " MIP: " << _energyMIP << " costh: " << costh << endl;

  return _energyMIP;
}


