#ifndef ECALGARLICGEOMETRYHELPER_HH_
#define ECALGARLICGEOMETRYHELPER_HH_


class ECALGarlicGeometryHelpers {

public:

  ECALGarlicGeometryHelpers () {}

  ~ECALGarlicGeometryHelpers() {}

public:

  double Get2dProjDistance(const float* a, const float* b);
  double Get3dDistance(const float* a, const float* b);
  void getFrontFaceProj(const float* a, float* vec);

  void GetGeneralPointPlaneProjection(const float* point, const float* pointOnPlane, const float* planeNormal, const float* projDirection, float* projection);

  double GetDistToLine(const float* a, const float* l1, const float* l2);

protected:

  void   cross(const float* a, const float* b, float* vec);
  float  dot(const float* a, const float* b);
  float  mag(const float* a);
  void   diff(const float* a, const float* b, float* vec);
  void   sum(const float* a, const float* b, float* vec);
  float  angle(const float* a, const float* b);
  float  cosangle(const float* a, const float* b);
  float  phi(const float* a);
  float  costheta(const float* a);
  float  theta(const float* a);
  void   norm(const float* a, float* norma);
  void   scale(const float* a, float factor, float* aprime);

  float getPseudoLayerRadius(int psLayer);
  float getPseudoLayerDistNormal(int ps1, int ps2);

  int   getPseudoStave(const float* pos);
  float getPseudoStavePhi(int istave);

  float det(const float* abcd);

  void  get3dLinePlaneIntersection(const float* planePoint, const float* planeNormal, const float* linePoint, const float* lineDirection, float* intersection);
  void  get2dLineIntersection(const float* l1a, const float* l1b, const float* l2a, const float* l2b, float* intersection);
  void  getPointsOnBarrelPseudoLayer(int pslayer, int pseudostave, float* point1, float* point2);
  void  getLineBarrelPseudoLayerIntersection(const float* point1, const float* point2, int pseudolayer, float* intersection);

};


#endif
