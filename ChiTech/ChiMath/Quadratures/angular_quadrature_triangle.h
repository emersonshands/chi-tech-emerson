#ifndef ANGULAR_QUADRATURE_TRIANGLE_H
#define ANGULAR_QUADRATURE_TRIANGLE_H

#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/chi_math_04_Matrix_operations.h"

namespace chi_math
{
  class AngularQuadratureTriangle;
}
/** Derived class for this specific special quadrature */
class chi_math::AngularQuadratureTriangle :
  public chi_math::AngularQuadrature
{
protected:
  const unsigned int method;
  const unsigned int sn;
  const unsigned int moments;
public:
  explicit
  AngularQuadratureTriangle(unsigned int in_method,unsigned int in_sn);
  explicit
  AngularQuadratureTriangle(unsigned int in_method,unsigned int in_sn,
                            unsigned int in_moments);
protected:
  void TriangleInit();
  void MakeHarmonicIndices(unsigned int scattering_order, int dimension)
  override;

public:
  void BuildDiscreteToMomentOperator(unsigned int scattering_order,
                                     int dimension) override;
  void BuildMomentToDiscreteOperator(unsigned int scattering_order,
                                     int dimension) override;

  void FilterMoments(unsigned int scattering_order);
};

#endif //ANGULAR_QUADRATURE_TRIANGLE_H