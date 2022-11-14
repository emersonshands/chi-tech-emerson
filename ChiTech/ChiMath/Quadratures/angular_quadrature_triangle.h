#ifndef ANGULAR_QUADRATURE_TRIANGLE_H
#define ANGULAR_QUADRATURE_TRIANGLE_H

#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/chi_math.h"
#include "ChiMath/chi_math_04_gaussian_elimination_pivoting.h"

namespace chi_math
{
  class AngularQuadratureTriangle;
}
/** Derived class for this specific special quadrature */
class chi_math::AngularQuadratureTriangle :
  public chi_math::AngularQuadrature
{
protected:
  const unsigned int sn;
  const unsigned int method;
public:
  explicit
  AngularQuadratureTriangle(unsigned int sn, unsigned int in_method);
protected:
  void TriangleInit(unsigned int sn);
  void MakeHarmonicIndices(unsigned int l_max);

public:
  void BuildDiscreteToMomentOperator();
  void BuildMomentToDiscreteOperator();
};

#endif //ANGULAR_QUADRATURE_TRIANGLE_H