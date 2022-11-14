#ifndef CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H
#define CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H
#include "ChiMath/chi_math.h"
#include "ChiMath/Quadratures/quadrature.h"
#include "ChiMath/Quadratures/product_quadrature.h"
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "ChiMath/chi_math_04_gaussian_elimination_pivoting.h"
namespace chi_math
{
  class ProductQuadratureOp;
}
class chi_math::ProductQuadratureOp : public chi_math::ProductQuadrature
{
protected:
  chi_math::ProductQuadrature quad;
  const int method;
  int sn;
private:
void CheckInputs();
public:
  ProductQuadratureOp(const ProductQuadrature& inquad, int inmethod,
                      int order);

  double InnerProduct(const VecDbl& f, const VecDbl& g, const VecDbl& wt);
  void MakeHarmonicIndices();
  void BuildDiscreteToMomentOperator();
  void BuildMomentToDiscreteOperator();
  // Had to match up to dr. morel's symmetry;
  void OptimizeForPolarSymmetry(double normalization) override;
};
#endif //CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H