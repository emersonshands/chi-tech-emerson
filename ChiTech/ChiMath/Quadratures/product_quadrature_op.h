#ifndef CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H
#define CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H
#include "ChiMath/chi_math.h"
#include "ChiMath/Quadratures/quadrature.h"
#include "ChiMath/Quadratures/angular_product_quadrature.h"
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "ChiMath/chi_math_04_Matrix_operations.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"

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
  const int moments;
private:
void CheckInputs();
public:
  ProductQuadratureOp(const ProductQuadrature &inquad,
                     int inmethod, int inorder);
  ProductQuadratureOp(const ProductQuadrature& inquad,
                      int inmethod, int inorder, int inmoments);

  double InnerProduct(const VecDbl& f, const VecDbl& g, const VecDbl& wt);
  void MakeHarmonicIndices(unsigned int scattering_order, int dimension) override;
  void BuildDiscreteToMomentOperator(unsigned int scattering_order,
                                     int dimension) override;
  void BuildMomentToDiscreteOperator(unsigned int scattering_order,
                                     int dimension) override;
//  void OptimizeForPolarSymmetry(double normalization) override;
  void FilterMoments(unsigned int scattering_order);
};
#endif //CHITECHMATRIXWORK_PRODUCT_QUADRATURE_OP_H