#ifndef CHITECHMATRIXWORK_CHI_MATH_04_MATRIX_OPERATIONS_H
#define CHITECHMATRIXWORK_CHI_MATH_04_MATRIX_OPERATIONS_H
#include <ChiMath/chi_math.h>
namespace chi_math{
  //04 Solver for Gaussian Elimination using pivoting
  VecDbl GaussEliminationPivot(const MatDbl& A,
                               const VecDbl& b);

  double l1Norm(const MatDbl& ArrayGiven);
  double linfNorm(const MatDbl& ArrayGiven);
}
#endif //CHITECHMATRIXWORK_CHI_MATH_04_MATRIX_OPERATIONS_H
