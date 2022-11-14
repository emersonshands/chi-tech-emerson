#ifndef CHITECHMATRIXWORK_CHI_MATH_04_GAUSSIAN_ELIMINATION_PIVOTING_H
#define CHITECHMATRIXWORK_CHI_MATH_04_GAUSSIAN_ELIMINATION_PIVOTING_H
#include <ChiMath/chi_math.h>
namespace chi_math{
  //04 Solver for Gaussian Elimination using pivoting
  VecDbl GaussEliminationPivot(const MatDbl& A,
                               const VecDbl& b);
}
#endif //CHITECHMATRIXWORK_CHI_MATH_04_GAUSSIAN_ELIMINATION_PIVOTING_H
