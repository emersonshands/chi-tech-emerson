
#ifndef CHITECHMATRIXWORK_CHI_MATH_06_CONDITION_H
#define CHITECHMATRIXWORK_CHI_MATH_06_CONDITION_H
#include "chi_math.h"
#include "PETScUtils/petsc_utils.h"
#include "petsc.h"
namespace chi_math
{
  void condition(const MatDbl &A);
}
#endif //CHITECHMATRIXWORK_CHI_MATH_06_CONDITION_H