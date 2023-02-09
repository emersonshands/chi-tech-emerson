#include "chi_math_06_condition.h"

void chi_math::condition(const MatDbl &A)
{
  //Grab the first maximum eigen value
  //Solve for eigen values
  MatDbl D = chi_math::MatMul(A,chi_math::Transpose(A));
  size_t R = D.size();
  size_t C = D[0].size();
//  chi_math::PrintMatrix(D);
  VecDbl eVec_temp(R, 1);
  auto eigenMax = chi_math::PowerIteration(D, eVec_temp, 2000, 1.0e-16);

  //Spectral shift to get the minimum
  MatDbl I(R, VecDbl(C, 0.0));
  for (size_t j = 0; j < R; ++j)
    for (size_t i = 0; i < C; ++i)
      if (i == j)
        I[j][i] = 1.0;
//  chi_math::PrintMatrix(I);

  //Subtract the eigen*I from A to get B and solve for
  MatDbl B = chi_math::MatSubtract(chi_math::MatMul(I,eigenMax),D);
//  chi_math::PrintMatrix(B);
  auto eigenMin = chi_math::PowerIteration(B, eVec_temp, 2000, 1.0e-16) + eigenMax;
  eigenMin =sqrt(eigenMin);
  eigenMax = sqrt(eigenMax);
  printf("Here is the max\n %f\nHere is the min\n%f\n",eigenMax,eigenMin);
  printf("THE CONDITION NUMBER: %f\n",eigenMax/eigenMin);
}