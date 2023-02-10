#include "chi_math_06_condition.h"

void chi_math::condition(const MatDbl &A)
{

  Mat B;
  MatCreateSeqDense(PETSC_COMM_WORLD,A.size(), A[0].size(),PETSC_NULL,&B);
  MatSetType(B,MATMPIAIJ);
  MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,
              A.size(), A[0].size());

  MatMPIAIJSetPreallocation(B,1, nullptr,
                            0, nullptr);

  MatSetOption(B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  for(int i=0;i<A.size();++i)
    for (int j = 0; j < A[0].size(); ++j)
    {

      MatSetValue(B,)
    }
  std::string in_solver_name = "SVD_SOLVER";
  std::string in_solver_type = "KSPGMRES";
  std::string in_preconditioner_type = "svd";
  double in_relative_residual_tolerance = 1e-8;
  int64_t in_maximum_iterations = 300;
  //  -pc_type svd -pc_svd_monitor
  auto SVD_SOLVER = chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    B,
  in_solver_name,
  in_solver_type,
  in_preconditioner_type,
  in_relative_residual_tolerance,
  in_maximum_iterations);
  PCSetOptionsPrefix(SVD_SOLVER.pc,"pc_svd_monitor");

//  //Grab the first maximum eigen value
//  //Solve for eigen values
//  MatDbl D = chi_math::MatMul(A,chi_math::Transpose(A));
//  size_t R = D.size();
//  size_t C = D[0].size();
////  chi_math::PrintMatrix(D);
//  VecDbl eVec_temp(R, 1);
//  auto eigenMax = chi_math::PowerIteration(D, eVec_temp, 2000, 1.0e-16);
//
//  //Spectral shift to get the minimum
//  MatDbl I(R, VecDbl(C, 0.0));
//  for (size_t j = 0; j < R; ++j)
//    for (size_t i = 0; i < C; ++i)
//      if (i == j)
//        I[j][i] = 1.0;
////  chi_math::PrintMatrix(I);
//
//  //Subtract the eigen*I from A to get B and solve for
//  MatDbl B = chi_math::MatSubtract(chi_math::MatMul(I,eigenMax),D);
////  chi_math::PrintMatrix(B);
//  auto eigenMin = chi_math::PowerIteration(B, eVec_temp, 2000, 1.0e-16) + eigenMax;
//  eigenMin =sqrt(eigenMin);
//  eigenMax = sqrt(eigenMax);
//  printf("Here is the max\n %f\nHere is the min\n%f\n",eigenMax,eigenMin);
//  printf("THE CONDITION NUMBER: %f\n",eigenMax/eigenMin);
}