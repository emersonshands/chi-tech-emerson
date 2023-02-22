#include "chi_math_06_condition.h"
//#include <slepcsvd.h>

void chi_math::condition(const MatDbl &A)
{
//  //Create the PETSC matrix
//  Mat B;
//  SVD svd;
//  SlepcIntialize(0,NULL,NULL,NULL);
//  //Allocate size based on the pointer to A
//  MatCreateSeqDense(PETSC_COMM_WORLD,A.size(), A[0].size(),PETSC_NULL,&B);
//  // Let PETSC decide on the local allocation
//  MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,
//              A.size(), A[0].size());
//  //Print out the matrix
//  MatView(B,PETSC_VIEWER_STDOUT_WORLD);
//  //Get PETSC ints for the global sizes
//  PetscInt Row;
//  PetscInt Col;
//  MatGetSize(B,&Row,&Col);
//  //Fill the petsc matrix then build it
//  for(int k=0;k<Row;++k)
//    for(int j=0;j<Col;++j)
//      MatSetValue(B,k,j,A[k][j],INSERT_VALUES);
//  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
//  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
//
//  SVDCreate(PETSC_COMM_WORLD,&svd);
//  SVDSetOperators(svd,B,NULL);
//
//  SVDSolve(svd);
//
//  PetscInt nsv;
//  SVDGetNumSingularValues(svd, &nsv);
//
//  PetscScalar* singularValues;
//  PetscMalloc1(nsv, &singularValues);
//  SVDGetSingularValues(svd, singularValues);
//  // Print the singular values
//
//  PetscPrintf(PETSC_COMM_SELF, "Singular values:\n");
//
//  for (int i = 0; i < nsv; i++) {
//
//    PetscPrintf(PETSC_COMM_SELF, "%f\n", singularValues[i]);
//
//  }
//  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
//  SVDConvergedReasonView(svd,PETSC_VIEWER_STDOUT_WORLD);
//  SVDGetConverged(SVD svd, int *nconv);
//  SVDGetSingularTriplet(SVD svd,int i,PetscReal *sigma,Vec u,Vec v);
//  PetscPrintf(PETSC_COMM_WORLD, "The singular value %f\n", sigma);
//  SVDDestroy(&svd);
//  MatDestroy(&B);
//  SlepcFinalize();
//
//#################################################
//  //Setting the options for the chi-tech KSP macro
//  std::string in_solver_name = "SVD_SOLVER";
//  std::string in_solver_type = "KSPGMRES";
//  //Precondition must be svd
//  std::string in_preconditioner_type = "svd";
//  //Arbitrary
//  double in_relative_residual_tolerance = 1e-8;
//  int64_t in_maximum_iterations = 300;
//  // Use the macro
//  auto SVD_SOLVER = chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
//    B,
//  in_solver_name,
//  in_solver_type,
//  in_preconditioner_type,
//  in_relative_residual_tolerance,
//  in_maximum_iterations);
//  //This is a preconditioner option, it should tell the svd
//  // to print out extreme eigen values
//  PCSetOptionsPrefix(SVD_SOLVER.pc,"pc_svd_monitor");
//
//  //set up the extra vectors we dont care about to solve for
//  Vec b;
//  VecSetSizes(b,PETSC_DECIDE,Row);
//  VecSet(b,0.0);
//  Vec x;
//  VecSetSizes(x,PETSC_DECIDE,Row);
//  VecSet(x,0.0);
//  //Solve for 0 vector with 0 vector matrix
//  //This should force the preconditioning and simply solve for 0
//  KSPSolve(SVD_SOLVER.ksp,b,x);
//  //De-allocate the resources
//  MatDestroy(&B);
//  VecDestroy(&b);
//  VecDestroy(&x);
//  KSPDestroy(&SVD_SOLVER.ksp);


//    //Grab the first maximum eigen value
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