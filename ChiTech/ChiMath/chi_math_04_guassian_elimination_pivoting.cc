#include "ChiMath/chi_math.h"
#include "chi_math_04_Matrix_operations.h"
#include <cassert>
//###################################################################
/**Given an A matrix and b matrix, uses gaussian elimination
 * with semi-pivot to solve for Ax=b */
VecDbl chi_math::GaussEliminationPivot(const MatDbl& A,
                             const VecDbl& b)
{
  //check that the sizes are the same.
  assert(!A[0].empty());
  assert(A[0].size() == A[1].size());
  assert(!b.empty());
  assert(A[0].size() == b.size());

  int N = static_cast<int>(A[0].size());
  VecDbl x(N,0.0);

  //form the hold matrix
  MatDbl Q;
  for(int i=0;i<N;++i)
  {
    VecDbl tempA = A[i];
    tempA.push_back(b[i]);
    Q.emplace_back(tempA);
  }
  //Checking Q uncomment for verbose
  /**
  chi_math::PrintMatrix(Q);
   */
  //going down the diagonal
    //Look for best pivot
  for (int k = 0; k < N; k++)
  {
    int p = k;
    double Pivot = 0.0;
    for (int j = k; j < N; j++)
    {
      if (fabs(Q[j][k]) > Pivot)
      {
        Pivot = fabs(Q[j][k]);
        p = j;
      }
    }
    //swap
    if (p != k)
    {
      for (int i = k; i < N + 1; i++)
        std::swap(Q[p][i], Q[k][i]);
    }
    //factor out
    for (int i = k + 1; i < N; i++)
    {
      double m = Q[i][k] / Q[k][k];
      for (int j = k; j < N + 1; j++)
      {
        Q[i][j] -= m * Q[k][j];
      }
    }
    //Checking Q uncomment for verbose
    /**
     * chi_math::PrintMatrix(Q);
     */
  } // Going down diagonal
  //Back sub for x
  for (int f = N-1; f >=0; --f)
  {
    double b_val = Q[f][N];
    for (int j = f+1; j < N; j++)
      b_val -= Q[f][j] * x[j];
    x[f] = b_val / Q[f][f];
    /** Verbose check on x
      chi_math::PrintVector(x);
    */
  }
  return x;
}
///###################################################################
/**Given a matrix solve for the L1 norm */
double chi_math::l1Norm(const MatDbl& ArrayGiven)
{
  double condition = -1.0;
  const size_t rows = ArrayGiven[0].size();
  const size_t cols = ArrayGiven.size();
  for (size_t j = 0; j<cols;++j)
  {
    double sum = 0.0;
    for (size_t i=0; i<rows;++i)
    {
      sum+=abs(ArrayGiven[i][j]);
    }
    if (sum>condition) condition = sum;
  }
  return condition;
}
///###################################################################
/**Given a matrix solve for the L infinity norm */
double chi_math::linfNorm(const MatDbl& ArrayGiven)
{
  double condition = -1.0;
  const size_t rows = ArrayGiven[0].size();
  const size_t cols = ArrayGiven.size();
  for (size_t i=0; i<rows;++i)
  {
    double sum = 0.0;
    for (size_t j = 0; j<cols;++j)
    {
      sum+=abs(ArrayGiven[i][j]);
    }
    if (sum>condition) condition = sum;
  }
  return condition;
}