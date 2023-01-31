#include "legendrepoly.h"
#include <cmath>

#include <algorithm>
//###################################################################
/**Provides the function evaluation of the Legendre polynomial
 * P_N at value x.

 \param N int Order of the Legendre polynomial.
 \param x double The evaluation point.*/
double chi_math::Legendre(int N, double x)
{
  double Pnm1 = 1;
  double Pn   = x;
  double Pnp1 = 0;

  if (N==0) {return 1;}

  if (N==1) {return x;}

  for (int n=2;n<=N; n++)
  {
    int ns = n-1;
    Pnp1 = ((2.0*ns+1)/(ns+1.0))*x*Pn - (ns/(ns+1.0))*Pnm1;
    Pnm1 = Pn;
    Pn = Pnp1;
  }

  return Pnp1;
}

//###################################################################
/**Provides the function evaluation of the derivative of Pn at value x

 \param N int Order of the Legendre polynomial.
 \param x double The evaluation point.*/
double chi_math::dLegendredx(int N, double x)
{
  if (N==0) {return 0;}

  if (N==1) {return 1;}

  double retval = (N*x/(x*x-1))*Legendre(N,x);
         retval-= (N/(x*x-1))*Legendre(N-1,x);

  return retval;
}

//###################################################################
/**Provides the function evaluation of the second derivative of Pn at value x

 \param N int Order of the Legendre polynomial.
 \param x double The evaluation point.*/
double chi_math::d2Legendredx2(int N, double x)
{
  double epsilon = 1.0e-8;
  if (N==0) {return 0.0;}

  if (N==1) {return 0.0;}

  double xpos = std::min(x+epsilon, 1.0-1.0e-10);
  double xneg = std::max(x-epsilon,-1.0+1.0e-10);
  double dx = xpos - xneg;

  double dPdx_pos = dLegendredx(N,xpos);
  double dPdx_neg = dLegendredx(N,xneg);

  return (dPdx_pos - dPdx_neg)/dx;
}
