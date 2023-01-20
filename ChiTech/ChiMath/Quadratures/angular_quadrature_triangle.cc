#include <algorithm>
#include <cmath>
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include <iomanip>
#include "chi_runtime.h"
#include "angular_quadrature_triangle.h"
#include "chi_log.h"
#include "ChiMath/chi_math.h"

chi_math::AngularQuadratureTriangle::
AngularQuadratureTriangle(unsigned int in_method,
                          unsigned int in_sn) :
  method(in_method),
  sn(in_sn),
  moments(0)
{
  TriangleInit();
}

chi_math::AngularQuadratureTriangle::
AngularQuadratureTriangle(unsigned int in_method,
                          unsigned int in_sn,
                          unsigned int in_moments) :
  method(in_method),
  sn(in_sn),
  moments(in_moments)
{
  TriangleInit();
}

void chi_math::AngularQuadratureTriangle::
TriangleInit()
{

  chi::log.Log0() << "Given the method "
  << method << "\nGiven sn " << sn;

  if (method != 1 and method != 2 and method !=3)
  {
    printf("The method given is not 1, 2, or 3.\n Please reorder your input.\n "
           "Given value %i\n",method);
    chi::log.Log0Error() << "Mismatch in method given";
    chi::Exit(510);
  }
  if (moments<0 or moments>sn)
  {
    printf("The moments user is asking for is greater than or less "
           "than the allowable moments for sn\n"
           "Given value %i, sn=%i\n",moments,sn);
    chi::log.Log0Error() << "Mismatch between given moments to cut and sn used";
    chi::Exit(510);
  }
  if (moments == sn)
  {
    printf("The moments user is asking for is equal to the given sn\n"
           "All values used\n"
           "Given value %i, sn=%i\n",moments,sn);
  }
  chi_mesh::Vector3 new_omega;
  // grab the gauss points for the z axis, the number of points for GL is twice
  // that of the quadrature
  const auto old_omega = chi_math::QuadratureGaussLegendre(sn);
  // formulate the triangular quadrature

  int num_div = 1;
  int weightPos =0;
  for(auto u : old_omega.qpoints)
  {
    double deltaVPhi = M_PI/(2.0*(double)num_div);
    // When the QuadratureGaussLegendre gives us the x points we will
    // use for our z positions on the unit sphere, they are ordered in
    // descending order largest magnitude towards 0 on the negative side,
    // and ascending order on the positive side of 0
    // The weights are defined using the weights given by the GL quadrature
    // which is the position in the weights that weightPos keeps track of.
    // This will ignore the positive x values, and just use the descending
    // order given by the quadrature and use the absolute value of the x values.
    if (u.x >= 0.0)
    {
      break;
    }
    for(int v=0; v<num_div;++v)
    {
      double new_z_value = abs(u.x);
      double phi = deltaVPhi/2.0 + (double)v*deltaVPhi;
      double theta = acos(new_z_value);
      new_omega.x = sin(theta)*cos(phi);
      new_omega.y = sin(theta)*sin(phi);
      new_omega.z = cos(theta);
      weights.push_back(old_omega.weights[weightPos]/num_div);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);
      chi::log.Log0()<< "Phi value "<< phi << " Theta value "<<theta;
      chi::log.Log0()<< "OMEGA x "<< new_omega.x <<
      " OMEGA Y "<< new_omega.y << " OMEGA Z " << new_omega.z;
      chi::log.Log0()<< "WEIGHT "<< weights.back();
    }
    weightPos++;
    num_div++;
  }
  //This will loop through the other 3 parts of the unit circle
  //The order is x,y(done above); -x,y; -x,-y; x,-y
  double xsign = -1.0;
  double ysign = 1.0;
  size_t sizew = weights.size();
  for(int k=1;k<=3;++k)
  {
    if(k>1)
    {
      ysign=-1.0;
      if(k>2) xsign =1.0;
    }
    for(size_t l=0;l<sizew;++l)
    {
      double phi = abscissae[l].phi+k*(M_PI/2.0);
      double theta = abscissae[l].theta;
      new_omega.x = omegas[l].x*xsign;
      new_omega.y = omegas[l].y*ysign;
      new_omega.z = omegas[l].z;
      weights.push_back(weights[l]);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);
      chi::log.Log0()<< "Phi value "<< phi << " Theta value "<<theta;
      chi::log.Log0()<< "OMEGA x "<< new_omega.x <<
                     " OMEGA Y "<< new_omega.y << " OMEGA Z " << new_omega.z;
      chi::log.Log0()<< "WEIGHT "<< weights.back();
    }
  }
  //Now we need to call optimize for polar symmetry to normalize
  // the weights to 4pi to correctly integrate over the sphere
  chi_math::AngularQuadrature::OptimizeForPolarSymmetry(4.0*M_PI);
}
//--S4 Product JAN p0
//[0]  Balance table:
//[0]   Absorption rate          = 4.60886e-01
//[0]   Production rate          = 0.00000e+00
//[0]   In-flow rate             = 1.00000e+00
//[0]   Out-flow rate            = 5.39114e-01
//[0]   Integrated scalar flux   = 1.00000e+00
//[0]   Net Gain/Loss            = 1.73195e-12
//[0]   Net Gain/Loss normalized = 1.73195e-12
//                                 [0]
//                                 [0]  ********** Done computing balance
//[0]  XMax
//[0]  0.0072069542459793
//[0]  XMin
//[0]  0.33086004584623
//[0]  YMax
//[0]  0.15131901513752
//[0]  Ymin
//[0]  0.049727744282467


//--s2 product JAN
//[0]  Balance table:
//[0]   Absorption rate          = 4.28510e-01
//[0]   Production rate          = 0.00000e+00
//[0]   In-flow rate             = 1.00000e+00
//[0]   Out-flow rate            = 5.71490e-01
//[0]   Integrated scalar flux   = 1.00000e+00
//[0]   Net Gain/Loss            = -3.85219e-11
//[0]   Net Gain/Loss normalized = -3.85219e-11
//[0]
//[0]  ********** Done computing balance
//[0]  XMax
//[0]  0.0035282501384168
//[0]  XMin
//[0]  0.35647061065017
//[0]  YMax
//[0]  0.17163425902382
//[0]  Ymin
//[0]  0.039856762456773
//[0]


//p1 --JAN S4
//[0]  Balance table:
//[0]   Absorption rate          = 4.63607e-01
//[0]   Production rate          = 0.00000e+00
//[0]   In-flow rate             = 1.00000e+00
//[0]   Out-flow rate            = 5.36393e-01
//[0]   Integrated scalar flux   = 1.00000e+00
//[0]   Net Gain/Loss            = 3.09568e-11
//[0]   Net Gain/Loss normalized = 3.09568e-11
//                                 [0]
//                                 [0]  ********** Done computing balance
//[0]  XMax
//[0]  0.0060785969814323
//[0]  XMin
//[0]  0.32154187484683
//[0]  YMax
//[0]  0.14213211214103
//[0]  Ymin
//[0]  0.066640728697553


//P1 JAN S8
//[0]  Balance table:
//[0]   Absorption rate          = 4.72729e-01
//[0]   Production rate          = 0.00000e+00
//[0]   In-flow rate             = 1.00000e+00
//[0]   Out-flow rate            = 5.27271e-01
//[0]   Integrated scalar flux   = 1.00000e+00
//[0]   Net Gain/Loss            = 4.20567e-11
//[0]   Net Gain/Loss normalized = 4.20567e-11
//                                 [0]
//                                 [0]  ********** Done computing balance
//[0]  XMax
//[0]  0.0062171195870433
//[0]  XMin
//[0]  0.31083716971508
//[0]  YMax
//[0]  0.14123099123847
//[0]  Ymin
//[0]  0.068985327127354

p2 JAN S8
//[0]  Balance table:
//[0]   Absorption rate          = 4.74070e-01
//[0]   Production rate          = 0.00000e+00
//[0]   In-flow rate             = 1.00000e+00
//[0]   Out-flow rate            = 5.25930e-01
//[0]   Integrated scalar flux   = 1.00000e+00
//[0]   Net Gain/Loss            = 3.97460e-11
//[0]   Net Gain/Loss normalized = 3.97460e-11
//                                 [0]
//                                 [0]  ********** Done computing balance
//[0]  XMax
//[0]  0.0062173774560831
//[0]  XMin
//[0]  0.30671201423743
//[0]  YMax
//[0]  0.14208960239754
//[0]  Ymin
//[0]  0.070910848725567
