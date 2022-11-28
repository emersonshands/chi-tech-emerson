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
                          unsigned int sn_in) :
  method(in_method),
  sn(sn_in),
  moments(0)
{
  TriangleInit(sn);
}

chi_math::AngularQuadratureTriangle::
AngularQuadratureTriangle(unsigned int in_method,
                          unsigned int sn_in,
                          unsigned int inmoments) :
  method(in_method),
  sn(sn_in),
  moments(inmoments)
{
  TriangleInit(sn);
}

void chi_math::AngularQuadratureTriangle::
MakeHarmonicIndices(unsigned int l_max)
{
  const int L = static_cast<int>(l_max);
  //Figure out what to do in dimensions other than 2d
  if (m_to_ell_em_map.empty())
  {
    for (int ell = 0; ell <= L; ++ell)
      for (int m = -ell; m <= ell; m += 2)
        if (ell == L and m >= 0) break;
        else m_to_ell_em_map.emplace_back(ell, m);
  }
}

void chi_math::AngularQuadratureTriangle::
TriangleInit(unsigned int sn)
{
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
  //Need to form the harmonics first/ these change based on the method
  //Clear the old m_to_ell mapping and redo it based on the method
  m_to_ell_em_map.clear();
  MakeHarmonicIndices(sn);
  chi_mesh::Vector3 new_omega;
  // grab the gauss points for the z axis, the number of points for GL is twice
  // that of the quadrature
  const auto old_omega = chi_math::QuadratureGaussLegendre(sn);
  // formulate the triangular quadrature
  //chi_math::PrintVector(old_omega.weights);

  int num_div = 1;
  int weightPos =0;
  for(auto u : old_omega.qpoints)
  {
    double deltaVPhi = M_PI/(2.0*(double)num_div);
    //By polar symmetry we want to remove any negative polar angles
    if (u.x <= 0.0)
    {
      weightPos++;
      continue;
    }
    for(int v=0; v<num_div;++v)
    {
      double phi = deltaVPhi/2.0 + (double)v*deltaVPhi;
      double theta = acos(u.x);
      new_omega.x = sin(theta)*cos(phi);
      new_omega.y = sin(theta)*sin(phi);
      new_omega.z = cos(theta);
      weights.push_back(old_omega.weights[weightPos]);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);

    }
    weightPos++;
    num_div++;
  }
  //This will loop through the other 3 parts of the unit circle
  //The order is x,y(done above); -x,y; -x,-y; x,-y
  double xsign = -1.0;
  double ysign = 1.0;
  size_t size = weights.size();
  for(int k=1;k<=3;++k)
  {
    if(k>1)
    {
      ysign=-1.0;
      if(k>2) xsign =1.0;
    }
    for(size_t l=0;l<size;++l)
    {
      double phi = abscissae[l].phi+k*(M_PI/2.0);
      double theta = abscissae[l].theta;
      new_omega.x = omegas[l].x*xsign;
      new_omega.y = omegas[l].y*ysign;
      new_omega.z = omegas[l].z;
      weights.push_back(weights[l]);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);

    }
  }
  //Now we need to call optimize for polar symmetry to normalize
  // the weights to 4pi to correctly integrate over the sphere
  chi_math::AngularQuadrature::OptimizeForPolarSymmetry(4.0*M_PI);
}