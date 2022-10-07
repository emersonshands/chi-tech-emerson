#include <algorithm>
#include <cmath>
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include <iomanip>
#include "chi_runtime.h"
#include "angular_quadrature_triangle.h"
#include "chi_log.h"
#include "ChiMath/chi_math.h"

chi_math::AngularQuadratureTriangle::
AngularQuadratureTriangle(unsigned int sn_in,
                          unsigned int in_method) :
  sn(sn_in),
  method(in_method)
{
  TriangleInit(sn);
}

void chi_math::AngularQuadratureTriangle::
MakeHarmonicIndices(unsigned int l_max)
{
  const int L = static_cast<int>(l_max);
  //Figure out what to do in dimensions other than 2d
  if(method == 1)
  {
    if (m_to_ell_em_map.empty())
    {
      for (int ell = 0; ell <= L; ++ell)
        for (int m = -ell; m <= ell; m += 2)
          if (ell == L and m > 0) break;
          else m_to_ell_em_map.emplace_back(ell, m);
    }
  }
}

void chi_math::AngularQuadratureTriangle::
TriangleInit(unsigned int sn)
{
  //Need to form the harmonics first/ these change based on the method
  //Clear the old m_to_ell mapping and redo it based on the method
  m_to_ell_em_map.clear();
  MakeHarmonicIndices(sn);
  chi_mesh::Vector3 new_omega;
  // grab the gauss points for the z axis, the number of points for GL is twice
  // that of the quadrature
  const auto old_omega = chi_math::QuadratureGaussLegendre(sn);
  // formulate the triangular quadrature
  int num_div = 1;
  for(auto u : old_omega.qpoints){
    double deltaVPhi = M_PI/(2.0*(double)num_div);
    //By polar symmetry we want to remove any negative polar angles
    if (u.x < 0) continue;
    for(int v=0; v<num_div;++v){
      double phi = deltaVPhi/2.0 + (double)v*deltaVPhi;
      double theta = acos(u.x);
      new_omega.x = sin(theta)*cos(phi);
      new_omega.y = sin(theta)*sin(phi);
      new_omega.z = cos(theta);
      weights.push_back(old_omega.weights[num_div-1]);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);
    }
    num_div++;
  }
  //This will loop through the other 3 parts of the unit circle
  //The order is x,y(done above); -x,y; -x,-y; x,-y
  double xsign = -1.0;
  double ysign = 1.0;
  size_t size = weights.size();
  for(int k=1;k<=3;++k){
    if(k>1){
      ysign=-1.0;
      if(k>2) xsign =1.0;
    }
    for(size_t l=0;l<size;++l){
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
}

void chi_math::AngularQuadratureTriangle::
BuildDiscreteToMomentOperator()
{
  //Build the d2m
  if (not d2m_op_built and method==1)
  {
    d2m_op.clear();
    unsigned int l_max = sn;
    MakeHarmonicIndices(l_max);
    unsigned int num_angles = abscissae.size();
    unsigned int num_moms = m_to_ell_em_map.size();
    for (const auto &ell_em: m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);

      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        //We need to change this ylm based on what method we use
        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        chi::log.Log0() << value << " Here is values "<< std::endl;
        double w = weights[n];
        cur_mom.push_back(value*w);
      }
      d2m_op.push_back(cur_mom);
    }
    d2m_op_built = true;
  }
}

void chi_math::AngularQuadratureTriangle::
BuildMomentToDiscreteOperator()
{
  if (not d2m_op_built) BuildDiscreteToMomentOperator();
  //method 1 is just the inverse of d2m
  if (method == 1)m2d_op = chi_math::Inverse(d2m_op);
  m2d_op_built = true;
}