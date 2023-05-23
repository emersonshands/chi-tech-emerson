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

  if (method != 1 and method != 2 and method !=3 and method!=0)
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
  VecDbl newZi;
  VecDbl newWeights;
  for(size_t pos =0;pos<old_omega.qpoints.size();++pos)
  {
    if (old_omega.qpoints[pos].x < 0) continue;
    newZi.push_back(old_omega.qpoints[pos].x);
    newWeights.push_back(old_omega.weights[pos]);
  }
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
      chi::log.Log0() << std::setprecision(16) << " Z COSINE " << new_z_value;
      double phi = deltaVPhi/2.0 + (double)v*deltaVPhi;
      double theta = acos(new_z_value);
      double sinTheta = sqrt(1-new_z_value*new_z_value);
      new_omega.x = sinTheta*cos(phi);
      new_omega.y = sinTheta*sin(phi);
      new_omega.z = new_z_value;
      double weightCurrent = old_omega.weights[weightPos]/(num_div);
      weights.push_back(old_omega.weights[weightPos]*(M_PI/(2.0*num_div)));
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

  //This is the number of points in 1 octant
  size_t octSize = weights.size();
  for(int octant=1;octant<=3;++octant)
  {
    //This is how much should be added to each phi to get the orientation  right of the first index
    double offset = M_PI_2*octant;
    for(size_t point=0;point<octSize;++point)
    {
      double phi = abscissae[point].phi+offset;
      double theta = abscissae[point].theta;
      double new_z_value = omegas[point].z;
      double sinTheta = sqrt(1-new_z_value*new_z_value);
      new_omega.x = sinTheta*cos(phi);//omegas[l].x*xsign;
      new_omega.y = sinTheta*sin(phi); //omegas[l].y*ysign;
      new_omega.z = omegas[point].z;
      weights.push_back(weights[point]);
      omegas.emplace_back(new_omega);
      abscissae.emplace_back(phi,theta);
      chi::log.Log0() << " Z COSINE " << new_omega.z;
      chi::log.Log0()<< "Phi value original "<< phi << "Theta value "<<theta;
      chi::log.Log0()<< "OMEGA x "<< new_omega.x <<
                     " OMEGA Y "<< new_omega.y << " OMEGA Z " << new_omega.z;
      chi::log.Log0()<< "WEIGHT "<< weights.back();

    }
  }
//  //This builds the negative zi values
//  size_t hemisize = weights.size();
//  for(int position=0;position<hemisize;++position)
//  {
//
//    new_omega.x = omegas[position].x;//omegas[l].x*xsign;
//    new_omega.y = omegas[position].y; //omegas[l].y*ysign;
//    new_omega.z = -omegas[position].z;
//    double new_z_value = -omegas[position].z;
//    double phi = abscissae[position].phi;
//    double theta = acos(new_z_value);
//    weights.push_back(weights[position]);
//    omegas.emplace_back(new_omega);
//    abscissae.emplace_back(phi,theta);
//    chi::log.Log0() << " Z COSINE " << new_omega.z;
//    chi::log.Log0()<< "Phi value " << phi << " Theta value "<<theta;
//    chi::log.Log0()<< "OMEGA x "<< new_omega.x <<
//                   " OMEGA Y "<< new_omega.y << " OMEGA Z " << new_omega.z;
//    chi::log.Log0()<< "WEIGHT "<< weights.back();
//    chi::log.Log0()<< "OMEGA Z NEGATIVE: " << new_omega.z << "; VS COS OF NEW THETA:" << cos(theta);
//  }

  chi_math::AngularQuadrature::OptimizeForPolarSymmetry(4.0*M_PI);
  //Now we need to call optimize for polar symmetry to normalize
  // the weights to 4pi to correctly integrate over the sphere
//  if (method!=0){
//    chi_math::AngularQuadrature::OptimizeForPolarSymmetry(4.0*M_PI);
//  }
//  else
//  {
//    double normal = 4.0*M_PI;
//    std::vector<chi_math::QuadraturePointPhiTheta> new_abscissae;
//    std::vector<double>                            new_weights;
//    std::vector<chi_mesh::Vector3>                 new_omegas;
//
//    const size_t num_dirs = omegas.size();
//    double weight_sum = 0.0;
//    for (size_t d=0; d<num_dirs; ++d)
//      if (omegas[d].y > 0.0)
//      {
//        new_abscissae.emplace_back(abscissae[d]);
//        new_weights.emplace_back(weights[d]);
//        new_omegas.emplace_back(omegas[d]);
//        std::cout<<"Current Omega: x: "<< omegas[d].x <<
//        ", y: "<< omegas[d].y<< ", z: "<< omegas[d].z << std::endl;
//        std::cout<<"Current Angle: phi: "<< abscissae[d].phi <<
//                 ", theta: "<< abscissae[d].theta << std::endl;
//        weight_sum += weights[d];
//      }
//
//    if (normal > 0.0)
//      for (double& w : new_weights)
//        w *= normal/weight_sum;
//
//    abscissae = std::move(new_abscissae);
//    weights   = std::move(new_weights);
//    omegas    = std::move(new_omegas);
//  }
}