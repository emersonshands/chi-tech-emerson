#include <algorithm>
#include <cmath>
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "angular_quadrature_triangle.h"
#include "ChiMath/chi_math.h"
#include <cassert>
#include <numeric>
#include <iomanip>

//##################################################
/** Local inner product function, only useful for method 3.
 * Not added into general chi_math, because no other method will use it.
 */
double InnerProduct(const VecDbl& f, const VecDbl& g, const VecDbl& wt)
{
  size_t fsize = f.size();
  double sum_val = 0.0;
  for (size_t i=0; i<fsize; ++i)
  {
    sum_val += f[i] * g[i] * wt[i];
  }
  return sum_val;
}

void chi_math::AngularQuadratureTriangle::FilterMoments(unsigned int scattering_order)
{
  if (m2d_op_built and d2m_op_built and moments!=0)
  {
    int s_order = static_cast<int>(scattering_order);
    int moments_to_keep = 1 +
      (s_order*3 + s_order*s_order)/2;
    auto m2d_transposed = m2d_op;
    chi::log.Log0() << "Size of m2ell preresize " << m_to_ell_em_map.size();
    m_to_ell_em_map.resize(moments_to_keep);
    chi::log.Log0() << "Size of m2ell postresize " << m_to_ell_em_map.size();
    std::vector<std::vector<double>> m2dworking;
    std::vector<std::vector<double>> d2mworking;
    for (size_t i {}; i<moments_to_keep;++i)
    {
      chi::log.Log0() << "L " << m_to_ell_em_map[i].ell << " M " << m_to_ell_em_map[i].m;
      d2mworking.push_back(d2m_op.at(i));
      m2dworking.push_back(m2d_transposed.at(i));
    }
    m2d_op = m2dworking;
    d2m_op = d2mworking;
  }
}

void chi_math::AngularQuadratureTriangle::
MakeHarmonicIndices(unsigned int scattering_order, int dimension)
{
  int L = static_cast<int>(sn);
  int L_max = static_cast<int>(scattering_order);
  if (method == 0 and m_to_ell_em_map.empty())
  {
    for (int ell=0; ell<=L; ell++)
      for (int m=-ell; m<=ell; m+=2)
      {
        m_to_ell_em_map.emplace_back(ell, m);
        chi::log.Log0() << "l " << ell << " and m "<< m << "\n";
      }
  }
  if (m_to_ell_em_map.empty() )
  {
    for (int ell = 0; ell <= L_max; ++ell)
      for (int m = -ell; m <= ell; m += 2)
      {
        if (ell == L and m >= 0 and ell != 0) break;
        else m_to_ell_em_map.emplace_back(ell, m);
        chi::log.Log0() << "l " << ell << " and m "<< m << "\n";
      }
  }

}

void chi_math::AngularQuadratureTriangle::
BuildDiscreteToMomentOperator
  (unsigned int scattering_order,
   int dimension)
{
  if (d2m_op_built) return;
  MakeHarmonicIndices(scattering_order,dimension);
  //Standard Sn method
  if (method == 0)
  {
    d2m_op.clear();
    const size_t num_angles = abscissae.size();
    const size_t num_moms = m_to_ell_em_map.size();
//    //####################################################
//    //To solve for the weights
//    chi::log.Log0() << "Solving Weights \n";
//    std::vector<std::vector<double>> cmt;
//    for (const auto &ell_em: m_to_ell_em_map)
//    {
//      std::vector<double> cur_mom;
//      cur_mom.reserve(num_angles);
//
//      for (int n = 0; n < num_angles; n++)
//      {
//        const auto &cur_angle = abscissae[n];
//        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
//                                     cur_angle.phi,
//                                     cur_angle.theta);
//        cur_mom.push_back(value);
//      }
//      cmt.push_back(cur_mom);
//    }
//    //     solve for the weights
//    std::vector<double> wt = {4.0*M_PI};
//    for (size_t i = 1; i < weights.size(); ++i)
//      wt.emplace_back(0.0);
//
//    chi::log.Log0() << "THE WT \n";
//    chi_math::PrintVector(wt);
////    chi::log.Log0() << "The original cmt \n";
////    chi_math::PrintMatrix(cmt);
////    chi::log.Log0() << "The inverse cmt \n";
////    chi_math::PrintMatrix(invt);
////    auto new_weights = chi_math::MatMul(invt, wt);
//    chi::log.Log0() << "The old weights \n";
//    chi_math::PrintVector(weights);
////    chi::log.Log0() << "The new weights \n";
////    chi_math::PrintVector(new_weights);
//    chi::log.Log0() << "NEW SOLVE \n";
//    std::vector<double> new_weights = chi_math::GaussEliminationPivot(cmt,wt);
//    chi_math::PrintVector(chi_math::GaussEliminationPivot(cmt,wt));
//    weights = new_weights;
//    //####################################################
    for (const auto& ell_em : m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      double integralYlm = 0.0;
      double integralProduct = 0.0;

      for (int n=0; n<num_angles; n++)
      {
        const auto& cur_angle = abscissae[n];
        double value = chi_math::Ylm(ell_em.ell,ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        double w = weights[n];
        cur_mom.push_back(value*w);
        integralProduct += value*value*w;
        integralYlm += value*w;
      }
      chi::log.Log0() << "$$$$$$$$$$\nHARMONIC l=" << ell_em.ell << " m=" << ell_em.m;
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value " << integralYlm;
      chi::log.Log0() << "Product Value " << integralProduct;
      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(integralProduct))/
                                              (4.0*M_PI/(2.0*ell_em.ell+1.0));

      d2m_op.push_back(cur_mom);
    }
    d2m_op_built = true;
  }
  if (method == 1)
  {
    d2m_op.clear();
    if (not m2d_op_built) BuildMomentToDiscreteOperator(scattering_order, dimension);
    chi::log.Log0() << "Building d2m";
    const size_t num_angles = abscissae.size();
    d2m_op = chi_math::Transpose(chi_math::Inverse(m2d_op));
    int i = 0;
    weights.clear();
    chi::log.Log0()<< "RESET WEIGHTS ";
    for(const auto& wt : d2m_op[0])
    {
      weights.push_back(wt);
      std::cout << std::setprecision(16) << weights[i] << ' ';
      ++i;
    }
    std::cout << std::endl;

    d2m_op_built = true;

  }
  if (method==2)
  {
    d2m_op.clear();
    std::vector<std::vector<double>> cmt;
    unsigned int num_angles = abscissae.size();
    unsigned int num_moms = m_to_ell_em_map.size();
    for (const auto &ell_em: m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);

      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        cur_mom.push_back(value);
      }
      cmt.push_back(cur_mom);
    }

//     solve for the weights
    std::vector<double> wt = {4.0*M_PI};
    for (size_t i = 1; i < weights.size(); ++i)
      wt.emplace_back(0.0);
    auto invt = chi_math::Inverse(cmt);
    chi::log.Log0() << "THE WT \n";
    chi_math::PrintVector(wt);
//    chi::log.Log0() << "The original cmt \n";
//    chi_math::PrintMatrix(cmt);
//    chi::log.Log0() << "The inverse cmt \n";
//    chi_math::PrintMatrix(invt);
    auto new_weights = chi_math::MatMul(invt, wt);
    chi::log.Log0() << "The old weights \n";
    chi_math::PrintVector(weights);
    chi::log.Log0() << "The new weights \n";
    chi_math::PrintVector(new_weights);
    chi::log.Log0() << "NEW SOLVE \n";
    chi_math::PrintVector(chi_math::GaussEliminationPivot(cmt,wt));
    weights = new_weights;

    for (const auto &ell_em: m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
//      double sumForMoments =0.0;
//      for (int m =0; m < num_angles;++m )
//      {
//        chi::log.Log0() << "THE L " << 0 << " THE M " << ell_em.m;
//        chi::log.Log0() << "DIRECTION " << m;
//        chi::log.Log0() << "OMEGA M YKM " << chi_math::Ylm(0,ell_em.m,
//                                                           abscissae[m].phi,abscissae[m].theta);
//        chi::log.Log0() << "OMEGA M YKm * Wm* " << weights[m]*chi_math::Ylm(0,ell_em.m,
//                                                           abscissae[m].phi,abscissae[m].theta);
//        sumForMoments += weights[m]*chi_math::Ylm(0,ell_em.m,abscissae[m].phi,abscissae[m].theta);
//      }
//      chi::log.Log0() << "THE YKM*Wm SUM " << sumForMoments;
      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        double w = weights[n];
        cur_mom.push_back(value * w);
      }
      d2m_op.push_back(cur_mom);
    }
    //Check the integral values of the harmonics
    for (const auto &ell_em: m_to_ell_em_map)
    {
      double integralVal = 0.0;
      double productIntegral = 0.0;
      chi::log.Log0() << "$$$$$$$$$$\nHARMONIC l=" << ell_em.ell << " m=" << ell_em.m;
      for (int j=0;j<num_angles;++j)
      {
        const auto &cur_angle = abscissae[j];
        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        integralVal += value*weights[j];
        productIntegral += value*value*weights[j];
}
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value " << integralVal;
//      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(integralVal))/
//                                              (4.0*M_PI/(2.0*ell_em.ell+1.0))<< "\n$$$$$";
      chi::log.Log0() << "Product Value " << productIntegral;
      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
//      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(productIntegral))/
//                                              (4.0*M_PI/(2.0*ell_em.ell+1.0));
    }

    d2m_op_built = true;
  }
////    Method 3 using gram-schmidtt orthogonalization
  if(method ==3)
  {
    d2m_op.clear();
    unsigned int num_angles = abscissae.size();
    unsigned int num_moms = m_to_ell_em_map.size();
    MatDbl cmt;
    //Make the coefficent matrix
    for (const auto &ell_em: m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);

      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
                                    cur_angle.phi,
                                    cur_angle.theta);
        cur_mom.push_back(value);
      }
      cmt.push_back(cur_mom);
    }

    // solve for the weights
    std::vector<double> wt = {4.0*M_PI};
    for (size_t i = 1; i < weights.size(); ++i)
      wt.emplace_back(0.0);
    chi::log.Log0() << "THE WT \n";
    chi_math::PrintVector(wt);
    auto new_weights = chi_math::GaussEliminationPivot(cmt,wt);
    chi::log.Log0() << "The old weights \n";
    for (auto& xi : weights)
      std::cout << std::setprecision(16) <<xi << ' ';
    std::cout << std::endl;
    chi::log.Log0() << "The new weights \n";
    for (auto& xi : new_weights)
      std::cout << std::setprecision(16) <<xi << ' ';
    std::cout << std::endl;
    weights = new_weights;


    //##############################################
    //CHECKING HARMONICS
    chi::log.Log0() << "CHECKING ORIGINAL HARMONICS";
    for (int i=0;i<m_to_ell_em_map.size();++i)
    {
      auto em_l_current = m_to_ell_em_map[i];
      double integralVal = 0.0;
      double productIntegral = 0.0;
      chi::log.Log0() << "$$$$$$$$$$\nHARMONIC l=" << em_l_current.ell
      << " m=" << em_l_current.m;
      for (int j=0;j<num_angles;++j)
      {
        double value = cmt[i][j];
//        const auto &cur_angle = abscissae[j];
//        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
//                                     cur_angle.phi,
//                                     cur_angle.theta);
        integralVal += value*weights[j];
        productIntegral += value*value*weights[j];
      }
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value         : " << std::setprecision(16) << integralVal;
      chi::log.Log0() << "Product Integral Value : " << std::setprecision(16) << productIntegral;
      chi::log.Log0() << "4pi/2l+1               : " << std::setprecision(16) << 4.0*M_PI/(2.0*em_l_current.ell+1.0);
      chi::log.Log0() << "Ratio of error         : " << std::setprecision(16) <<
      (4.0*M_PI/(2.0*em_l_current.ell+1.0)-abs(productIntegral))/
                                                        (4.0*M_PI/(2.0*em_l_current.ell+1.0));
    }

    //##############################################


    //Make the holder for the altered coefficients
    MatDbl cmt_hat=cmt;
    for (size_t i = 1;i<num_moms;++i)
    {
      VecDbl current_vec = cmt[i];
      VecDbl sum_val(num_angles);
      //Do GS- and go through every previous row to get the new row
      for (int j=0;j<i;++j)
      {
        //Do GS- and go through every previous row to get the new row
        VecDbl previous_vec_hat = cmt_hat[j];
        double multiplier = InnerProduct(previous_vec_hat, current_vec, weights)
                            / InnerProduct(previous_vec_hat,
                                           previous_vec_hat, weights);
        for (int o = 0; o < num_angles; ++o)
          sum_val[o] += previous_vec_hat[o]*multiplier;
      }
      cmt_hat[i] = current_vec - sum_val;
    }

    chi::log.Log() << "Printing the original cmt";
    chi_math::PrintMatrix(cmt);

    chi::log.Log() << "Printing cmt hat pre normalization matrix";
    chi_math::PrintMatrix(cmt_hat);


    //Now to normalize the values
    for (int i = 0; i<num_moms;++i)
    {
      auto ell = m_to_ell_em_map[i].ell;
      chi::log.Log0() << "Normalized to inner product of " << InnerProduct(cmt_hat[i],cmt_hat[i],weights);
      double normal = (4.0*M_PI)/(2.0*ell+1.0);
      double multiplier = sqrt(normal /
                               InnerProduct(cmt_hat[i],cmt_hat[i],weights));
      chi::log.Log0() << "Multiplier of " << multiplier;

      for (int k=0; k<num_angles;++k)
        cmt_hat[i][k] *= multiplier;
    }

    chi::log.Log() << "Printing cmt hat post normalization matrix";
    chi_math::PrintMatrix(cmt_hat);
  //##############################################
    //CHECKING HARMONICS
    chi::log.Log0() << "CHECKING APPROXIMATE HARMONICS";
    for (size_t i=0;i<m_to_ell_em_map.size();++i)
    {
      auto em_l_current = m_to_ell_em_map[i];
      double integralVal = 0.0;
      double productIntegral = 0.0;
      chi::log.Log0() << "$$$$$$$$$$\nHARMONIC l=" << em_l_current.ell
                      << " m=" << em_l_current.m;
      for (int j=0;j<num_angles;++j)
      {
        double value = cmt_hat[i][j];
//        const auto &cur_angle = abscissae[j];
//        double value = chi_math::Ylm(ell_em.ell, ell_em.m,
//                                     cur_angle.phi,
//                                     cur_angle.theta);
        integralVal += value*weights[j];
        productIntegral += value*value*weights[j];
      }
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value         : " << std::setprecision(16) << integralVal;
      chi::log.Log0() << "Product Integral Value : " << std::setprecision(16) << productIntegral;
      chi::log.Log0() << "4pi/2l+1               : " << std::setprecision(16) << 4.0*M_PI/(2.0*em_l_current.ell+1.0);
      chi::log.Log0() << "Ratio of error         : " << std::setprecision(16) <<
      (4.0*M_PI/(2.0*em_l_current.ell+1.0)-abs(productIntegral))/
                                                        (4.0*M_PI/(2.0*em_l_current.ell+1.0));
    }

    //##############################################


    chi::log.Log0() << "Check for orthogonalization";
    auto A = chi_math::MatMul(cmt_hat,chi_math::Transpose(cmt_hat));
    for (auto &row:A)
      for (auto &element:row)
        if (element<1e-10)
          element=0.0;
    chi_math::PrintMatrix(A);



    //Make the d2m matrix and m2d matrix
    MatDbl holder_m2d;
    for (size_t i = 0; i<num_moms;++i)
    {
      VecDbl temp_d2m;
      VecDbl temp_m2d;
      auto ell = m_to_ell_em_map[i].ell;
      for (int k=0; k<num_angles;++k)
      {
        temp_m2d.emplace_back(cmt_hat[i][k] * ((2.0*ell+1)/(4.0*M_PI)));
        temp_d2m.emplace_back(cmt_hat[i][k] * weights[k]);
      }
      d2m_op.push_back(temp_d2m);
      holder_m2d.push_back(temp_m2d);
    }
    //now we need to transpose the temporary m2d to get the actual m2d
//    m2d_op = chi_math::Transpose(holder_m2d); This is the normal method but due to a bug this must be entered without transpose
    m2d_op = holder_m2d;
    d2m_op_built = true;
  }
//  if (scattering_order < sn and method!=3 and method!=0)
//  {
//    FilterMoments(scattering_order);
//  }
}

void chi_math::AngularQuadratureTriangle::
BuildMomentToDiscreteOperator
  (unsigned int scattering_order,
   int dimension)
{
  if (m2d_op_built) return;
  MakeHarmonicIndices(scattering_order,dimension);
  if (method == 0)
  {
    const size_t num_angles = abscissae.size();
    const size_t num_moms = m_to_ell_em_map.size();

    const auto normalization =
      std::accumulate(weights.begin(), weights.end(), 0.0);

    for (const auto& ell_em : m_to_ell_em_map)
    {
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);

      for (int n=0; n<num_angles; n++)
      {
        const auto& cur_angle = abscissae[n];
        double value = ((2.0*ell_em.ell+1.0)/normalization)*
                       chi_math::Ylm(ell_em.ell,ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
        cur_mom.push_back(value);
      }

      m2d_op.push_back(cur_mom);
    }//for m
    m2d_op_built = true;
  }
  if (method == 1)
  {
    chi::log.Log0() << "Building m2d";
    const size_t num_angles = abscissae.size();
    const size_t num_moms = m_to_ell_em_map.size();

    const auto normalization = 4.0*M_PI;

    for (const auto& ell_em : m_to_ell_em_map)
    {
      double integral = 0.0;
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      chi::log.Log0() << "HARMONIC l=" << ell_em.ell << " m=" << ell_em.m;
      for (int n=0; n<num_angles; n++)
      {
        const auto& cur_angle = abscissae[n];
        double value = ((2.0*ell_em.ell+1.0)/normalization)*
                       chi_math::Ylm(ell_em.ell,ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
//        chi::log.Log0() << "ANGLE phi " <<  cur_angle.phi*180/M_PI << " theta " << cur_angle.theta*180/M_PI;
//        chi::log.Log0() << "Value " <<  value << " YLM " << chi_math::Ylm(ell_em.ell,ell_em.m,
//                                                                          cur_angle.phi,
//                                                                          cur_angle.theta);
//        integral += weights[n]*chi_math::Ylm(ell_em.ell,ell_em.m,
//                                  cur_angle.phi,
//                                  cur_angle.theta);

        cur_mom.push_back(value);
      }
//      chi::log.Log0() << "####################";
//      chi::log.Log0() << "Integral Value " << integral;
//      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
//      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(integral))/(4.0*M_PI/(2.0*ell_em.ell+1.0));

      m2d_op.push_back(cur_mom);
    }//for m
    m2d_op_built = true;
  }
  if (method == 2)
  {
    if (not d2m_op_built) BuildDiscreteToMomentOperator(sn,dimension);
    m2d_op = chi_math::Transpose(chi_math::Inverse(d2m_op));
    m2d_op_built = true;
  }
  if (method == 3)
  {
    if (not d2m_op_built)
      BuildDiscreteToMomentOperator(scattering_order,
                                    dimension);
    m2d_op_built = true;
  }
//  if (scattering_order < sn and method!=3 and method!=0)
//  {
//    FilterMoments(scattering_order);
//  }
  chi::log.Log0() << "The m2d op \n";
  chi_math::PrintMatrix(m2d_op);
  chi::log.Log0() << "Check Identity \n";
  auto B =chi_math::MatMul(m2d_op,chi_math::Transpose(d2m_op));
  for (auto &row:B)
    for (auto &element:row)
      if (element<1e-8)
        element=0.0;
  chi_math::PrintMatrix(B);
}