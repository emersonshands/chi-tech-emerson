#include <algorithm>
#include <cmath>
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "angular_quadrature_triangle.h"
#include "ChiMath/chi_math.h"
#include <cassert>
#include <numeric>

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
    int moments_to_keep = 1 +
      (scattering_order*3 + scattering_order*scattering_order)/2;
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
MakeHarmonicIndices(unsigned int, int dimension)
{
  int L = static_cast<int>(sn);
  if (m_to_ell_em_map.empty())
  {
    for (int ell = 0; ell <= L; ++ell)
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
  MakeHarmonicIndices(sn,dimension);
  if (method == 1)
  {
    d2m_op.clear();
    if (not m2d_op_built) BuildMomentToDiscreteOperator(sn, dimension);
    chi::log.Log0() << "Building d2m";
    const size_t num_angles = abscissae.size();
    d2m_op = chi_math::Transpose(chi_math::Inverse(m2d_op));
    int i = 0;
    weights.clear();
    for(const auto& wt : d2m_op[0])
    {
      weights.push_back(wt);
      chi::log.Log0()<< "RESET WEIGHTS " << weights[i];
      ++i;
    }

//    for (int ell_em = 0; ell_em <m_to_ell_em_map.size(); ++ell_em)
//    {
//      chi::log.Log0() << "L " << m_to_ell_em_map[ell_em].ell << " M " << m_to_ell_em_map[ell_em].m;
//
//      for (int n = 0; n < num_angles; n++)
//      {
//        const auto &cur_angle = abscissae[n];
//        double value = d2m_op[ell_em][n];
////        chi::log.Log0() << "ANGLE phi " << cur_angle.phi * 180 / M_PI << " theta " << cur_angle.theta * 180 / M_PI;
////        chi::log.Log0() << "Value " << value;
//      }

//    }
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

    // solve for the weights
//    std::vector<double> wt = {4.0*M_PI};
//    for (size_t i = 1; i < weights.size(); ++i)
//      wt.emplace_back(0.0);
//    auto invt = chi_math::Inverse(cmt);
//    chi::log.Log0() << "THE WT \n";
//    chi_math::PrintVector(wt);
////    chi::log.Log0() << "The original cmt \n";
////    chi_math::PrintMatrix(cmt);
////    chi::log.Log0() << "The inverse cmt \n";
////    chi_math::PrintMatrix(invt);
//    auto new_weights = chi_math::MatMul(invt, wt);
//    chi::log.Log0() << "The old weights \n";
//    chi_math::PrintVector(weights);
//    chi::log.Log0() << "The new weights \n";
//    chi_math::PrintVector(new_weights);
//    chi::log.Log0() << "NEW SOLVE \n";
//    chi_math::PrintVector(chi_math::GaussEliminationPivot(cmt,wt));
//    weights = new_weights;
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

    d2m_op_built = true;
  }
////    Method 3 using grah-schmidtt orthogonalization
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

    //Make the holder for the altered coefficients
    MatDbl cmt_hat=cmt;
    size_t ndir = cmt[0].size();
    size_t nmom = cmt.size();
    for (size_t i = 1;i<nmom;++i)
    {
      VecDbl current_vec = cmt[i];
      VecDbl sum_val(ndir);
      //Do GS- and go through every previous row to get the new row
      for (int j=0;j<i;++j)
      {
        //Do GS- and go through every previous row to get the new row
        VecDbl previous_vec_hat = cmt_hat[j];
        double multiplier = InnerProduct(previous_vec_hat, current_vec, weights)
                            / InnerProduct(previous_vec_hat,
                                           previous_vec_hat, weights);
//        chi::log.Log0() << "MULTIPLIER " << multiplier;
        for (int o = 0; o < ndir; ++o)
          previous_vec_hat[o] *= multiplier;
        auto temp = sum_val + previous_vec_hat;
        sum_val = temp;
//        chi::log.Log0() << "New vector \n";
//        chi_math::PrintVector(sum_val);
        //Set the new row
      }
      cmt_hat[i] = cmt[i] - sum_val;
    }

    //Now to normalize the values
    for (int i = 0; i<ndir;++i)
    {
      double normal = (4.0*M_PI)/(2.0*i+1.0);
      double multiplier = sqrt(normal /
                               InnerProduct(cmt_hat[i],cmt_hat[i],weights));
      for (int k=0; k<ndir;++k)
        cmt_hat[i][k] *= multiplier;
    }
    //Check the integral values of the harmonics
    for (const auto &ell_em: m_to_ell_em_map)
    {
      double integralVal = 0.0;
      chi::log.Log0() << "$$$$$$$$$$\nHARMONIC l=" << ell_em.ell << " m=" << ell_em.m;
      for (int j=0;j<ndir;++j)
      {
        integralVal += cmt_hat[ell_em.ell][j]*weights[j];
      }
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value " << integralVal;
      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(integralVal))/
      (4.0*M_PI/(2.0*ell_em.ell+1.0));
    }
    //Make the d2m matrix and m2d matrix
    MatDbl holder_m2d;
    for (int i = 0; i<nmom;++i)
    {
      VecDbl temp_d2m;
      VecDbl temp_m2d;
      for (int k=0; k<ndir;++k)
      {
        temp_m2d.emplace_back(cmt_hat[i][k] * ((2.0*i+1)/(4.0*M_PI)));
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
//  chi::log.Log0() << "The d2m op \n";
//  chi_math::PrintMatrix(d2m_op);
//  if (m2d_op_built)
//  {
//    chi::log.Log0() << "The m2d*d2m op \n";
//    auto check = chi_math::MatMul(chi_math::Transpose(m2d_op), d2m_op);
//    for (int i = 0; i < check.size(); ++i)
//      for (int j = 0; j < check[0].size(); ++j)
//        if (check[i][j] < 1e-8)
//          check[i][j] = 0.0;
//    chi_math::PrintMatrix(check);
//  }
  if (scattering_order < sn)
  {
    FilterMoments(scattering_order);
  }
//  chi::Exit(0);
}

void chi_math::AngularQuadratureTriangle::
BuildMomentToDiscreteOperator
  (unsigned int scattering_order,
   int dimension)
{
//  if (not d2m_op_built) BuildDiscreteToMomentOperator(sn,dimension);
  // Method 1 and 2 is just the inverse of d2m
  if (m2d_op_built) return;
  MakeHarmonicIndices(sn,dimension);
  if (method ==1)
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
        integral += weights[n]*chi_math::Ylm(ell_em.ell,ell_em.m,
                                  cur_angle.phi,
                                  cur_angle.theta);

        cur_mom.push_back(value);
      }
      chi::log.Log0() << "####################";
      chi::log.Log0() << "Integral Value " << integral;
      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-abs(integral))/(4.0*M_PI/(2.0*ell_em.ell+1.0));

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
//  chi::log.Log0() << "The m2d op \n";
//  chi_math::PrintMatrix(m2d_op);
//  if (d2m_op_built)
//  {
//    chi::log.Log0() << "The m2d*d2m op \n";
//    auto check = chi_math::MatMul(chi_math::Transpose(m2d_op), d2m_op);
//    for (int i = 0; i < check.size(); ++i)
//      for (int j = 0; j < check[0].size(); ++j)
//        if (check[i][j] < 1e-8)
//          check[i][j] = 0.0;
//    chi_math::PrintMatrix(check);
//  }
  if (scattering_order < sn)
  {
    FilterMoments(scattering_order);
  }
  chi::log.Log0() << "The m2d op \n";
  chi_math::PrintMatrix(m2d_op);

}