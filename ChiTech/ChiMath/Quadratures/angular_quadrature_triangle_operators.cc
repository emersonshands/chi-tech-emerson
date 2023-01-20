#include <algorithm>
#include <cmath>
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "angular_quadrature_triangle.h"
#include "ChiMath/chi_math.h"
#include <cassert>

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
  if (m2d_op_built and d2m_op_built)
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
  //We need to limit the scattering order to make sure this doesn't blow up
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

  chi::log.Log0() <<
  "$$$$$$$$$$$$$$$$$$$$$$$$$$\n Using the self made d2m\n$$$$$$$$$$$$$$$$$\n";
  chi::log.Log0() << "THIS IS SCATTERING ORDER " << scattering_order <<
  " AND DIMENSIONS " << dimension << "\n";
  if (d2m_op_built) return;
  MakeHarmonicIndices(scattering_order,dimension);
  if (method==1 or method==2)
  {
    d2m_op.clear();
    std::vector<std::vector<double>> cmt;
    unsigned int num_angles = abscissae.size();
    unsigned int num_moms = m_to_ell_em_map.size();
    // Changed this iteration scheme from the basic class
    // because the basic class would do the whole D2M matrix based
    // on the scattering order given and our method doesn't allow that
//    for (size_t mom_pos = 0; mom_pos<num_moms; ++mom_pos)
    for (const auto &ell_em: m_to_ell_em_map)
    {
//      const auto &ell_em =  m_to_ell_em_map[mom_pos];
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
      printf("This is the l %d and m %d \n",ell_em.ell,ell_em.m);
      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
                                    cur_angle.phi,
                                    cur_angle.theta);
        double w = weights[n];
        cur_mom.push_back(value*w);
      }
      if (method ==1)
      {
        d2m_op.push_back(cur_mom);
      }
      else
        cmt.push_back(cur_mom);
    }

    if (method == 2)
    {
      // solve for the weights
      //change this to change normalization of the weights
      double normalization = 4.0*M_PI;
      std::vector<double> wt = {normalization};
      for(size_t i = 1; i<weights.size(); ++i)
        wt.emplace_back(0.0);
      auto invt = chi_math::Inverse(cmt);
//      chi::log.Log0() << "THE WT \n";
//      chi_math::PrintVector(wt);
//      chi::log.Log0() << "The original cmt \n";
//      chi_math::PrintMatrix(cmt);
//      chi::log.Log0() << "The inverse cmt \n";
//      chi_math::PrintMatrix(invt);
      auto new_weights = chi_math::MatMul(invt,wt);
      weights = new_weights;
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
          double w = weights[n];
          cur_mom.push_back(value*w);
        }
        d2m_op.push_back(cur_mom);
      }
    }
    d2m_op_built = true;
  }
//    Method 3 using grah-schmidtt orthogonalization
  else if(not d2m_op_built and method ==3)
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
        double w = weights[n];
        cur_mom.push_back(value*w);
      }
      cmt.push_back(cur_mom);
    }
    //Make the holder for the altered coefficients
    MatDbl cmt_hat=cmt;
    size_t ndir = cmt[0].size();
    for (size_t i = 1;i<ndir;++i)
    {
      VecDbl current_vec = cmt[i];
      VecDbl sum_val(ndir);
      //Do GS- and go through every previous row to get the new row
      for (size_t j = 0; j<i; ++j)
      {
        VecDbl previous_vec_hat = cmt_hat[j];
        double multiplier = InnerProduct(previous_vec_hat,current_vec,weights)
                            / InnerProduct(previous_vec_hat,
                                           previous_vec_hat,weights);
        for (int o=0;o<ndir;++o)
          previous_vec_hat[o] *= multiplier;
        sum_val = previous_vec_hat+sum_val;
      }
      //Set the new row
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
    //Make the d2m matrix and m2d matrix
    MatDbl holder_m2d;
    for (int i = 0; i<ndir;++i)
    {
      VecDbl temp_d2m;
      VecDbl temp_m2d;
      for (int k=0; k<num_moms;++k)
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

//  if(moments!=0 and d2m_op_built and m2d_op_built and moments<scattering_order and method!=3)
//    FilterMoments();
}

void chi_math::AngularQuadratureTriangle::
BuildMomentToDiscreteOperator
  (unsigned int scattering_order,
   int dimension)
{
  if (not d2m_op_built) BuildDiscreteToMomentOperator(scattering_order,dimension);
  // Method 1 and 2 is just the inverse of d2m
  if (m2d_op_built) return;
  if (method == 1 or method == 2)
  {
    m2d_op = chi_math::Transpose(chi_math::Inverse(d2m_op));
  }
  //The M2D operator has to be built in the D2M build if method 3
  m2d_op_built = true;
  //Now filter by the moment number given
//  if(moments!=0 and d2m_op_built and m2d_op_built and moments<scattering_order)
//    FilterMoments();
  chi::log.Log0() << "The d2m op \n";
  chi_math::PrintMatrix(d2m_op);
  chi::log.Log0() << "The m2d op \n";
  chi_math::PrintMatrix(m2d_op);
  chi::log.Log0() << "The m2d*d2m op \n";
  chi_math::PrintMatrix(chi_math::MatMul(chi_math::Transpose(m2d_op),d2m_op));
  if (scattering_order < sn )
    FilterMoments(scattering_order);
  chi::log.Log0() << "AFTER FILTERING FOR SCATTERING ORDER \n";
  chi::log.Log0() << "The d2m op \n";
  chi_math::PrintMatrix(d2m_op);
  chi::log.Log0() << "The m2d op \n";
  chi_math::PrintMatrix(m2d_op);
}