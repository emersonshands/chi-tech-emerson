#include "product_quadrature_op.h"
#include <cassert>
#include "chi_runtime.h"
#include "chi_log.h"
#include "algorithm"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

chi_math::ProductQuadratureOp::ProductQuadratureOp(const
chi_math::ProductQuadrature&
 inquad, const int inmethod, const int inorder, const int inmoments) :
 quad(inquad), method(inmethod), sn(inorder), moments(inmoments)
{
  //Check for the sizes of inquad to see if it matches for the method
  weights = quad.weights;
  omegas = quad.omegas;
  abscissae = quad.abscissae;
  azimu_ang = quad.azimu_ang;
  polar_ang = quad.polar_ang;
  map_directions = quad.GetDirectionMap();
  CheckInputs();
}

chi_math::ProductQuadratureOp::ProductQuadratureOp(const
                                                   chi_math::ProductQuadrature&
inquad, const int inmethod, const int inorder) :
  quad(inquad), method(inmethod), sn(inorder), moments(0)
{
  //Check for the sizes of inquad to see if it matches for the method
  weights = quad.weights;
  omegas = quad.omegas;
  abscissae = quad.abscissae;
  azimu_ang = quad.azimu_ang;
  polar_ang = quad.polar_ang;
  map_directions = quad.GetDirectionMap();
  CheckInputs();
}


void chi_math::ProductQuadratureOp::CheckInputs()
{
  //check the given inputs
  if ( azimu_ang.size()*polar_ang.size() != sn*sn)
  {
    printf("Product of Azimuthal=%zu and Polar=%zu "
           "does not match sn*sn=%i\n",azimu_ang.size(),
           polar_ang.size(),sn*sn);
    chi::log.Log0Error() << "Mismatch between polar*azimuthal and sn given";
    chi::Exit(510);
  }
  //check if the moments are more than possible
  if (moments<0 or moments>(sn))
  {
    printf("The moments user is asking for is greater than or less "
           "than the allowable moments for sn\n"
           "Given value %i, sn=%i\n",moments,sn);
    chi::log.Log0Error() << "Mismatch between given moments to cut and sn used";
    chi::Exit(510);
  }
//  OptimizeForPolarSymmetry(4.0*M_PI);
}

void chi_math::ProductQuadratureOp::FilterMoments(unsigned int scattering_order)
{
  if (m2d_op_built and d2m_op_built and moments!=0)
  {
    int moments_to_keep = 0;
    if (scattering_order >= sn)
    {
      int sn_adj = sn-1;
      moments_to_keep = 1 + (sn_adj*3 + sn_adj*sn_adj)/2;
      chi::log.Log0() << "MOMENTS TO KEEP " << moments_to_keep;
      for (int i =0;i<=scattering_order-sn;++i)
      {
        int next_layer = sn - i - 1;
        moments_to_keep = moments_to_keep + next_layer;
      }
    }
    else
      moments_to_keep = 1 + (scattering_order*3 +
                             scattering_order*scattering_order)/2;
    auto m2d_transposed = m2d_op;
    chi::log.Log0() << "Size of m2ell preresize " << m_to_ell_em_map.size();
    m_to_ell_em_map.resize(moments_to_keep);
    chi::log.Log0() << "Size of m2ell postresize " << m_to_ell_em_map.size();
    std::vector<std::vector<double>> m2dworking;
    std::vector<std::vector<double>> d2mworking;
    for (size_t i =0; i<moments_to_keep;++i)
    {
      chi::log.Log0() << "L " << m_to_ell_em_map[i].ell << " m "<< m_to_ell_em_map[i].m;
      d2mworking.push_back(d2m_op.at(i));
      m2dworking.push_back(m2d_transposed.at(i));
    }
    m2d_op = m2dworking;
    d2m_op = d2mworking;
  }
}

double chi_math::ProductQuadratureOp::InnerProduct(const VecDbl& f,
                                                   const VecDbl& g,
                                                   const VecDbl& wt)
{
  assert(!f.empty());
  assert(!g.empty());
  assert(!wt.empty());
  size_t fsize = wt.size();
  double sum_val = 0.0;
  for (size_t i=0; i<fsize; ++i)
  {
    sum_val += f[i] * g[i] * wt[i];
  }
  return sum_val;
}


void chi_math::ProductQuadratureOp::MakeHarmonicIndices
(unsigned int,  int dimension)
{
  if(not m_to_ell_em_map.empty()) return;
  printf("MAKING CUSTOM HARMONICS\n");
  const int nquad = sn;
  int Lmax = 2 * (nquad - 1);
  m_to_ell_em_map.clear();
  for (int ell = 0; ell <= Lmax+1; ++ell)
    for (int m = -ell; m <= ell; m += 1)
    {
      if (ell >= nquad)
      {
        if (abs(m)>ell - nquad and abs(m)<=nquad and m <nquad
            and (ell + abs(m)) % 2 == 0)
        {
          m_to_ell_em_map.emplace_back(ell, m);
          printf("Harmonics found l=%i m=%i \n", ell, m);
        } else
          continue;
      }
      else
        if ((ell+abs(m))%2==0)
      {
        m_to_ell_em_map.emplace_back(ell, m);
        printf("Harmonics found l=%i m=%i \n",ell,m);
      }
//      m_to_ell_em_map.emplace_back(ell, m);
//      printf("Harmonics found l=%i m=%i \n", ell, m);
    }
}

void chi_math::ProductQuadratureOp::BuildDiscreteToMomentOperator
(unsigned int scattering_order,int dimension)
{
//  chi::log.Log0() <<
//                  "$$$$$$$$$$$$$$$$$$$$$$$$$$\n Using the self made d2m\n$$$$$$$$$$$$$$$$$\n";
//  chi::log.Log0() << "THIS IS SCATTERING ORDER " << scattering_order <<
//                  " AND DIMENSIONS " << dimension << "\n";
  if (d2m_op_built) return;
  MakeHarmonicIndices(sn,dimension);
  if (method == 1)
  {
    if (not m2d_op_built) BuildMomentToDiscreteOperator(sn,dimension);
    chi::log.Log0() << "Building d2m";
    d2m_op = chi_math::Transpose(chi_math::Inverse(m2d_op));
////    std::ofstream DwriteFile("/home/grads/e/emersonshands01/CLASSwork/d2m.txt",std::ofstream::out);
//
//    for (int k=0;k<d2m_op[0].size();++k)
//    {
//      for (auto & j : d2m_op)
//        DwriteFile << std::setprecision(16) << j[k] << " ";
//      DwriteFile << "\n";
//    }
//    DwriteFile.close();

    int i =0;
    auto w_copy = weights;
    weights.clear();
//    chi::log.Log0() << "SIZE OF THE WEIGHTS " << weights.size();
//    std::cout.precision(18);
    for(const double& wt : d2m_op[0])
    {
      weights.push_back(wt);
//      std::cout << "\nCheck " << " " << weights[i] << " OLD VALUE " << w_copy[i];
      ++i;
    }

//    chi::log.Log0() << "SIZE OF THE WEIGHTS " << weights.size();
    d2m_op_built = true;
//    chi::Exit(99);
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
        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
                                    cur_angle.phi,
                                    cur_angle.theta);
        cur_mom.push_back(value);
      }
        cmt.push_back(cur_mom);
    }
    // solve for the weights
    //change this to change normalization of the weights
    double normalization = 4.0*M_PI;
    std::vector<double> wt = {normalization};
    for(size_t i = 1; i<weights.size(); ++i)
      wt.emplace_back(0.0);
    auto invt = chi_math::Inverse(cmt);
    auto new_weights = chi_math::MatMul(invt,wt);
//    chi::log.Log0() << "The old weights \n";
//    chi_math::PrintVector(weights);
//    chi::log.Log0() << "THE SUM OF WEIGHTS";
//    double sum = 0.0;
//    for (auto & k : weights)
//      sum += k;
//    chi::log.Log0() << sum;
//    chi::log.Log0() << "The new weights \n";
//    chi_math::PrintVector(new_weights);
    weights.clear();
    weights = new_weights;
//    sum = 0.0;
//    for (auto & k : weights)
//      sum += k;
//    chi::log.Log0() << sum;
//    chi::Exit(99);
    for (const auto &ell_em: m_to_ell_em_map)
    {
      std::vector<double> cur_moment_new;
      cur_moment_new.reserve(num_angles);

      for (int n = 0; n < num_angles; n++)
      {
        const auto &cur_angle = abscissae[n];
        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
                                    cur_angle.phi,
                                    cur_angle.theta);
        double w = weights[n];
        cur_moment_new.push_back(value*w);
      }
      d2m_op.push_back(cur_moment_new);
    }
    d2m_op_built = true;
  }
//    Method 3 using GS orthogonalization
  if(method ==3)
  {
    std::vector<std::vector<double>> cmt;
    unsigned int num_angles = abscissae.size();
    unsigned int num_moms = m_to_ell_em_map.size();
//    for (const auto &ell_em: m_to_ell_em_map)
//    {
//      std::vector<double> cur_mom;
//      cur_mom.reserve(num_angles);
//
//      for (int n = 0; n < num_angles; n++)
//      {
//        const auto &cur_angle = abscissae[n];
//        double value =chi_math::Ylm(ell_em.ell, ell_em.m,
//                                    cur_angle.phi,
//                                    cur_angle.theta);
//        double w = weights[n];
//        cur_mom.push_back(value*w);
//      }
//      cmt.push_back(cur_mom);
//    }
//    chi_math::PrintMatrix(cmt);
//    chi::log.Log0() << "WEIGHTS GIVEN$$$$$$$$$$$$$$$$$$$";
//    chi_math::PrintVector(weights);
//    chi_math::condition(cmt);
    // solve for the weights
    //change this to change normalization of the weights
//    double normalization = 4.0*M_PI;
//    std::vector<double> wt = {normalization};
//    for(size_t i = 1; i<weights.size(); ++i)
//      wt.emplace_back(0.0);
//    auto invt = chi_math::Inverse(cmt);
//    auto new_weights = chi_math::MatMul(invt,wt);
//    weights.clear();
//    weights = new_weights;
//
//    d2m_op.clear();
//    printf("the weights\n");
//    chi::log.Log0() << "WEIGHTS GIVEN$$$$$$$$$$$$$$$$$$$";
//    chi_math::PrintVector(weights);
//    chi_math::condition(cmt);
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
//        printf("Weights %0.4f and points %0.4f \n",w,value);
        cur_mom.push_back(value);
      }
      cmt.push_back(cur_mom);
    }

    //Make the holder for the altered coefficients
    MatDbl cmt_hat=cmt;
    size_t ndir = cmt[0].size();
    size_t nmom = cmt.size();
    for (size_t i = 1;i<nmom;++i)
    {
      VecDbl current_vec = cmt[i];
      VecDbl sum_val(ndir);
//      chi::log.Log0() << "SUM VECTOR \n";
//      chi_math::PrintVector(sum_val);
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
//    chi::log.Log0() << "$$$$$$$$$$$$$$$$$$$$$$$$CMT MATRIX";
//    chi_math::PrintMatrix(cmt);
//    chi::log.Log0() << "WEIGHTS GIVEN$$$$$$$$$$$$$$$$$$$";
//    chi_math::PrintVector(weights);
//    chi_math::condition(cmt);
//    chi::log.Log0() << "$$$$$$$$$$$$ CMT_HAT MATRIX";
//    chi_math::PrintMatrix(cmt_hat);
    //Now to normalize the values
    for (int i = 0; i<nmom;++i)
    {
      double normal = (4.0*M_PI)/(2.0*i+1.0);
      double multiplier = sqrt(normal /
                               InnerProduct(cmt_hat[i],cmt_hat[i],weights));
      chi::log.Log0() << "MULTIPLIER " << multiplier;
      for (int k=0; k<ndir;++k)
        cmt_hat[i][k] *= multiplier;
    }
//    chi::log.Log0() << "$$$$$$$$$$$$ CMT_HAT MATRIX";
//    chi_math::PrintMatrix(cmt_hat);
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
    m2d_op = holder_m2d;
    d2m_op_built = true;
//    chi::log.Log0() << "############D2M#####################";
//    chi_math::PrintMatrix(d2m_op);
//    chi::log.Log0() << "############m2d#####################";
//    chi_math::PrintMatrix(m2d_op);

  }
//  chi_math::condition(d2m_op);
  if (scattering_order < 2 * (sn - 1))
  {
    FilterMoments(scattering_order);
  }
}

void chi_math::ProductQuadratureOp::BuildMomentToDiscreteOperator
(unsigned int scattering_order,int dimension)
{
  std::cout.precision(16);
  if (m2d_op_built) return;
  MakeHarmonicIndices(sn,dimension);
  if (method ==1)
  {
    chi::log.Log0() << "Building m2d";
    const size_t num_angles = abscissae.size();
    const size_t num_moms = m_to_ell_em_map.size();

    const double normalization = 4.0*M_PI;
//    std::ofstream writeFile("/home/grads/e/emersonshands01/CLASSwork/M2D.txt",std::ofstream::out);

    for (const auto& ell_em : m_to_ell_em_map)
    {
//      double integral = 0.0;
      std::vector<double> cur_mom;
      cur_mom.reserve(num_angles);
//      chi::log.Log0() << "HARMONIC l=" << ell_em.ell << " m=" << ell_em.m;
      for (int n=0; n<num_angles; n++)
      {
        const auto& cur_angle = abscissae[n];
        double value = ((2.0*ell_em.ell+1.0)/normalization)*
                       chi_math::Ylm(ell_em.ell,ell_em.m,
                                     cur_angle.phi,
                                     cur_angle.theta);
//        integral += weights[n]*(chi_math::Ylm(ell_em.ell,ell_em.m,
//                                  cur_angle.phi,
//                                  cur_angle.theta)*
//                                    chi_math::Ylm(ell_em.ell,ell_em.m,
//                                                                 cur_angle.phi,
//                                                                 cur_angle.theta));
        cur_mom.push_back(value);
//        writeFile << std::setprecision(16) << value << " ";
      }
//      chi::log.Log0() << "####################";
//      chi::log.Log0() << "Integral Value " << integral;
//      chi::log.Log0() << "4pi/2l+1 " << 4.0*M_PI/(2.0*ell_em.ell+1.0);
//      chi::log.Log0() << "Ratio of error " << (4.0*M_PI/(2.0*ell_em.ell+1.0)-integral)/4.0*M_PI/(2.0*ell_em.ell+1.0);
//      writeFile << "\n";
      m2d_op.push_back(cur_mom);
    }//for m
    m2d_op_built = true;
//    writeFile.close();
//    chi::Exit(99);
  }
  if (method == 2)
  {
    if (not d2m_op_built)
      BuildDiscreteToMomentOperator(scattering_order,
                                    dimension);
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
  if (scattering_order < 2 * (sn - 1))
  {
    FilterMoments(scattering_order);
  }
//  chi::log.Log0() << "The m2d op \n";
//  chi_math::PrintMatrix(m2d_op);
}