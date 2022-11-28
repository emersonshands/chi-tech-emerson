#include "product_quadrature_op.h"
#include <cassert>
#include "chi_runtime.h"
#include "chi_log.h"

chi_math::ProductQuadratureOp::ProductQuadratureOp(const
chi_math::ProductQuadrature&
 inquad, const int inmethod, const int order, const int inmoments) :
 quad(inquad), method(inmethod), sn(order), moments(inmoments)
{
  //Check for the sizes of inquad to see if it matches for the method
  weights = quad.weights;
  omegas = quad.omegas;
  abscissae = quad.abscissae;
  azimu_ang = quad.azimu_ang;
  polar_ang = quad.polar_ang;
  map_directions = quad.GetDirectionMap();
  CheckInputs();
  chi_math::ProductQuadratureOp::MakeHarmonicIndices();
}

chi_math::ProductQuadratureOp::ProductQuadratureOp(const
                                                   chi_math::ProductQuadrature&
inquad, const int inmethod, const int order) :
  quad(inquad), method(inmethod), sn(order), moments(0)
{
  //Check for the sizes of inquad to see if it matches for the method
  weights = quad.weights;
  omegas = quad.omegas;
  abscissae = quad.abscissae;
  azimu_ang = quad.azimu_ang;
  polar_ang = quad.polar_ang;
  map_directions = quad.GetDirectionMap();
  CheckInputs();
  chi_math::ProductQuadratureOp::MakeHarmonicIndices();
}


void chi_math::ProductQuadratureOp::CheckInputs()
{
  //check the given inputs
  if ( azimu_ang.size()*polar_ang.size() != 2*sn*sn)
  {
    printf("Product of Azimuthal=%zu and Polar=%zu "
           "does not match 2*sn*sn=%i\n",azimu_ang.size(),
           polar_ang.size(),2*sn*sn);
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
  OptimizeForPolarSymmetry(4.0*M_PI);
}

void chi_math::ProductQuadratureOp::
OptimizeForPolarSymmetry(const double normalization)
{
  std::vector<chi_math::QuadraturePointPhiTheta> new_abscissae;
  std::vector<double>                            new_weights;
  std::vector<chi_mesh::Vector3>                 new_omegas;
  std::vector<double>                            new_polar_ang;
  std::vector<double>                            new_azimu_ang;

  const size_t num_pol = polar_ang.size();
  const size_t num_azi = azimu_ang.size();

  //Change this to azimuthal angle between 0 and pi
  // Grab all polar angles
  // we want only the azimuthal positive in y
  std::vector<unsigned int> new_azimu_map;
  new_polar_ang = polar_ang;
  for (size_t a = 0; a<num_azi; ++a)
  {
    if (azimu_ang[a]<M_PI)
    {
      new_azimu_ang.push_back(azimu_ang[a]);
      new_azimu_map.push_back(a);
    }
  }

  const size_t new_num_azimu = new_azimu_ang.size();
  double weight_sum = 0.0;
  for (size_t a=0; a<new_num_azimu; ++a)
    for (size_t p=0; p<num_pol; ++p)
    {
      const auto amap = new_azimu_map[a];
      const auto dmap = GetAngleNum(p,amap);
      new_weights.push_back(weights[dmap]);
      weight_sum += weights[dmap];
    }

  if (normalization > 0.0)
    for (double& w : new_weights)
      w *= normalization/weight_sum;

  chi_math::ProductQuadrature::AssembleCosines(new_azimu_ang,
                                               new_polar_ang,new_weights,false);
  polar_ang = new_polar_ang;
  azimu_ang = new_azimu_ang;
}


void chi_math::ProductQuadratureOp::FilterMoments()
{
  if (m2d_op_built and d2m_op_built)
  {
    int moments_to_keep = 1 + (moments*3 + moments*moments)/2;
    auto m2d_transposed = chi_math::Transpose(m2d_op);
    std::vector<std::vector<double>> m2dworking;
    std::vector<std::vector<double>> d2mworking;
    for (size_t i {}; i<moments_to_keep;++i)
    {
      d2mworking.push_back(d2m_op.at(i));
      m2dworking.push_back(m2d_transposed.at(i));
    }
    m2d_op = chi_math::Transpose(m2dworking);
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
  size_t fsize = f.size();
  double sum_val = 0.0;
  for (size_t i=0; i<fsize; ++i)
  {
    sum_val += f[i] * g[i] * wt[i];
  }
  return sum_val;
}


void chi_math::ProductQuadratureOp::MakeHarmonicIndices()
{
  const int nquad = static_cast<int>(sn);
  const int Lmax = 2 * (nquad - 1);
  for (int ell = 0; ell <= Lmax; ++ell)
    for (int m = 0; m <= ell; m += 1)
    {
      if (ell >= nquad)
      {
        if (m <= (ell - nquad)) continue;
        if (m >= nquad) continue;
        m_to_ell_em_map.emplace_back(ell, m);
      }
      else
      {
        m_to_ell_em_map.emplace_back(ell, m);
      }
    }
}

void chi_math::ProductQuadratureOp::BuildDiscreteToMomentOperator
()
{
  if (not d2m_op_built and (method==1 or method==2))
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
        double w = weights[n];
        cur_mom.push_back(value*w);
      }
      if (method ==1)
        d2m_op.push_back(cur_mom);
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
      auto new_weights = chi_math::MatMul(invt,wt);
      auto g = chi_math::GaussEliminationPivot(cmt,wt);
      new_weights = chi_math::GaussEliminationPivot(cmt,wt);
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
    //Method 3 using grah-schmidtt orthogonalization
  else if(not d2m_op_built and method ==3)
  {
    d2m_op.clear();
    unsigned int num_angles = abscissae.size();
    if (moments!=0 and moments!=sn)
      unsigned int num_moms = 1 + (moments*3 + moments*moments)/2;
    else
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
      for (int k=0; k<ndir;++k)
      {
        temp_m2d.emplace_back(cmt_hat[i][k] * ((2.0*i+1)/(4.0*M_PI)));
        temp_d2m.emplace_back(cmt_hat[i][k] * weights[k]);
      }
      d2m_op.push_back(temp_d2m);
      holder_m2d.push_back(temp_m2d);
    }
    //now we need to transpose the temporary m2d to get the actual m2d
    m2d_op = chi_math::Transpose(holder_m2d);
    d2m_op_built = true;
  }
  if(moments!=0 and d2m_op_built and m2d_op_built and moments!=sn and method!=3)
    FilterMoments();
}

void chi_math::ProductQuadratureOp::BuildMomentToDiscreteOperator
()
{
  if (not d2m_op_built) BuildDiscreteToMomentOperator();
  // Method 1 and 2 is just the inverse of d2m
  if (method == 1 or method == 2)
  {
    m2d_op = chi_math::Inverse(d2m_op);
  }
  //The M2D operator has to be built in the D2M build if method 3
  m2d_op_built = true;
  //Now filter by the moment number given
  if(moments!=0 and d2m_op_built and m2d_op_built and moments!=sn)
    FilterMoments();
}




