#include "pwl_hexahedron.h"

//###################################################################
/**Computes cell volume and surface integrals.*/
void HexahedronPWLFEValues::
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data)
{
  size_t num_vol_qpoints = default_volume_quadrature.qpoints.size();
  size_t num_srf_qpoints = default_surface_quadrature.qpoints.size();

  //=========================================== Lambda for jacobian
  auto Jacobian = [this](double xi, double eta, double zeta)
  {
    chi_mesh::Matrix3x3 retJ;

    double dx_dxi   = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dx_dxi += _1_8th*xi_i[i]*(1.0+eta*eta_i[i])*(1.0+zeta*zeta_i[i])*x_i[i];
    double dx_deta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dx_deta += _1_8th*eta_i[i]*(1.0+xi*xi_i[i])*(1.0+zeta*zeta_i[i])*x_i[i];
    double dx_dzeta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dx_dzeta += _1_8th*zeta_i[i]*(1.0+xi*xi_i[i])*(1.0+eta*eta_i[i])*x_i[i];

    double dy_dxi   = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dy_dxi += _1_8th*xi_i[i]*(1.0+eta*eta_i[i])*(1.0+zeta*zeta_i[i])*y_i[i];
    double dy_deta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dy_deta += _1_8th*eta_i[i]*(1.0+xi*xi_i[i])*(1.0+zeta*zeta_i[i])*y_i[i];
    double dy_dzeta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dy_dzeta += _1_8th*zeta_i[i]*(1.0+xi*xi_i[i])*(1.0+eta*eta_i[i])*y_i[i];

    double dz_dxi   = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dz_dxi += _1_8th*xi_i[i]*(1.0+eta*eta_i[i])*(1.0+zeta*zeta_i[i])*z_i[i];
    double dz_deta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dz_deta += _1_8th*eta_i[i]*(1.0+xi*xi_i[i])*(1.0+zeta*zeta_i[i])*z_i[i];
    double dz_dzeta  = 0.0;
    for (unsigned i=0; i<num_nodes; ++i)
      dz_dzeta += _1_8th*zeta_i[i]*(1.0+xi*xi_i[i])*(1.0+eta*eta_i[i])*z_i[i];

    retJ.SetIJ(0,0,dx_dxi  ); retJ.SetIJ(0,0,dx_deta ); retJ.SetIJ(0,0,dx_dzeta);
    retJ.SetIJ(0,0,dy_dxi  ); retJ.SetIJ(0,0,dy_deta ); retJ.SetIJ(0,0,dy_dzeta);
    retJ.SetIJ(0,0,dz_dxi  ); retJ.SetIJ(0,0,dz_deta ); retJ.SetIJ(0,0,dz_dzeta);

    return retJ;
  };

  auto GetShape = [this](unsigned int i, const chi_mesh::Vector3& qpoint)
  {
    return HexShape(i,qpoint.x,qpoint.y,qpoint.z);
  };

  auto GetGradShape_xi_eta_zeta = [this](unsigned int i, const chi_mesh::Vector3& qpoint)
  {
    auto grad_xi_eta_zeta = chi_mesh::Vector3(
      HexGradShape_x(i,qpoint.x,qpoint.y,qpoint.z),
      HexGradShape_y(i,qpoint.x,qpoint.y,qpoint.z),
      HexGradShape_z(i,qpoint.x,qpoint.y,qpoint.z));

    return grad_xi_eta_zeta;
  };

  //============================================= Prefetch J and JTinv at QPs
  std::vector<chi_mesh::Matrix3x3> J_at_qp;
  std::vector<chi_mesh::Matrix3x3> JTinv_at_qp;
  std::vector<double>              detJ_at_qp;

  J_at_qp.reserve(num_vol_qpoints);
  JTinv_at_qp.reserve(num_vol_qpoints);
  detJ_at_qp.resize(num_vol_qpoints);

  for (const auto& qp : default_volume_quadrature.qpoints)
  {
    auto J     = Jacobian(qp.x,qp.y,qp.z);
    auto JT    = J.Transpose();
    auto JTinv = JT.Inverse();


    J_at_qp.push_back(J);
    JTinv_at_qp.push_back(JTinv);
    detJ_at_qp.push_back(J.Det());
  }

  //============================================= Prefetch shape, grad-shape and J
  std::vector<std::vector<double>>            shapes_i_at_qp;
  std::vector<std::vector<chi_mesh::Vector3>> grad_shapes_i_at_qp;

  shapes_i_at_qp.resize(num_nodes,std::vector<double>(num_vol_qpoints));
  grad_shapes_i_at_qp.resize(num_nodes,std::vector<chi_mesh::Vector3>(num_vol_qpoints));

  for (unsigned i=0; i<num_nodes; ++i)
    for (unsigned int q=0; q<num_vol_qpoints; ++q)
    {
      const auto qp = default_volume_quadrature.qpoints[q];

      shapes_i_at_qp[i][q] = GetShape(i,qp);

      auto grad_shape_xi_eta_zeta = GetGradShape_xi_eta_zeta(i,qp);
      const auto JTinv = JTinv_at_qp[q];

      grad_shapes_i_at_qp[i][q] = JTinv*grad_shape_xi_eta_zeta;
    }


  //============================================= Volume integrals
  ui_data.IntV_gradShapeI_gradShapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI_gradshapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI_shapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI.reserve(num_nodes);

  for (int i=0; i < num_nodes; i++)
  {
    std::vector<double>            gradijvalue_i(num_nodes, 0.0);
    std::vector<chi_mesh::Vector3> varphi_i_gradj(num_nodes, chi_mesh::Vector3());
    std::vector<double>            varphi_i_varphi_j(num_nodes, 0);

    // Computing
    // GradVarphi_i*GradVarphi_j and
    // Varphi_i*GradVarphi_j_xyz and
    // Varphi_i*Varphi_j
    for (int j=0; j < num_nodes; j++)
    {
      gradijvalue_i[j] = 0.0;
      varphi_i_varphi_j[j] = 0.0;

      for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
      {
        double weight = default_volume_quadrature.weights[qp];

        gradijvalue_i[j]
          += grad_shapes_i_at_qp[i][qp].Dot(grad_shapes_i_at_qp[j][qp]) *
             detJ_at_qp[qp] * weight;
      }//for qp

      for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
      {
        double weight = default_volume_quadrature.weights[qp];

        varphi_i_gradj[j] += shapes_i_at_qp[i][qp] *
                             grad_shapes_i_at_qp[j][qp] *
                             detJ_at_qp[qp] * weight;

        varphi_i_varphi_j[j] += shapes_i_at_qp[i][qp] *
                                shapes_i_at_qp[j][qp] *
                                detJ_at_qp[qp] * weight;
      }// for qp
    }// for j

    //Computing Varphi_i
    double  valuei_i = 0.0;
    for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
    {
      double weight = default_volume_quadrature.weights[qp];

      valuei_i += shapes_i_at_qp[i][qp] * detJ_at_qp[qp] * weight;
    }// for gp

    ui_data.IntV_gradShapeI_gradShapeJ.push_back(std::move(gradijvalue_i));
    ui_data.IntV_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradj));
    ui_data.IntV_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j));
    ui_data.IntV_shapeI.push_back(valuei_i);
  }
}