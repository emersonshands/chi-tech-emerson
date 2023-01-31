#include "pwl_slab.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

void chi_math::SlabMappingFE_PWL::InitializeVolumeQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data) const
{
  //=================================== Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = volume_quadrature.qpoints.size();

  //=================================== Declare necessary vars
  std::vector<unsigned int>     V_quadrature_point_indices;
  VecVec3                       V_qpoints_xyz;
  std::vector<VecDbl>           V_shape_value;
  std::vector<VecVec3>          V_shape_grad;
  VecDbl                        V_JxW;
  size_t                        V_num_nodes;

  //=================================== Init volumetric quadrature
  V_quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    V_quadrature_point_indices.push_back(qp);

  V_shape_value.reserve(m_num_nodes);
  V_shape_grad.reserve(m_num_nodes);
  for (size_t i=0; i < m_num_nodes; i++)
  {
    VecDbl  node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (const auto& qpoint : volume_quadrature.qpoints)
    {
      node_shape_value.push_back(SlabShape(i,qpoint));
      node_shape_grad.emplace_back(0.0,                //x
                                   0.0,                //y
                                   SlabGradShape(i));  //z
    }//for qp

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  }//for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  const double J = h;
  for (size_t qp=0; qp<ttl_num_vol_qpoints; ++qp)
  {
    const double w = volume_quadrature.weights[qp];
    V_JxW.push_back(J * w);

    const double qp_xyz_tilde = volume_quadrature.qpoints[qp][0];
    V_qpoints_xyz.push_back(v0 + J * chi_mesh::Vector3(0.0,0.0,qp_xyz_tilde));
  }//for qp

  V_num_nodes = m_num_nodes;

  internal_data.InitializeData(V_quadrature_point_indices,
                               V_qpoints_xyz,
                               V_shape_value,
                               V_shape_grad,
                               V_JxW,
                               face_node_mappings,
                               V_num_nodes);
}

void chi_math::SlabMappingFE_PWL::InitializeFaceQuadraturePointData(unsigned int face,
                                                                    chi_math::finite_element::FaceQuadraturePointData& faces_qp_data) const
{
  const bool ON_SURFACE = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = 1;

  unsigned int f=face;
  {
    //=================================== Declare necessary vars
    std::vector<unsigned int>     F_quadrature_point_indices;
    VecVec3                       F_qpoints_xyz;
    std::vector<VecDbl>           F_shape_value;
    std::vector<VecVec3>          F_shape_grad;
    VecDbl                        F_JxW;
    VecVec3                       F_normals;
    size_t                        F_num_nodes;

    size_t ttl_num_face_qpoints = num_srf_qpoints;

    F_quadrature_point_indices.reserve(ttl_num_face_qpoints);
    for (unsigned int qp=0; qp<ttl_num_face_qpoints; ++qp)
      F_quadrature_point_indices.push_back(qp);

    F_normals.reserve(ttl_num_face_qpoints);
    for (size_t qp=0; qp<ttl_num_face_qpoints; ++qp)
      F_normals.push_back(normals[f]);

    F_shape_value.reserve(m_num_nodes);
    F_shape_grad.reserve(m_num_nodes);
    for (size_t i=0; i < m_num_nodes; i++)
    {
      VecDbl  node_shape_value;
      VecVec3 node_shape_grad;

      node_shape_value.reserve(ttl_num_face_qpoints);
      node_shape_grad.reserve(ttl_num_face_qpoints);

      for (const auto& qpoint : {chi_mesh::Vector3(0.0,0.0,0.0)})
      {
        node_shape_value.push_back(SlabShape(i,qpoint,ON_SURFACE,f));
        node_shape_grad.emplace_back(0.0,                //x
                                     0.0,                //y
                                     SlabGradShape(i));  //z
      }//for qp
      F_shape_value.push_back(node_shape_value);
      F_shape_grad.push_back(node_shape_grad);
    }//for i

    F_JxW.reserve(ttl_num_face_qpoints);
    F_qpoints_xyz.reserve(ttl_num_face_qpoints);
    const double JxW = 1.0;
    for (size_t qp=0; qp<num_srf_qpoints; ++qp)
    {
      F_JxW.push_back(JxW);

      F_qpoints_xyz.push_back(chi_mesh::Vector3(0.0,0.0,f));
    }

    F_num_nodes = 1;

    faces_qp_data.InitializeData(F_quadrature_point_indices,
                                 F_qpoints_xyz,
                                 F_shape_value,
                                 F_shape_grad,
                                 F_JxW,
                                 F_normals,
                                 face_node_mappings,
                                 F_num_nodes);
  }//face
}
