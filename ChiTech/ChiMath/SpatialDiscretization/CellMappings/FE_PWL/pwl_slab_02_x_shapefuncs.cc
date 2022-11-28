#include "pwl_slab.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

/**Shape function i evaluated at given point for the slab.*/
double chi_math::SlabMappingFE_PWL::
  ShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  const auto& p0 = m_grid_ptr->vertices[v0i];
  const auto& p1 = m_grid_ptr->vertices[v1i];
  chi_mesh::Vector3 xyz_ref = xyz - p0;

  chi_mesh::Vector3 v01 = p1 - p0;

  double xi   = v01.Dot(xyz_ref)/v01.Norm()/h;

  if ((xi>=-1.0e-6) and (xi<=1.0+1.0e-6))
  {
    if (i==0)
      return 1.0 - xi;
    else
      return xi;
  }//if in cell

  return 0.0;
}

//#################################################################
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void chi_math::SlabMappingFE_PWL::
  ShapeValues(const chi_mesh::Vector3& xyz,
              std::vector<double>& shape_values) const
{
  shape_values.resize(m_num_nodes, 0.0);
  const auto& p0 = m_grid_ptr->vertices[v0i];
  const auto& p1 = m_grid_ptr->vertices[v1i];
  chi_mesh::Vector3 xyz_ref = xyz - p0;

  chi_mesh::Vector3 v01 = p1 - p0;

  double xi   = v01.Dot(xyz_ref)/v01.Norm()/h;

  if ((xi>=-1.0e-6) and (xi<=1.0+1.0e-6))
  {
    for (int i=0; i < m_num_nodes; i++)
    {
      if (i==0)
        shape_values[i] = 1.0 - xi;
      else
        shape_values[i] = xi;
    }//for dof

    return;
  }//if in cell

}

//###################################################################
/**Returns the evaluation of grad-shape function i at the supplied point.*/
chi_mesh::Vector3 chi_math::SlabMappingFE_PWL::
  GradShapeValue(const int i, const chi_mesh::Vector3& xyz) const
{
  if (i==0)
    return chi_mesh::Vector3(0.0, 0.0, -1.0 / h);
  else
    return chi_mesh::Vector3(0.0, 0.0, 1.0 / h);
}

//###################################################################
/**Populates shape_values with the value of each shape function's
 * value evaluate at the supplied point.*/
void chi_math::SlabMappingFE_PWL::
  GradShapeValues(const chi_mesh::Vector3& xyz,
                  std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  gradshape_values.clear();
  gradshape_values.emplace_back(GradShapeValue(0,xyz));
  gradshape_values.emplace_back(GradShapeValue(1,xyz));
}