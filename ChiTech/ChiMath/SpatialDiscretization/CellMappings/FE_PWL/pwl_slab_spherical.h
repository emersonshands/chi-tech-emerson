#ifndef SLAB_MAPPING_FE_PWL_SPHERICAL_H
#define SLAB_MAPPING_FE_PWL_SPHERICAL_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_slab.h"
#include "pwl_slab.h"


namespace chi_math
{
  /** Object for handling slab shaped 1D cells in point-symmetric
     *  spherical coordinates. */
  class SlabMappingFE_PWL_Spherical : public chi_math::SlabMappingFE_PWL
  {
  //  Methods
  public:
    SlabMappingFE_PWL_Spherical(
      const chi_mesh::Cell& slab_cell,
      const chi_mesh::MeshContinuumConstPtr& ref_grid,
      const QuadratureLine& volume_quadrature)
    : SlabMappingFE_PWL(slab_cell, ref_grid, volume_quadrature)
    {}
  private:
    double SpatialWeightFunction(const chi_mesh::Vector3& pt) const override
    { return pt[2]*pt[2]; }
  public:
    void
    ComputeUnitIntegrals(finite_element::UnitIntegralData& ui_data) const override
    { ComputeWeightedUnitIntegrals(ui_data); }
  };
}

#endif // SLAB_MAPPING_FE_PWL_SPHERICAL_H
