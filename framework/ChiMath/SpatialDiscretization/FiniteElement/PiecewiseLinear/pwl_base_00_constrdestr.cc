#include "pwl_base.h"

//###################################################################
/**Constructor*/
chi_math::SpatialDiscretization_PWLBase::
  SpatialDiscretization_PWLBase(chi_mesh::MeshContinuumPtr &in_grid,
                                finite_element::SetupFlags in_setup_flags,
                                QuadratureOrder qorder,
                                SDMType in_type,
                                CoordinateSystemType in_cs_type
                                ) :
    SpatialDiscretization_FE(in_grid, in_cs_type, in_type,in_setup_flags),
    line_quad_order_arbitrary(qorder),
    tri_quad_order_arbitrary(qorder),
    quad_quad_order_arbitrary(qorder),
    tet_quad_order_arbitrary(qorder)
{}