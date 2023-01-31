#include "LBSAdjointSolver/lbsadj_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"

//###################################################################
/**Computes the inner product of the flux and the material source.*/
double lbs_adjoint::AdjointSolver::ComputeInnerProduct()
{
  double local_integral = 0.0;

  auto pwl =
      std::dynamic_pointer_cast<chi_math::SpatialDiscretization_FE>(discretization);

  //============================================= Material sources
  for (const auto& cell : grid->local_cells)
  {
    if (matid_to_src_map.count(cell.material_id) == 0) continue; //Skip if no src

    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& source = matid_to_src_map[cell.material_id];
    const auto& fe_values = pwl->GetUnitIntegrals(cell);

    for (const auto& group : groups)
    {
      const int g = group.id;
      const double Q = source->source_value_g[g];

      if (Q > 0.0)
      {
        const int num_nodes = transport_view.NumNodes();
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t dof_map = transport_view.MapDOF(i, 0, g); //unknown map

          const double phi = phi_old_local[dof_map];

          local_integral += Q * phi * fe_values.IntV_shapeI(i);
        }//for node
      }//check source value >0
    }//for group
  }//for cell

  //============================================= Point sources
  for (const auto& point_source : point_sources)
  {
    const auto& info_list = point_source.ContainingCellsInfo();
    for (const auto& info : info_list)
    {
      const auto& cell = grid->local_cells[info.cell_local_id];
      const auto& transport_view = cell_transport_views[cell.local_id];
      const auto& source_strength = point_source.Strength();
      const auto& shape_values = info.shape_values;
      const auto& fe_values = pwl->GetUnitIntegrals(cell);

      for (const auto& group : groups)
      {
        const int g = group.id;
        const double S = source_strength[g] * info.volume_weight;

        if (S > 0.0)
        {
          const int num_nodes = transport_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            const size_t dof_map = transport_view.MapDOF(i, 0, g); //unknown map

            const double phi_i = phi_old_local[dof_map];

            local_integral += S * phi_i * shape_values[i];
          }//for node
        }//check source value >0
      }//for group
    }//for cell
  }//for point source

  double global_integral = 0.0;

  MPI_Allreduce(&local_integral,     //sendbuf
                &global_integral,    //recvbuf
                1, MPI_DOUBLE,       //count, datatype
                MPI_SUM,             //op
                MPI_COMM_WORLD);     //comm

  return global_integral;
}