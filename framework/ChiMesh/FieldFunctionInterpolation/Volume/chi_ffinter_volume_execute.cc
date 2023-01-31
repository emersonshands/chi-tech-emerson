#include "chi_ffinter_volume.h"

#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Executes the volume interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Execute()
{
  const auto& ref_ff = *field_functions.front();
  const auto& sdm    = ref_ff.SDM();
  const auto& grid   = *sdm.ref_grid;

  const auto& uk_man = ref_ff.UnkManager();
  const auto uid = 0;
  const auto cid = m_ref_component;

  using namespace chi_mesh::ff_interpolation;
  const auto field_data = ref_ff.GetGhostedFieldVector();

  double local_volume = 0.0;
  double local_sum = 0.0;
  double local_max = 0.0;
  double local_min = 0.0;
  for (const uint64_t cell_local_id : cell_local_ids_inside_logvol)
  {
    const auto& cell = grid.local_cells[cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto qp_data = cell_mapping.MakeVolumeQuadraturePointData();

    std::vector<double> node_dof_values(num_nodes, 0.0);
    for (size_t i=0; i<num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell,i,uk_man,uid,cid);
      node_dof_values[i] = field_data[imap];
    }//for i

    if (cell_local_id == cell_local_ids_inside_logvol.front())
    {
      local_max = node_dof_values.front();
      local_min = node_dof_values.front();
    }

    for (size_t i=0; i<num_nodes; ++i)
    {
      local_max = std::fmax(node_dof_values[i], local_max);
      local_min = std::fmin(node_dof_values[i], local_min);
    }

    for (const size_t qp : qp_data.QuadraturePointIndices())
    {
      double ff_value = 0.0;
      for (size_t j=0; j<num_nodes; ++j)
        ff_value += qp_data.ShapeValue(j,qp) * node_dof_values[j];

      double function_value = ff_value;
      if (op_type >= Operation::OP_SUM_LUA and
          op_type <= Operation::OP_MAX_LUA)
        function_value = CallLuaFunction(ff_value, cell.material_id);

      local_volume += qp_data.JxW(qp);
      local_sum += function_value * qp_data.JxW(qp);
      local_max = std::fmax(ff_value, local_max);
      local_min = std::fmin(ff_value, local_min);
    }//for qp
  }//for cell-id

  if (op_type == Operation::OP_SUM or op_type == Operation::OP_SUM_LUA)
  {
    double global_sum;
    MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    op_value = global_sum;
  }
  if (op_type == Operation::OP_AVG or op_type == Operation::OP_AVG_LUA)
  {
    double local_data[] = {local_volume, local_sum};
    double global_data[] = {0.0,0.0};

    MPI_Allreduce(&local_data,&global_data,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    double global_volume = global_data[0];
    double global_sum = global_data[1];
    op_value = global_sum/global_volume;
  }
  if (op_type == Operation::OP_MAX or op_type == Operation::OP_MAX_LUA)
  {
    double global_value;
    MPI_Allreduce(&local_max,&global_value,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    op_value = global_value;
  }
}