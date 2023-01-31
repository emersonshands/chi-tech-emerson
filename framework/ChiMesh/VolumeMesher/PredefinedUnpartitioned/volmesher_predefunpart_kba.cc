#include "volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"


#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

;


//###################################################################
/** Applies KBA-style partitioning to the mesh.*/
std::vector<int64_t> chi_mesh::VolumeMesherPredefinedUnpartitioned::
  KBA(const chi_mesh::UnpartitionedMesh& umesh)
{
  chi::log.Log() << "Partitioning mesh KBA-style.";

  const size_t num_raw_cells = umesh.raw_cells.size();

  //======================================== Lambda to get partition-id
  //                                         from centroid
  auto GetPIDFromCentroid = [](const chi_mesh::Vertex& centroid)
  {
    auto& handler = chi_mesh::GetCurrentHandler();

    int Px = handler.volume_mesher->options.partition_x;
    int Py = handler.volume_mesher->options.partition_y;

    chi_mesh::Cell temp_cell(CellType::GHOST, CellType::GHOST);
    temp_cell.centroid = centroid;

    auto xyz = GetCellXYZPartitionID(&temp_cell);

    int nxi = std::get<0>(xyz);
    int nyi = std::get<1>(xyz);
    int nzi = std::get<2>(xyz);

    return nzi*Px*Py + nyi*Px + nxi;
  };

  //======================================== Determine cell partition-IDs
  //                                         only on home location
  std::vector<int64_t> cell_pids(num_raw_cells, 0);
  if (chi::mpi.location_id == 0)
  {
    uint64_t cell_id = 0;
    for (auto& raw_cell : umesh.raw_cells)
      cell_pids[cell_id++] = GetPIDFromCentroid(raw_cell->centroid);
  }//if home location

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),                 //buffer [IN/OUT]
            static_cast<int>(num_raw_cells),  //count
            MPI_LONG_LONG_INT,                //data type
            0,                                //root
            MPI_COMM_WORLD);                  //communicator
  chi::log.Log() << "Done partitioning mesh.";

  return cell_pids;
}