#include "volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include "ChiTimer/chi_timer.h"
#include "ChiConsole/chi_console.h"

//###################################################################
/**Executes the predefined3D mesher.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::Execute()
{
  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " VolumeMesherPredefinedUnpartitioned executing. Memory in use = "
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB"
    << std::endl;

  //======================================== Check partitioning params
  if (options.partition_type == KBA_STYLE_XYZ)
  {
    int Px = this->options.partition_x;
    int Py = this->options.partition_y;
    int Pz = this->options.partition_z;

    int desired_process_count = Px*Py*Pz;

    if (desired_process_count != chi::mpi.process_count)
    {
      chi::log.LogAllError()
        << "ERROR: Number of processors available ("
        << chi::mpi.process_count <<
        ") does not match amount of processors "
        "required by partitioning parameters ("
        << desired_process_count << ").";
     chi::Exit(EXIT_FAILURE);
    }
  }

  //======================================== Get unpartitioned mesh
  auto umesh = m_umesh;

  chi::log.Log() << "Computed centroids";
  MPI_Barrier(MPI_COMM_WORLD);


  //======================================== Apply partitioning scheme
  std::vector<int64_t> cell_pids;
  auto grid = chi_mesh::MeshContinuum::New();

  if (options.partition_type == PartitionType::KBA_STYLE_XYZ)
    cell_pids = KBA(*umesh);
  else
    cell_pids = PARMETIS(*umesh);

  //======================================== Load up the cells
  auto& vertex_subs = umesh->vertex_cell_subscriptions;
  size_t cell_globl_id = 0;
  for (auto raw_cell : umesh->raw_cells)
  {
    if (CellHasLocalScope(*raw_cell, cell_globl_id, vertex_subs, cell_pids))
    {
      auto cell = MakeCell(*raw_cell, cell_globl_id,
                           cell_pids[cell_globl_id], umesh->vertices);

      for (uint64_t vid : cell->vertex_ids)
        grid->vertices.Insert(vid, umesh->vertices[vid]);

      grid->cells.push_back(std::move(cell));
    }

    ++cell_globl_id;
  }//for raw_cell

  grid->SetGlobalVertexCount(umesh->vertices.size());

  chi::log.Log() << "Cells loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  SetContinuum(grid);
  SetGridAttributes(umesh->attributes, {umesh->mesh_options.ortho_Nx,
                                        umesh->mesh_options.ortho_Ny,
                                        umesh->mesh_options.ortho_Nz});

  //======================================== Concluding messages
  chi::log.LogAllVerbose1()
    << "### LOCATION[" << chi::mpi.location_id
    << "] amount of local cells="
    << grid->local_cell_glob_indices.size();

  size_t total_local_cells = grid->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(&total_local_cells,
                &total_global_cells,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi::log.Log()
    << "VolumeMesherPredefinedUnpartitioned: Cells created = "
    << total_global_cells
    << std::endl;

}