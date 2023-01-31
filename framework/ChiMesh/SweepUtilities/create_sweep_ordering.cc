#include "../chi_mesh.h"
#include "sweep_namespace.h"

#include "../MeshHandler/chi_meshhandler.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"

#include "chi_mpi.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiConsole/chi_console.h"
#include "ChiTimer/chi_timer.h"


;

#include "ChiGraph/chi_directed_graph.h"



//###################################################################
/**Develops a sweep ordering for a given angle for locally owned
 * cells.*/
std::shared_ptr<chi_mesh::sweep_management::SPDS>
chi_mesh::sweep_management::
  CreateSweepOrder(const chi_mesh::Vector3& omega,
                   const chi_mesh::MeshContinuumPtr& grid,
                   bool cycle_allowance_flag)
{
  auto sweep_order  = std::make_shared<chi_mesh::sweep_management::SPDS>();
  sweep_order->grid = grid;

  size_t num_loc_cells = grid->local_cell_glob_indices.size();

  //============================================= Assign direction vector
  sweep_order->omega = omega;

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Building sweep ordering for Omega = "
    << omega.PrintS();

  //============================================= Populate Cell Relationships
  chi::log.Log0Verbose1() << "Populating cell relationships";
  std::vector<std::set<std::pair<int,double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(grid,
                            omega,
                            location_dependencies,
                            location_successors,
                            cell_successors);

  sweep_order->location_successors.reserve(location_successors.size());
  sweep_order->location_dependencies.reserve(location_dependencies.size());

  for (auto v : location_successors)
    sweep_order->location_successors.push_back(v);

  for (auto v : location_dependencies)
    sweep_order->location_dependencies.push_back(v);

  //============================================= Build graph
  chi_graph::DirectedGraph local_DG;

  // Add vertex for each local cell
  for (int c=0; c<num_loc_cells; ++c)
    local_DG.AddVertex();

  // Create graph edges
  for (int c=0; c<num_loc_cells; c++)
    for (auto& successor : cell_successors[c])
      local_DG.AddEdge(c, successor.first, successor.second);

  //============================================= Remove local cycles if allowed
  if (cycle_allowance_flag)
  {
    chi::log.Log0Verbose1()
      << chi::program_timer.GetTimeString()
      << " Removing inter-cell cycles.";
    RemoveLocalCyclicDependencies(sweep_order,local_DG);
  }

  //============================================= Generate topological sorting
  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Generating topological sorting for local sweep ordering";
  auto so_temp = local_DG.GenerateTopologicalSort();
  sweep_order->spls.item_id.clear();
  for (auto v : so_temp)
    sweep_order->spls.item_id.emplace_back(v);

  if (sweep_order->spls.item_id.empty())
  {
    chi::log.LogAllError()
      << "Topological sorting for local sweep-ordering failed. "
      << "Cyclic dependencies detected. Cycles need to be allowed"
      << " by calling application.";
    chi::Exit(EXIT_FAILURE);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Task
  //                                                        Dependency Graphs
  //All locations will gather other locations' dependencies
  //so that each location has the ability to build
  //the global task graph.

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Communicating sweep dependencies.";

  auto& global_dependencies = sweep_order->global_dependencies;
  global_dependencies.resize(chi::mpi.process_count);

  CommunicateLocationDependencies(sweep_order->location_dependencies,
                                  global_dependencies);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Build task
  //                                                        dependency graph
  sweep_order->BuildTaskDependencyGraph(cycle_allowance_flag);

  MPI_Barrier(MPI_COMM_WORLD);

  chi::log.Log0Verbose1()
    << chi::program_timer.GetTimeString()
    << " Done computing sweep ordering.\n\n";

  return sweep_order;
}
