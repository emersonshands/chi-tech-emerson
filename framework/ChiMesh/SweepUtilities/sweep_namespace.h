#ifndef CHI_SWEEP_H
#define CHI_SWEEP_H

#include "../chi_mesh.h"
#include <set>

#include <memory>

namespace chi_graph
{
  class DirectedGraph;
}

//###################################################################
namespace chi_mesh
{
namespace sweep_management
{
  struct STDG;           ///< Global Sweep Plane Ordering
  struct SPLS;           ///< Sweep Plane Local Subgrid
  class  PRIMARY_FLUDS;  ///< Primary Flux Data Structure
  class  AUX_FLUDS;      ///< Auxiliary Flux Data Structure
  struct SPDS;           ///< Sweep Plane Data Structure

  class  SweepBuffer;
  class AngleSet;
  class AngleSetGroup;
  class  AngleAggregation;

  class SweepChunk;

  class SweepScheduler;

  void PopulateCellRelationships(
    chi_mesh::MeshContinuumPtr grid,
    const chi_mesh::Vector3& omega,
    std::set<int>& location_dependencies,
    std::set<int>& location_successors,
    std::vector<std::set<std::pair<int,double>>>& cell_successors);

  void CommunicateLocationDependencies(
    const std::vector<int>& location_dependencies,
    std::vector<std::vector<int>>& global_dependencies);

  void RemoveGlobalCyclicDependencies(
    chi_mesh::sweep_management::SPDS* sweep_order,
    chi_graph::DirectedGraph& TDG);

  void RemoveLocalCyclicDependencies(
    std::shared_ptr<SPDS> sweep_order,
    chi_graph::DirectedGraph& local_DG);

  std::shared_ptr<SPDS> CreateSweepOrder(const chi_mesh::Vector3& omega,
                                         const chi_mesh::MeshContinuumPtr& grid,
                                         bool cycle_allowance_flag=false);

  void PrintSweepOrdering(SPDS* sweep_order,
                          MeshContinuumPtr vol_continuum);

  enum class AngleSetStatus{
    NOT_FINISHED = 0,
    FINISHED = 1,
    RECEIVING = 2,
    READY_TO_EXECUTE = 3,
    EXECUTE = 4,
    NO_EXEC_IF_READY = 5,
    MESSAGES_SENT = 6,
    MESSAGES_PENDING = 7
  };
  typedef AngleSetStatus ExecutionPermission;
}
}


#endif //CHI_SWEEP_H
