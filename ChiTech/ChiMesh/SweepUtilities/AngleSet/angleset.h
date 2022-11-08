#ifndef CHI_ANGLESET_H
#define CHI_ANGLESET_H

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/SweepUtilities/SweepBuffer/sweepbuffer.h"
#include "ChiMesh/SweepUtilities/SweepBoundary/sweep_boundaries.h"
#include "ChiMesh/SweepUtilities/FLUDS/FLUDS.h"

#include <chi_mpi.h>

typedef chi_mesh::sweep_management::BoundaryBase SweepBndry;

#include <memory>

//###################################################################
/**Manages the workstages of a single angle set.*/
class chi_mesh::sweep_management::AngleSet
{
private:
  size_t                            num_grps;
  std::shared_ptr<SPDS>             spds;
  bool                              executed = false;

  chi_mesh::sweep_management::SweepBuffer sweep_buffer;

public:
  FLUDS*                            fluds;
  std::vector<size_t>               angles;
  std::vector<std::shared_ptr<SweepBndry>>&         ref_boundaries;
  size_t                            ref_subset;

  //FLUDS
  std::vector<std::vector<double>>  local_psi;
  std::vector<double>               delayed_local_psi;
  std::vector<double>               delayed_local_psi_old;
  std::vector<std::vector<double>>  deplocI_outgoing_psi;
  std::vector<std::vector<double>>  prelocI_outgoing_psi;
  std::vector<std::vector<double>>  boundryI_incoming_psi;

  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi;
  std::vector<std::vector<double>>  delayed_prelocI_outgoing_psi_old;

  AngleSet(size_t in_numgrps,
           size_t in_ref_subset,
           std::shared_ptr<SPDS>& in_spds,
           FLUDS* in_fluds,
           std::vector<size_t>& angle_indices,
           std::vector<std::shared_ptr<SweepBndry>>& sim_boundaries,
           int sweep_eager_limit,
           chi_objects::ChiMPICommunicatorSet* in_comm_set);

  void InitializeDelayedUpstreamData();

  std::shared_ptr<chi_mesh::sweep_management::SPDS> GetSPDS();

  int GetMaxBufferMessages() const;

  void SetMaxBufferMessages(int new_max);

  size_t GetNumGrps() const;

  AngleSetStatus AngleSetAdvance(
             SweepChunk& sweep_chunk,
             int angle_set_num,
             const std::vector<size_t>& timing_tags,
             ExecutionPermission permission = ExecutionPermission::EXECUTE);
  AngleSetStatus FlushSendBuffers();
  void ResetSweepBuffers();
  void ReceiveDelayedData(int angle_set_num);

  double* PsiBndry(uint64_t bndry_map,
                   int angle_num,
                   uint64_t cell_local_id,
                   int face_num,
                   int fi,
                   int g,
                   int gs_ss_begin,
                   bool surface_source_active);
  double* ReflectingPsiOutBoundBndry(uint64_t bndry_map,
                                     int angle_num,
                                     uint64_t cell_local_id,
                                     int face_num,
                                     int fi,
                                     int gs_ss_begin);
};

#endif //CHI_ANGLESET_H