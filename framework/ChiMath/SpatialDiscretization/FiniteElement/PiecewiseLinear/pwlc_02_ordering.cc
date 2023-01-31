#include "pwlc.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include "ChiTimer/chi_timer.h"
#include "ChiMPI/chi_mpi_utils_map_all2all.h"

#include <algorithm>

//###################################################################
/**Reorders the nodes for parallel computation in a Continuous
 * Finite Element calculation.*/
void chi_math::SpatialDiscretization_PWLC::OrderNodes()
{
  const std::string fname = __FUNCTION__;
  //============================================= Build set of local scope nodes
  // ls_node_id = local scope node id
  std::set<uint64_t> ls_node_ids_set;
  for (const auto& cell : ref_grid->local_cells)
    for (uint64_t node_id : cell.vertex_ids)
      ls_node_ids_set.insert(node_id);

  //============================================ Build node partition
  //                                             subscriptions
  // psub = partition subscription
  // Multiple partitions can subscribe to a given
  // node. We build this list here.
  // We start by adding the current location id
  // as the first subscription
  typedef std::set<uint64_t> PSUBS;
  std::map<uint64_t, PSUBS> ls_node_ids_psubs;
  for (const uint64_t node_id : ls_node_ids_set)
    ls_node_ids_psubs[node_id] = {static_cast<uint64_t>(chi::mpi.location_id)};

  // Now we add the partitions assoc. with the
  // ghost cells.
  const auto ghost_cell_ids = ref_grid->cells.GetGhostGlobalIDs();
  for (const uint64_t ghost_id : ghost_cell_ids)
  {
    const auto& ghost_cell = ref_grid->cells[ghost_id];
    for (const uint64_t vid : ghost_cell.vertex_ids)
      ls_node_ids_psubs[vid].insert(ghost_cell.partition_id);
  }//for ghost_id

  //============================================= Build lists of local- and
  //                                              non-local nodes
  // The lowest partition-# owns a node.
  std::vector<uint64_t> local_node_ids;
  std::map<uint64_t, std::vector<uint64_t>> nonlocal_node_ids_map;
  for (const uint64_t node_id : ls_node_ids_set)
  {
    uint64_t smallest_partition_id = chi::mpi.location_id;
    for (const uint64_t pid : ls_node_ids_psubs[node_id]) //pid = partition id
      smallest_partition_id = std::min(smallest_partition_id, pid);

    if (smallest_partition_id == chi::mpi.location_id)
      local_node_ids.push_back(node_id);
    else
      nonlocal_node_ids_map[smallest_partition_id].push_back(node_id);
  }

  //============================================= Communicate node counts
  const uint64_t local_num_nodes = local_node_ids.size();
  locJ_block_size.assign(chi::mpi.process_count,0);
  MPI_Allgather(&local_num_nodes,       // sendbuf
                1, MPI_UINT64_T,        // sendcount, sendtype
                locJ_block_size.data(), // recvbuf
                1, MPI_UINT64_T,        // recvcount, recvtype
                MPI_COMM_WORLD);        // comm

  //============================================= Build block addresses
  locJ_block_address.assign(chi::mpi.process_count,0);
  uint64_t global_num_nodes = 0;
  for (int j=0; j<chi::mpi.process_count; ++j)
  {
    locJ_block_address[j] = global_num_nodes;
    global_num_nodes += locJ_block_size[j];
  }

  local_block_address = locJ_block_address[chi::mpi.location_id];

  local_base_block_size = local_num_nodes;
  globl_base_block_size = global_num_nodes;

  //============================================= Build node mapping for local
  //                                              nodes
  node_mapping.clear();
  for (uint64_t i=0; i<local_num_nodes; ++i)
    node_mapping[local_node_ids[i]] =
      static_cast<int64_t>(local_block_address + i);

  //============================================= Communicate nodes in need
  //                                              of mapping
  std::map<uint64_t, std::vector<uint64_t>> query_node_ids =
    chi_mpi_utils::MapAllToAll(nonlocal_node_ids_map, MPI_UINT64_T);

  //============================================= Map the query nodes
  std::map<uint64_t, std::vector<int64_t>> mapped_node_ids;
  for (const auto& key_value : query_node_ids)
  {
    const uint64_t& pid                    = key_value.first;
    const auto& node_list = key_value.second;

    for (const uint64_t node_id : node_list)
      if (node_mapping.count(node_id) == 0)
        throw std::logic_error("Error mapping query node.");
      else
      {
        const int64_t mapping = node_mapping.at(node_id);
        mapped_node_ids[pid].push_back(mapping);
      }
  }//for query location and nodes

  //============================================= Communicate back the mappings
  std::map<uint64_t, std::vector<int64_t>> nonlocal_node_ids_map_mapped =
    chi_mpi_utils::MapAllToAll(mapped_node_ids, MPI_INT64_T);

  //============================================= Processing the mapping for
  //                                              non-local nodes
  m_ghost_node_mapping.clear();
  try
  {
    for (const auto& pid_node_ids : nonlocal_node_ids_map)
    {
      const uint64_t& pid       = pid_node_ids.first;
      const auto&     node_list = pid_node_ids.second;
      const auto&     mappings  = nonlocal_node_ids_map_mapped.at(pid);

      if (mappings.size() != node_list.size())
        throw std::logic_error("mappings.size() != node_list.size()");

      const size_t num_nodes = node_list.size();
      for (size_t i=0; i<num_nodes; ++i)
      {
        node_mapping[node_list[i]] = mappings[i];
        m_ghost_node_mapping[node_list[i]] = mappings[i];
      }
    }//for pid and non-local id
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(fname + ": Processing non-local mapping failed.");
  }
  catch (const std::logic_error& lerr)
  {
    throw std::logic_error(fname + ": Processing non-local mapping failed." +
                           lerr.what());
  }


}