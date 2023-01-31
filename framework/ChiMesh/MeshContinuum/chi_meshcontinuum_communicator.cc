#include "chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_runtime.h"
#include "chi_mpi.h"

//###################################################################
/***/
chi_objects::ChiMPICommunicatorSet& chi_mesh::MeshContinuum::GetCommunicator()
{
  //================================================== Check if already avail
  if (communicators_available)
    return communicator_set;

  //================================================== Build the communicator
  chi::log.Log0Verbose1() << "Building communicator.";
  std::set<int>    local_graph_edges;
  std::vector<int> local_connections;

  //================================================== Loop over local cells
  //Populate local_graph_edges
  local_graph_edges.insert(chi::mpi.location_id); //add current location
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor)
        if (not face.IsNeighborLocal(*this))
          local_graph_edges.insert(face.GetNeighborPartitionID(*this));
    }//for f
  }//for local cells

  //============================================= Convert set to vector
  //This is just done for convenience
  std::set<int>::iterator graph_edge;
  for (graph_edge =  local_graph_edges.begin();
       graph_edge != local_graph_edges.end();
       graph_edge++)
  {
    local_connections.push_back(*graph_edge);
  }

  //============================================= Broadcast local connection size
  chi::log.Log0Verbose1()
    << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(chi::mpi.process_count,
                                             std::vector<int>());
  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    int locI_num_connections = local_connections.size();

    //If chi::mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(&locI_num_connections,1,MPI_INT,locI,MPI_COMM_WORLD);

    if (chi::mpi.location_id != locI)
    {global_graph[locI].resize(locI_num_connections,-1);}
    else
    {
      std::copy(local_connections.begin(),
                local_connections.end(),
                std::back_inserter(global_graph[locI]));
    }
  }



  //============================================= Broadcast local connections
  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    //If chi::mpi.location_id == locI then this call will
    //act like a send instead of receive. Otherwise
    //It receives the count.
    MPI_Bcast(global_graph[locI].data(),
              global_graph[locI].size(),
              MPI_INT,locI,MPI_COMM_WORLD);
  }

  chi::log.Log0Verbose1()
    << "Done communicating local connections.";


  //============================================= Build groups
  MPI_Comm_group(MPI_COMM_WORLD,&communicator_set.world_group);
  communicator_set.location_groups.resize(chi::mpi.process_count, MPI_Group());

  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    MPI_Group_incl(communicator_set.world_group,
                   global_graph[locI].size(),
                   global_graph[locI].data(),
                   &communicator_set.location_groups[locI]);
  }

  //============================================= Build communicators
  chi::log.Log0Verbose1()
    << "Building communicators.";
  communicator_set.communicators.resize(chi::mpi.process_count, MPI_Comm());

  for (int locI=0;locI<chi::mpi.process_count; locI++)
  {
    int err = MPI_Comm_create_group(MPI_COMM_WORLD,
                                    communicator_set.location_groups[locI],
                                    0, //tag
                                    &communicator_set.communicators[locI]);

    if (err != MPI_SUCCESS)
    {
      chi::log.Log0Verbose1()
        << "Communicator creation failed.";
    }
  }

  chi::log.Log0Verbose1()
    << "Done building communicators.";

  communicators_available = true;
  return communicator_set;
}
