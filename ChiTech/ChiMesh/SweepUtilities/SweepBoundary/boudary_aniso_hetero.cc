#include "sweep_boundaries.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Returns a pointer to a heterogenous flux storage location.*/
double* chi_mesh::sweep_management::BoundaryIncidentHeterogenous::
  HeterogenousPsiIncoming(uint64_t cell_local_id,
                          int face_num,
                          int fi,
                          int angle_num,
                          int group_num,
                          int gs_ss_begin)
{
  if (local_cell_data.empty())
  {
    chi::log.LogAllError()
      << "HeterogenousPsiIncoming call made to a heterogeneous boundary "
         "with that information not yet set up.";
    exit(EXIT_FAILURE);
  }

  const size_t dof_offset = num_groups*angle_num + group_num;

//  return &local_cell_data[cell_local_id][face_num][fi][dof_offset];
  return &local_cell_data.at(cell_local_id).at(face_num).at(fi).at(dof_offset);
}

//###################################################################
/**Performs the setup for a particular quadrature.*/
void chi_mesh::sweep_management::BoundaryIncidentHeterogenous::
  Setup(const chi_mesh::MeshContinuum &grid,
        const chi_math::AngularQuadrature &quadrature)
{
  const size_t num_local_cells = grid.local_cells.size();
  local_cell_data.clear();
  local_cell_data.reserve(num_local_cells);

  std::vector<bool> cell_bndry_flags(num_local_cells,false);
  for (const auto& cell : grid.local_cells)
    for (const auto& face : cell.faces)
      if (not face.has_neighbor)
      {
        cell_bndry_flags[cell.local_id] = true;
        break;
      }

  size_t num_angles = quadrature.omegas.size();

  typedef std::pair<double, double> PhiTheta;

  std::vector<int>               angle_indices;
  std::vector<chi_mesh::Vector3> angle_vectors;
  std::vector<PhiTheta>          phi_theta_angles;
  std::vector<int>               group_indices;

  angle_indices.reserve(num_angles);
  angle_vectors.reserve(num_angles);
  phi_theta_angles.reserve(num_angles);
  group_indices.reserve(num_groups);

  int num_angles_int = static_cast<int>(num_angles);
  for (int n=0; n<num_angles_int; ++n)
    angle_indices.emplace_back(n);
  for (int n=0; n<num_angles_int; ++n)
    angle_vectors.emplace_back(quadrature.omegas[n]);
  for (int n=0; n<num_angles_int; ++n)
  {
    auto& abscissae = quadrature.abscissae[n];
    double phi   = abscissae.phi;
    double theta = abscissae.theta;
    phi_theta_angles.emplace_back(std::make_pair(phi, theta));
  }
  for (int g=0; g<static_cast<int>(num_groups); ++g)
    group_indices.emplace_back(g);

  for (const auto& cell : grid.local_cells)
  {
    if (cell_bndry_flags[cell.local_id])
    {
      CellData cell_data(cell.faces.size());

      for (size_t f=0; f<cell.faces.size(); ++f)
      {
        auto& face = cell.faces[f];
        size_t face_num_nodes = face.vertex_ids.size();
        FaceData face_data;

        if (not face.has_neighbor and face.neighbor_id == ref_boundary_id)
        {
          face_data.reserve(face_num_nodes);
          for (size_t i=0; i<face_num_nodes; ++i)
          {
            std::vector<double> face_node_data =
              boundary_function->Evaluate(cell.global_id,
                                          cell.material_id,
                                          f, i,
                                          grid.vertices[face.vertex_ids[i]],
                                          face.normal,
                                          angle_indices,
                                          angle_vectors,
                                          phi_theta_angles,
                                          group_indices);

            face_data.push_back(std::move(face_node_data));
          }//for face node-i
        }//bndry face

        cell_data[f] = std::move(face_data);
      }//for face f

      local_cell_data.push_back(std::move(cell_data));
    }//if bndry cell
    else
      local_cell_data.emplace_back();

  }//for cell
}