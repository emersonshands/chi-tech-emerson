#include "lbs_linear_boltzmann_solver.h"

/**Computes the angular flux based leakage from boundary surfaces.
\param groupset_id The groupset for which to compute the leakage.
\param boundary_id uint64_t The boundary-id for which to perform the integration.

\return The leakage as a value.
*/
std::vector<double> lbs::SteadySolver::
  ComputeLeakage(const int groupset_id, const uint64_t boundary_id) const
{
  const std::string fname = "lbs::SteadySolver::ComputeLeakage";

  //================================================== Perform checks
  if (groupset_id<0 or groupset_id>=groupsets.size())
    throw std::invalid_argument(fname + ": Invalid groupset_id specified.");

  if (not options.save_angular_flux)
    throw std::logic_error(fname + ": Requires options.save_angular_flux to be"
                                   " true.");

  //================================================== Get info
  const auto& sdm = *discretization;
  const auto& groupset = groupsets.at(groupset_id);
  const auto& psi_uk_man = groupset.psi_uk_man;
  const auto& quadrature = groupset.quadrature;

  const size_t num_angles = quadrature->omegas.size();

  const int gsi = groupset.groups.front().id;
  const int gsf = groupset.groups.back().id;
  const int gs_num_groups = gsf+1-gsi;

  //================================================== Start integration
  std::vector<double> leakage(gs_num_groups, 0.0);
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto& fe_values = unit_cell_matrices[cell.local_id];

    size_t f=0;
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor and face.neighbor_id == boundary_id)
      {
        const auto& IntF_shapeI = fe_values.face_Si_vectors[f];
        const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
        for (size_t fi=0; fi<num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          for (size_t n=0; n<num_angles; ++n)
          {
            const auto& omega  = quadrature->omegas[n];
            const auto& weight = quadrature->weights[n];
            const double mu = omega.Dot(face.normal);
            if (mu > 0.0)
            {
              for (int gi=0; gi<gs_num_groups; ++gi)
              {
                const int g = gi+gsi;
                const int64_t imap = sdm.MapDOFLocal(cell, i, psi_uk_man, n, g);

                const double psi = psi_new_local[groupset_id][imap];

                leakage[gi] += weight * mu * psi * IntF_shapeI[i];
              }//for g
            }//outgoing
          }//for n
        }//for face node
      }//if right bndry
      ++f;
    }//for face
  }//for cell

  return leakage;
}