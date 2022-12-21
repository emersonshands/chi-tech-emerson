#include "lbts_transient_solver.h"

//###################################################################
/**Performs a timestep of the precursors.*/
void lbs::TransientSolver::StepPrecursors()
{
  const auto& BackwardEuler = chi_math::SteppingMethod::BACKWARD_EULER;
  const auto& CrankNicolson = chi_math::SteppingMethod::CRANK_NICHOLSON;

  double theta;
  if      (method == BackwardEuler) theta = 1.0;
  else if (method == CrankNicolson) theta = 0.5;
  else                              theta = 0.7;

  const double eff_dt = theta*dt;

  //============================================= Clear destination vector
  precursor_new_local.assign(precursor_new_local.size(), 0.0);

  //================================================== Loop over local cells
  // Uses phi_new and precursor_prev_local to compute
  // precursor_new_local(theta-flavor)
  for (auto& cell : grid->local_cells)
  {
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    const auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.Volume();

    //==================== Obtain xs
    const auto& xs = matid_to_xs_map.at(cell.material_id);

    //======================================== Compute delayed fission rate
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_view.NumNodes(); ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.Vi_vectors[i]/cell_volume;

      for (int g = 0; g < groups.size(); ++g)
        delayed_fission += xs->nu_delayed_sigma_f[g] *
                           phi_new_local[uk_map + g] *
                           node_V_fraction;
    }

    //========================================= Loop over precursors
    const auto& J = max_precursors_per_material;
    for (int j = 0; j < xs->num_precursors; ++j)
    {
      const size_t dof_map = cell.local_id * J + j;
      const double coeff = 1.0 / (1.0 + eff_dt*xs->precursor_lambda[j]);

      //contribute last time step precursors
      precursor_new_local[dof_map] = coeff * precursor_prev_local[dof_map];

      //contribute delayed fission production
      precursor_new_local[dof_map] +=
        coeff * eff_dt * xs->precursor_yield[j] * delayed_fission;
    }
  }//for cell

  //======================================== Compute t^{n+1} value
  {
    auto& Cj = precursor_new_local;
    const auto& Cj_prev = precursor_prev_local;

    const double inv_theta = 1.0/theta;
    for (size_t i = 0; i < Cj.size(); ++i)
      Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
  }
}