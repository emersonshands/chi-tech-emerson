#include "lbs_linear_boltzmann_solver.h"

//###################################################################
/**Initializes common groupset items.*/
void lbs::SteadySolver::InitializeGroupsets()
{
  for (auto& groupset : groupsets)
  {
    //================================================== Build groupset angular
    //                                                   flux unknown manager
    groupset.psi_uk_man.unknowns.clear();
    const size_t gs_num_angles = groupset.quadrature->abscissae.size();
    const size_t gs_num_groups = groupset.groups.size();
    auto& grpset_psi_uk_man = groupset.psi_uk_man;

    for (unsigned int n=0; n<gs_num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N, gs_num_groups);

    groupset.BuildDiscMomOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildMomDiscOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildSubsets();
  }//for groupset
}