#include "lbs_linear_boltzmann_solver.h"
#include "Tools/boundary_func_lua.h"

//###################################################################
/**Initializes transport related boundaries. */
void lbs::SteadySolver::InitializeBoundaries()
{
  //================================================== Initialize default
  //                                                   incident boundary
  using namespace chi_mesh::sweep_management;
  std::vector<std::vector<double>>& flux_vec = incident_P0_mg_boundaries;

  // ================================================= Populate boundaries
  const size_t G = groups.size();
  if (sweep_boundaries.empty())
  {
    chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
    chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
    chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
    int bndry_id=0;
    for (auto bndry_type : boundary_types)
    {
      const int vec_index = bndry_type.second;

      if (bndry_type.first == lbs::BoundaryType::VACUUM)
        sweep_boundaries.emplace_back(
          new BoundaryIsotropicHomogenous(G,std::vector<double>(G,0.0)));
      else if (bndry_type.first == lbs::BoundaryType::INCIDENT_ISOTROPIC)
        sweep_boundaries.emplace_back(
          new BoundaryIsotropicHomogenous(G,flux_vec[vec_index]));
      else if (bndry_type.first == lbs::BoundaryType::REFLECTING)
      {
        chi_mesh::Normal normal;
        if (bndry_id == 0) normal = ihat;
        if (bndry_id == 1) normal = ihat*-1.0;
        if (bndry_id == 2) normal = jhat;
        if (bndry_id == 3) normal = jhat*-1.0;
        if (bndry_id == 4) normal = khat;
        if (bndry_id == 5) normal = khat*-1.0;

        sweep_boundaries.emplace_back(new BoundaryReflecting(G, normal));
      }
      else if (bndry_type.first == BoundaryType::INCIDENT_ANISTROPIC_HETEROGENOUS)
      {
        const auto& bndry_pref = boundary_preferences.at(bndry_id);
        const auto& lua_fname = bndry_pref.source_function;

        sweep_boundaries.emplace_back(new BoundaryIncidentHeterogenous(G,
                 std::make_unique<BoundaryFunctionToLua>(lua_fname), bndry_id));
      }

      ++bndry_id;
    }


  }//if empty

}