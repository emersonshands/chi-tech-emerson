#include "diffusion_solver.h"

#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/PhysicsMaterial/material_property_scalarvalue.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "chi_runtime.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"


//###################################################################
/**Gets material properties various sources.*/
void chi_diffusion::Solver::GetMaterialProperties(const chi_mesh::Cell& cell,
                                                  int cell_dofs,
                                                  std::vector<double>& diffCoeff,
                                                  std::vector<double>& sourceQ,
                                                  std::vector<double>& sigmaa,
                                                  int group,
                                                  int moment)
{
  uint64_t cell_glob_index = cell.global_id;
  bool cell_is_local = (cell.partition_id == chi::mpi.location_id);
  uint64_t cell_local_id = cell.local_id;
  int mat_id = cell.material_id;

  if (mat_id<0)
  {
    chi::log.Log0Error()
      << "Cell encountered with no material id. ";
    chi::Exit(EXIT_FAILURE);
  }

  if (mat_id>=chi::material_stack.size())
  {
    chi::log.Log0Error()
      << "Cell encountered with material id pointing to "
         "non-existing material.";
    chi::Exit(EXIT_FAILURE);
  }

  auto property_map_D     = basic_options("property_map_D").IntegerValue();
  auto property_map_q     = basic_options("property_map_q").IntegerValue();
  auto property_map_sigma = basic_options("property_map_sigma").IntegerValue();

  auto material = chi::GetStackItemPtr(chi::material_stack, mat_id, __FUNCTION__);

  //====================================== Process material properties
  diffCoeff.resize(cell_dofs,1.0);
  sourceQ.resize(cell_dofs,0.0);
  sigmaa.resize(cell_dofs,0.0);

  //####################################################### REGULAR MATERIAL
  if (material_mode == DIFFUSION_MATERIALS_REGULAR)
  {
    //We absolutely need the diffusion coefficient so process error
    if ((property_map_D < 0) || (property_map_D >= material->properties.size()))
    {
      chi::log.Log0Error()
        << "Solver diffusion coefficient mapped to property index "
        << property_map_D << " is not a valid index for material \""
        << material->name <<"\" id " << mat_id;
      chi::Exit(EXIT_FAILURE);
    }

    //For now we can only support scalar values so lets check that
    if (std::dynamic_pointer_cast<chi_physics::ScalarValue>
        (material->properties[property_map_D]))
    {
      diffCoeff.assign(cell_dofs,
                       material->properties[property_map_D]->GetScalarValue());
    }
    else
    {
      chi::log.Log0Error()
        << "Solver diffusion coefficient mapped to property index "
        << property_map_D << " is not a valid property type"
        << " for material \""
        << material->name <<"\" id " << mat_id
        << ". Currently SCALAR_VALUE and THERMAL_CONDUCTIVITY are the "
        << "only supported types.";
      chi::Exit(EXIT_FAILURE);
    }


    if ((property_map_q < material->properties.size()) &&
        (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<chi_physics::ScalarValue>
          (material->properties[property_map_q]))
      {
        sourceQ.assign(cell_dofs,
                       material->properties[property_map_q]->GetScalarValue());
      }
      else
      {
        chi::log.Log0Error()
          << "Source value mapped to property index "
          << property_map_q << " is not a valid property type"
          << " for material \""
          << material->name <<"\" id " << mat_id
          << ". Currently SCALAR_VALUE is the "
          << "only supported type.";
        chi::Exit(EXIT_FAILURE);
      }
    }

    if (not ((property_map_sigma < 0) ||
             (property_map_sigma >= material->properties.size())))
    {
      sigmaa.assign(cell_dofs,
                    material->properties[property_map_sigma]->GetScalarValue());
    }
  }//regular

  //####################################################### TRANSPORT XS D
  //                                                        TRANSPORT XS SIGA
  //                                                        SCALAR       Q
  else if (material_mode == DIFFUSION_MATERIALS_FROM_TRANSPORTXS_TTR)
  {
    //====================================== Setting D and Sigma_a
    bool transportxs_found = false;
    for (int p=0; p<material->properties.size(); p++)
    {
      if (std::dynamic_pointer_cast<chi_physics::TransportCrossSections>
          (material->properties[p]))
      {
        auto xs = std::static_pointer_cast<chi_physics::TransportCrossSections>(
          material->properties[p]);

        if (!xs->diffusion_initialized)
          xs->ComputeDiffusionParameters();

        diffCoeff.assign(cell_dofs,xs->diffusion_coeff[group]);
        sigmaa.assign(cell_dofs,xs->sigma_removal[group]);
        transportxs_found = true;
      }
    }//for properties

    if (!transportxs_found)
    {
      chi::log.LogAllError()
        << "Diffusion Solver: Material encountered with no tranport xs"
           " yet material mode is DIFFUSION_MATERIALS_FROM_TRANSPORTXS.";
      chi::Exit(EXIT_FAILURE);
    }

    //====================================== Setting Q
    if ((property_map_q < material->properties.size()) &&
        (property_map_q >= 0))
    {
      if (std::dynamic_pointer_cast<chi_physics::ScalarValue>
          (material->properties[property_map_q]))
      {
        sourceQ.assign(cell_dofs,
                       material->properties[property_map_q]->GetScalarValue());
      }
      else
      {
        chi::log.Log0Error()
          << "Source value mapped to property index "
          << property_map_q << " is not a valid property type"
          << " for material \""
          << material->name <<"\" id " << mat_id
          << ". Currently SCALAR_VALUE is the "
          << "only supported type.";
        chi::Exit(EXIT_FAILURE);
      }
    }
  }//transport xs TTR
  else
  {
    chi::log.Log0Error()
      << "Diffusion Solver: Invalid material mode.";
    chi::Exit(EXIT_FAILURE);
  }


}


//###################################################################
/**Update the field functions with the latest data.*/
void chi_diffusion::Solver::UpdateFieldFunctions()
{
  chi::log.LogAll() << "Updating field functions" << std::endl;
  auto& ff = *field_functions.front();
  const auto& OneDofPerNode = discretization->UNITARY_UNKNOWN_MANAGER;

  std::vector<double> data_vector;
  discretization->LocalizePETScVector(x, data_vector, OneDofPerNode);

  ff.UpdateFieldVector(data_vector);
}
