#include "chi_ffinter_volume.h"
#include "ChiMesh/Cell/cell.h"


#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Initializes the volume field function interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Initialize()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing volume interpolator.";
  //================================================== Check grid available
  if (field_functions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unassigned field function in volume field function interpolator.";
    exit(EXIT_FAILURE);
  } else
  {
    this->grid_view = field_functions[0]->grid;
  }

  //================================================== Find cell inside volume
  for (const auto& cell : grid_view->local_cells)
  {
    int cell_local_index = cell.local_id;

    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell.centroid);

    if (inside_logvolume)
    {
      for (int i=0; i < cell.vertex_ids.size(); i++)
      {
        cfem_local_nodes_needed_unmapped.push_back(cell.vertex_ids[i]);
        pwld_local_nodes_needed_unmapped.push_back(i);
        pwld_local_cells_needed_unmapped.push_back(cell_local_index);
      }//for dof
    }//if inside logicalVol

  }//for local cell
}