#include "chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_log.h"

//###################################################################
/**Obtains a pointer to the last created grid. This method will
 * get a smart-pointer to a grid object. If a volume-mesher has not
 * been created, or if a grid is not available, this method will
 * throw `std::logic_error`.*/
chi_mesh::MeshContinuumPtr& chi_mesh::MeshHandler::GetGrid() const
{
  if (volume_mesher == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetGrid: Volume mesher "
                           "undefined. This usually means a grid is not defined"
                           " or is incomplete.");

  auto& grid_ptr = volume_mesher->GetContinuum();

  if (grid_ptr == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetGrid: Volume mesher has "
                           "no grid available. Make sure the volume mesher has "
                           "been executed.");

  return grid_ptr;
}