#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "CellViews/fv_cellbase.h"

void chi_math::SpatialDiscretization_FV::CreateCellMappings()
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization_FV::"
                                     "CreateCellMappings";

  auto MakeCellMapping = [this, fname](const chi_mesh::Cell& cell)
  {
    using namespace std;
    using namespace chi_math;
    std::unique_ptr<chi_math::CellMapping> mapping;

    switch (cell.Type())
    {
      case chi_mesh::CellType::SLAB:
      case chi_mesh::CellType::POLYGON:
      case chi_mesh::CellType::POLYHEDRON:
      {
        typedef std::vector<std::vector<int>> FaceDofMapping;
        mapping = make_unique<CellFVValues>(
          ref_grid,cell,cell.centroid,
          FaceDofMapping(cell.faces.size(),{-1}));
        break;
      }
      default:
        throw std::logic_error(std::string(fname) +
        std::string(": Invalid cell type encountered."));
    }
    return mapping;
  };

  for (const auto& cell : ref_grid->local_cells)
    cell_mappings.push_back(MakeCellMapping(cell));

  const auto ghost_ids = ref_grid->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto ghost_mapping = MakeCellMapping(ref_grid->cells[ghost_id]);
    nb_cell_mappings.insert(std::make_pair(ghost_id, std::move(ghost_mapping)));
  }
}