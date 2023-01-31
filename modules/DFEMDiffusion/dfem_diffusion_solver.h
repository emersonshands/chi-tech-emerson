#ifndef DFEM_DIFFUSION_SOLVER_H
#define DFEM_DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "dfem_diffusion_bndry.h"
#include "ChiTimer/chi_timer.h"

#include "ChiConsole/chi_console.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "ChiMesh/chi_mesh.h"

#include <map>

// forward declaration
namespace chi_mesh
{
class MeshContinuum; 
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
}
namespace chi_math
{
class SpatialDiscretization; 
typedef std::shared_ptr<SpatialDiscretization> SDMPtr ;
}

namespace dfem_diffusion
{
/** DFEM diffusion solver
 * 
*/
class Solver : public chi_physics::Solver
{
public:
  chi_mesh::MeshContinuumPtr grid_ptr = nullptr;

  chi_math::SDMPtr sdm_ptr = nullptr;

  size_t num_local_dofs = 0;
  size_t num_globl_dofs = 0;

  std::vector<double> field;

  Vec            x = nullptr;            // approx solution
  Vec            b = nullptr;            // RHS
  Mat            A = nullptr;            // linear system matrix

  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::map<uint, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences      boundary_preferences;
  std::vector<Boundary>   boundaries;

  explicit Solver(const std::string& in_solver_name);
  ~Solver() override;

  // void Initialize() override;
  void Initialize() override;

  void Execute() override;

  double HPerpendicular(const chi_mesh::Cell& cell, unsigned int f);

  int MapFaceNodeDisc(const chi_mesh::Cell& cur_cell,
                      const chi_mesh::Cell& adj_cell,
                      const std::vector<chi_mesh::Vector3>& cc_node_locs,
                      const std::vector<chi_mesh::Vector3>& ac_node_locs,
                      size_t ccf, size_t acf,
                      size_t ccfi,
                      double epsilon=1.0e-12);

  static
  double CallLua_iXYZFunction(lua_State* L,
                              const std::string&,
                              int,
                              const chi_mesh::Vector3&);

  void UpdateFieldFunctions();
};

} // namespace dfem_diffusion


#endif //DFEM_DIFFUSION_SOLVER_H

