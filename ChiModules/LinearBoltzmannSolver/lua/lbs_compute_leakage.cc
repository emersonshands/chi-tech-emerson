#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

//###################################################################
/**Computes the leakage for the specified groupset and boundary id.
 *
\param SolverIndex int Handle to the solver.
\param GroupSetHandle int Handle to the groupset.
\param BoundaryID int Id of the boundary for which leakage is to be computed.

\return The leakage on a per group basis.

\ingroup LuaLBS
\author Jan*/
int chiLBSComputeLeakage(lua_State *L)
{
  const std::string fname = "chiLBSComputeLeakage";
  const int num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = chi::GetStackItem<lbs::SteadySolver>(chi::solver_stack,
                                                          solver_handle,
                                                          fname);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int groupset_id = lua_tonumber(L, 2);
  const int boundary_id = lua_tonumber(L, 3);

  const auto leakage = lbs_solver.ComputeLeakage(groupset_id, boundary_id);

  //============================================= Push up the table
  lua_newtable(L);

  for (int i=0; i<static_cast<int>(leakage.size()); ++i)
  {
    lua_pushinteger(L, i+1);
    lua_pushnumber(L, leakage[i]);
    lua_settable(L, -3);
    std::cout << leakage[i] << " ";
  }
  std::cout << "\n";

  return 1;
}
