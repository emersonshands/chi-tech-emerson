#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

//###################################################################
/**Computes and prints the balance tables.
 *
\param SolverIndex int Handle to the solver.

\ingroup LuaLBS
\author Jan*/
int chiLBSComputeBalance(lua_State *L)
{
  const std::string fname = "chiLBSComputeBalance";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = chi::GetStackItem<lbs::SteadySolver>(chi::solver_stack,
                                                          solver_handle,
                                                          fname);

  lbs_solver.ComputeBalance();

  return 0;
}
