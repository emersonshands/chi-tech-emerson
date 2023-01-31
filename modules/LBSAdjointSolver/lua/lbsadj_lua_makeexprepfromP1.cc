#include "ChiLua/chi_lua.h"

#include "LBSAdjointSolver/lbs_adjoint.h"

#include <stdexcept>

namespace lbs_adjoint::lua_utils
{
int chiAdjointSolverMakeExpRepFromP1Moments(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 1);

  std::vector<double> P1;
  bool verbose = false;
  LuaPopulateVectorFrom1DArray(fname, L, 1, P1);

  if (P1.size() != 4)
    throw std::invalid_argument(fname + ": Supplied table argument must"
                                        " have 4 entries.");

  if (num_args == 2)
  {
    LuaCheckNilValue(fname, L, 2);
    LuaCheckBoolValue(fname, L, 2);

    verbose = lua_toboolean(L, 2);
  }

  auto solution =
    lbs_adjoint::MakeExpRepFromP1({P1[0],P1[1],P1[2],P1[3]}, verbose);

  lua_pushnumber(L, solution[0]);
  lua_pushnumber(L, solution[1]);
  return 2;
}
}//namespace lbs_adjoint