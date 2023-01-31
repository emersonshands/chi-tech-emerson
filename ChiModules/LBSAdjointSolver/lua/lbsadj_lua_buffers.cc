#include "ChiLua/chi_lua.h"

#include "lbsadj_lua_utils.h"

namespace lbs_adjoint::lua_utils
{

/**Reads flux-moments file to a buffer and returns a handle to that buffer.

\param SolverHandle int Handle to the relevant solver.
\param FileBaseName string The base-name of the file(s) from which to read
                           the flux moments.

\return handle int A handle that can be used with
                   `chiAdjointSolverApplyFluxMomentBuffer`.*/
int chiAdjointSolverReadFluxMomentsToBuffer(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle     = lua_tointeger(L, 1);

  auto& solver = chi::GetStackItem<lbs_adjoint::AdjointSolver>(
    chi::solver_stack, solver_handle, fname);

  const std::string file_basename = lua_tostring(L,2);

  std::vector<double> moments;
  solver.ReadFluxMoments(file_basename, moments);

  solver.m_moment_buffers.push_back(std::move(moments));

  const size_t handle = solver.m_moment_buffers.size()-1;
  lua_pushinteger(L, static_cast<lua_Integer>(handle));

  return 1;
}

/**Applies buffered flux-moments data to the current phi-old.

\param SolverHandle int Handle to the relevant solver.
\param BufferHandle int The handle to the buffer-position to be applied.
 */
int chiAdjointSolverApplyFluxMomentBuffer(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle     = lua_tointeger(L, 1);

  auto& solver = chi::GetStackItem<lbs_adjoint::AdjointSolver>(
    chi::solver_stack, solver_handle, fname);

  const int buffer_handle = lua_tointeger(L,2);

  if (buffer_handle < 0 or
      buffer_handle >= solver.m_moment_buffers.size())
    throw std::invalid_argument(fname + ": Invalid buffer handle.");

  solver.phi_old_local = solver.m_moment_buffers[buffer_handle];

  return 0;
}

}//namespace lbs_adjoint::lua_utils