#include "ChiLua/chi_lua.h"
#include "ChiMesh/FieldFunctionInterpolation/Point/chi_ffinter_point.h"
#include "ChiMesh/FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "ChiMesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define dcastPoint(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationPoint&>(x)
#define dcastLine(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationLine&>(x)
#define dcastVolume(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationVolume&>(x)

//#############################################################################
/** Gets the value(s) associated with an interpolation provided the
 * interpolation type has an associated value.
 *
\param FFIHandle int Handle to the field function interpolation.

###Note:
Currently only the POINT, LINE and VOLUME interpolation supports obtaining a
value. For the POINT and VOLUME types a single value is returned. For the LINE
type a table of tables is returned with the first index being the field function
(in the order it was assigned) and the second index being the point index.

\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationGetValue(lua_State *L)
{
  const std::string fname = __FUNCTION__;

  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiFFInterpolationGetValue",1,num_args);

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, fname);

  if (p_ffi->Type() == chi_mesh::ff_interpolation::Type::POINT)
  {
    auto& cur_ffi_point = dcastPoint(*p_ffi);
    double value = cur_ffi_point.GetPointValue();

    lua_pushnumber(L,value);
    return 1;
  }
  else if (p_ffi->Type() == chi_mesh::ff_interpolation::Type::LINE)
  {
    auto& cur_ffi_line = dcastLine(*p_ffi);

    lua_newtable(L);

    for (int ff=0; ff<cur_ffi_line.field_functions.size(); ff++)
    {
      lua_pushnumber(L,ff+1);

      lua_newtable(L);
      const auto& ff_ctx = cur_ffi_line.ff_contexts[ff];

      for (int p=0; p<cur_ffi_line.interpolation_points.size(); p++)
      {
        lua_pushnumber(L,p+1);
        lua_pushnumber(L,ff_ctx.interpolation_points_values[p]);
        lua_settable(L,-3);
      }

      lua_settable(L,-3);
    }

    return 1;
  }
  else if ( p_ffi->Type() == chi_mesh::ff_interpolation::Type::VOLUME)
  {
    auto& cur_ffi_volume = dcastVolume(*p_ffi);
    double value = cur_ffi_volume.op_value;

    lua_pushnumber(L,value);
    return 1;
  }
  else
  {
    chi::log.Log0Warning()
      << "chiFFInterpolationGetValue is currently only supported for "
      << " POINT, LINE and VOLUME interpolator types.";
  }

  return 0;
}
