#include "ChiLua/chi_lua.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/FieldFunctionInterpolation/Point/chi_ffinter_point.h"
#include "ChiMesh/FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "ChiMesh/FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "ChiMesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define dcastPoint(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationPoint&>(x)
#define dcastLine(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationLine&>(x)
#define dcastSlice(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationSlice&>(x)
#define dcastVolume(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationVolume&>(x)

//#############################################################################
/** Creates a new field function interpolation.
 *
\param FFIHandle int Handle to the field function interpolation.
\param PropertyIndex int Type of property to set.

##_

###PropertyIndex\n
ADD_FIELDFUNCTION     = Add field function to interpolation.\n
SET_FIELDFUNCTIONS    = Sets the field functions to interpolate using a lua table.\n
PROBEPOINT            = Reference point for POINT type FFIs.\n
SLICE_POINT           = Reference point for SLICE type FFIs.\n
SLICE_NORMAL          = Normal of slice plane.\n
SLICE_TANGENT         = Tangent vector of slice plane.\n
SLICE_BINORM          = Binorm vector of slice plane.\n
LINE_FIRSTPOINT   = Line start point.\n
LINE_SECONDPOINT  = Line end point.\n
LINE_NUMBEROFPOINTS = Number of points to put in the line interpolator.
                          Minimum 2.\n
LINE_CUSTOM_ARRAY = Adds custom array to line interpolator.\n
OPERATION  =  Some interpolations support operation types. See OpTypes.\n
LOGICAL_VOLUME = To be followed by a handle to a logical volume to be
                 used by the interpolator.\n

###OpTypes
Basic operations are volume sum, volume average or volume max. The volume
sum is computed from
\f[
Sum = \sum_k \sum_i u_i \int_V N_i .dV.
\f]
The volume average is computed from
\f[
Avg = \frac{\sum_k \sum_i u_i \int_V N_i .dV}{\sum_k \sum_i \int_V N_i .dV}
\f]
The volume max is simply \f$ max_{k,i}(u_i) \f$.\n\n

OP_SUM\n
For volume interpolations, computes the volume integral.\n
\n
OP_AVG\n
For volume interpolations, computes the volume average.\n
OP_MAX\n
For volume interpolations, computes the volume max.\n
\n
A modified version of these operations are also available. Instead of OP_SUM,
OP_AVG and OP_MAX, the user may supply OP_SUM_LUA, OP_AVG_LUA and OP_MAX_LUA
which then needs to be followed by a string value `LuaFunctionName` of a lua
function of the following form:

\code
function LuaFunctionName(ff_value, mat_id)
 ret_val = 0.0;   --Or some computation
 return ret_val
end
\endcode

This code will be called to return a value \f$ f(u_i) \f$ to be used instead of
the field function \f$ u_i \f$.\n

Example:
\code
xwing=2.0
function IntegrateMaterialVolume(ff_value,mat_id)
    return xwing
end
ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM_LUA,"IntegrateMaterialVolume")
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
print(chiFFInterpolationGetValue(curffi))
\endcode

The code above will return 2.0 times the volume of cells included in the logical
volume `vol0`.


\return Handle int Handle to the created interpolation.
\ingroup LuaFFInterpol
\author Jan*/
int chiFFInterpolationSetProperty(lua_State *L)
{
  const std::string fname = "chiFFInterpolationSetProperty";
  int numArgs = lua_gettop(L);
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Get handle to field function
  const size_t ffihandle = lua_tonumber(L,1);

  auto p_ffi = chi::GetStackItemPtr(chi::field_func_interpolation_stack,
                                    ffihandle, fname);

  //================================================== Process properties
  using namespace chi_mesh::ff_interpolation;
  auto property = static_cast<Property>(lua_tonumber(L,2));
  //======================================== Check point properties
  if (property == Property::PROBEPOINT)
    if (p_ffi->Type() != chi_mesh::ff_interpolation::Type::POINT)
      throw std::logic_error(
        "Point property" + std::to_string(static_cast<int>(property)) +
        " used in chiFFInterpolationSetProperty but FFI is not a point-probe.");

  //======================================== Check slice properties
  if ((property >= Property::SLICEPOINT) && (property <= Property::SLICEBINORM))
    if (p_ffi->Type() != chi_mesh::ff_interpolation::Type::SLICE)
      throw std::logic_error(
        "Slice property" + std::to_string(static_cast<int>(property)) +
        " used in chiFFInterpolationSetProperty but FFI is not a slice.");

  //======================================== Check Line properties
  if ((property >= Property::FIRSTPOINT) && (property <= Property::NUMBEROFPOINTS))
    if (p_ffi->Type() != chi_mesh::ff_interpolation::Type::LINE)
      throw std::logic_error("Line property " +
      std::to_string(static_cast<int>(property)) +
        " used in chiFFInterpolationSetProperty but FFI is not a line.");

  //========================================= Generic
  if (property == Property::ADD_FIELD_FUNCTION)
  {
    int ffhandle = lua_tonumber(L,3);
    auto cur_ff = chi::GetStackItemPtr(chi::field_function_stack, ffhandle, fname);

    p_ffi->field_functions.push_back(cur_ff);
  }
  else if (property == Property::SET_FIELD_FUNCTIONS)
  {
    LuaCheckTableValue(fname, L, 3);
    std::vector<double> handle_array;
    LuaPopulateVectorFrom1DArray(fname, L, 3, handle_array);

    for (double handle_d : handle_array)
    {
      const auto ffhandle = static_cast<int>(handle_d);
      auto cur_ff = chi::GetStackItemPtr(chi::field_function_stack, ffhandle, fname);

      p_ffi->field_functions.push_back(cur_ff);
    }//for handle
  }
  else if (property == Property::PROBEPOINT)
  {
    auto& cur_ffi = dcastPoint(*p_ffi);

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi.m_point_of_interest = chi_mesh::Vector3(x, y, z);
  }
  else if (property == Property::SLICEPOINT)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.plane_point = chi_mesh::Vector3(x, y, z);
  }
  else if (property == Property::SLICENORMAL)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.normal = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.normal = cur_ffi_slice.normal/
                            cur_ffi_slice.normal.Norm();
  }
  else if (property == Property::SLICETANGENT)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.tangent = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.tangent = cur_ffi_slice.tangent/
                            cur_ffi_slice.tangent.Norm();
  }
  else if (property == Property::SLICEBINORM)
  {
    auto& cur_ffi_slice = dcastSlice(*p_ffi);

    double x = lua_tonumber(L,3);
    double y = lua_tonumber(L,4);
    double z = lua_tonumber(L,5);

    cur_ffi_slice.binorm = chi_mesh::Vector3(x, y, z);
    cur_ffi_slice.binorm = cur_ffi_slice.binorm/
                            cur_ffi_slice.binorm.Norm();
  }
  else if (property == Property::FIRSTPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    cur_ffi_line.pi.x = lua_tonumber(L,3);
    cur_ffi_line.pi.y = lua_tonumber(L,4);
    cur_ffi_line.pi.z = lua_tonumber(L,5);

  }
  else if (property == Property::SECONDPOINT)
  {
    if (numArgs!=5)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",5,numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    cur_ffi_line.pf.x = lua_tonumber(L,3);
    cur_ffi_line.pf.y = lua_tonumber(L,4);
    cur_ffi_line.pf.z = lua_tonumber(L,5);
  }
  else if (property == Property::NUMBEROFPOINTS)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    int num_points = lua_tonumber(L,3);

    if (num_points<2)
    {
      chi::log.LogAllError()
        << "Line property FFI_LINE_NUMBEROFPOINTS"
        << " used in chiFFInterpolationSetProperty. Number of points must"
        << " be greater than or equal to 2.";
      chi::Exit(EXIT_FAILURE);
    }
    cur_ffi_line.number_of_points = num_points;
  }
  else if (property == Property::CUSTOM_ARRAY)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    auto& cur_ffi_line = dcastLine(*p_ffi);

    if (not lua_istable(L, 3))
    {
      chi::log.LogAllError()
        << "Line property FFI_LINE_CUSTOM_ARRAY"
        << " used in chiFFInterpolationSetProperty. Argument 3 is expected "
           "to be an array.";
      chi::Exit(EXIT_FAILURE);
    }

    const size_t table_len = lua_rawlen(L,3);

    std::vector<double> new_array(table_len,0.0);
    for (int k=0; k<table_len; ++k)
    {
      lua_pushnumber(L,k+1);
      lua_gettable(L,3);
      new_array[k] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }

    cur_ffi_line.custom_arrays.push_back(new_array);
  }
  else if (property == Property::OPERATION)
  {
    if (numArgs!=3 and numArgs!=4)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    if (p_ffi->Type() != chi_mesh::ff_interpolation::Type::VOLUME)
      throw std::logic_error(
        "Volume property FFI_PROP_OPERATION"
        " used in chiFFInterpolationSetProperty can only be used with "
        "Volume type interpolations.");

    auto& cur_ffi_volume = dcastVolume(*p_ffi);

    int op_type = lua_tonumber(L,3);

    int OP_SUM = static_cast<int>(Operation::OP_SUM);
    int OP_MAX_LUA = static_cast<int>(Operation::OP_MAX_LUA);\
    int OP_SUM_LUA = static_cast<int>(Operation::OP_SUM_LUA);

    if (!((op_type>=OP_SUM) && (op_type<=OP_MAX_LUA)))
    {
      chi::log.LogAllError()
        << "Volume property FFI_PROP_OPERATION"
        << " used in chiFFInterpolationSetProperty. Unsupported OPERATON."
        << " Supported types are OP_AVG and OP_SUM. " << op_type;
      chi::Exit(EXIT_FAILURE);
    }

    if ((op_type >= OP_SUM_LUA) and (op_type <= OP_MAX_LUA))
    {
      if (numArgs != 4)
        LuaPostArgAmountError("chiFFInterpolationSetProperty",4,numArgs);

      const char* func_name = lua_tostring(L,4);
      cur_ffi_volume.op_lua_func = std::string(func_name);
    }

    cur_ffi_volume.op_type = static_cast<Operation>(op_type);
  }
  else if (property == Property::LOGICAL_VOLUME)
  {
    if (numArgs!=3)
      LuaPostArgAmountError("chiFFInterpolationSetProperty",3,numArgs);

    int logvol_hndle = lua_tonumber(L,3);

    auto p_logical_volume = chi::GetStackItemPtr(chi::logicvolume_stack,
                                                 logvol_hndle, fname);

    if (p_ffi->Type() != chi_mesh::ff_interpolation::Type::VOLUME)
      throw std::logic_error(
        "Volume property FFI_PROP_LOGICAL_VOLUME"
        " used in chiFFInterpolationSetProperty can only be used with "
        "Volume type interpolations.");

    auto& cur_ffi_volume = dcastVolume(*p_ffi);

    cur_ffi_volume.logical_volume = p_logical_volume;
  }
  else                                              //Fall back
  {
    chi::log.LogAllError()
      << "Invalid PropertyIndex used in chiFFInterpolationSetProperty.";
    chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
