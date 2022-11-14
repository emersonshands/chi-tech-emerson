#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "ChiMath/Quadratures/product_quadrature_op.h"

#include <memory>

int chiCreateProductQuadratureOperator(lua_State *L)
{
  const std::string  fname = __FUNCTION__;
  //we should be given a product quadrature, sn, and method
  const int num_args = lua_gettop(L);
  if (num_args<3)
    LuaPostArgAmountError(fname,3,num_args);
  //Given the handle from the stack
  int handle = lua_tonumber(L,1);
  const int order = lua_tointeger(L,2);
  const int method = lua_tointeger(L,3);
  //find the pointer from the stack given the handle
  auto& quad = chi::GetStackItem<chi_math::ProductQuadrature>(chi::angular_quadrature_stack,handle,fname);

  //Run the actual method to form the new quadrature.
  auto new_quad = std::make_shared<chi_math::ProductQuadratureOp>
    (quad, order,method);


  return 1;
}