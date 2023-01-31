#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "ChiMath/Quadratures/angular_quadrature_triangle.h"

#include <memory>

int chiCreateAngularQuadratureTriangle(lua_State* L){
  const std::string  fname = __FUNCTION__;
  //we should be given a sn and method
  const int num_args = lua_gettop(L);
  if (num_args<2)
    LuaPostArgAmountError(fname,2,num_args);
  const int method = lua_tointeger(L,1);
  const int order = lua_tointeger(L,2);
  int moments = 0;
  if (num_args==3)
    moments = lua_tointeger(L,3);
  //constructor should make the matrix automatically given the order and method
  auto new_quad = std::make_shared<chi_math::AngularQuadratureTriangle>
    (method, order, moments);

  //push onto stack
  chi::angular_quadrature_stack.push_back(new_quad);

  //push back onto the stack and return location
  const size_t index = chi::angular_quadrature_stack.size() - 1;
  lua_pushnumber(L,static_cast<lua_Number>(index));

  return 1;
}