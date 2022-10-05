#include "chi_runtime.h"
#include "ChiLua/chi_lua.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "chi_log.h"

int chiPrintD2M(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError(fname ,1,num_args);
  int handle = lua_tonumber(L,1);
  auto& quad = chi::GetStackItem<chi_math::AngularQuadrature>(chi::angular_quadrature_stack,handle,fname);
  auto op = quad.GetDiscreteToMomentOperator();
  chi::log.Log() << "Now printing the D2M Matrix";
  for (std::vector<double> u : op)
  {
    std::string vec = "";
    for (auto i: u)
      vec += std::to_string(i) + " ";
    chi::log.Log() << vec;
  }
  return 0;
}
int chiPrintM2D(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError(fname ,1,num_args);
  int handle = lua_tonumber(L,1);
  std::cout << handle << std::endl;
  auto& quad = chi::GetStackItem<chi_math::AngularQuadrature>(chi::angular_quadrature_stack,handle,fname);
  auto op = quad.GetMomentToDiscreteOperator();
  chi::log.Log() << "Now printing the M2D Matrix";
  for (std::vector<double> u : op)
  {
    std::string vec = "";
    for (auto i: u)
      vec += std::to_string(i) + " ";
    chi::log.Log() << vec;
  }
  return 0;
}