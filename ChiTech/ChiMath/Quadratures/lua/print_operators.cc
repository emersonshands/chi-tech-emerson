#include "chi_runtime.h"
#include "ChiLua/chi_lua.h"
#include "ChiMath/Quadratures/angular_quadrature_base.h"
#include "ChiMath/Quadratures/angular_product_quadrature.h"
#include "chi_log.h"

int chiPrintD2M(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError(fname ,1,num_args);
  int handle = lua_tonumber(L,1);
  std::shared_ptr<chi_math::ProductQuadrature> quad;
  try
  {
    auto ang_quad = chi::angular_quadrature_stack.at(handle);
    if (ang_quad->type == chi_math::AngularQuadratureType::ProductQuadrature)
      quad = std::static_pointer_cast<chi_math::ProductQuadrature>(ang_quad);
    else
    {
      chi::log.LogAllError()
        << "chiGetProductQuadrature: Provided quadrature handle points to "
           "a quadrature that is not a product quadrature.";
      chi::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o){
    chi::log.LogAllError()
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    chi::Exit(EXIT_FAILURE);
  }
  auto op = quad->GetDiscreteToMomentOperator();
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