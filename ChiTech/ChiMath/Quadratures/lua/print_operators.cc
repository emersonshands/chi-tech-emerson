#include "chi_runtime.h"
#include "ChiLua/chi_lua.h"
#include "chi_log.h"
#include "ChiMath/Quadratures/angular_quadrature_triangle.h"

int chiPrintD2M(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError(fname ,1,num_args);
  int handle = lua_tonumber(L,1);
//  std::shared_ptr<chi_math::ProductQuadrature> quad;
  std::shared_ptr<chi_math::AngularQuadratureTriangle> an_quad;
  try
  {
    auto ang_quad = chi::angular_quadrature_stack.at(handle);
    an_quad = std::static_pointer_cast<chi_math::AngularQuadratureTriangle>(ang_quad);
//    if (ang_quad->type == chi_math::AngularQuadratureType::ProductQuadrature)
//      quad = std::static_pointer_cast<chi_math::ProductQuadrature>(ang_quad);
//    else
//    {
//      an_quad = std::static_pointer_cast<chi_math::AngularQuadratureTriangle>(ang_quad);
//    }
//    else
//    {
//      chi::log.LogAllError()
//        << "chiGetProductQuadrature: Provided quadrature handle points to "
//           "a quadrature that is not a product quadrature.";
//      chi::Exit(EXIT_FAILURE);
//    }
  }
  catch (const std::out_of_range& o){
    chi::log.LogAllError()
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    chi::Exit(EXIT_FAILURE);
  }
  //CHECKING VALUES
  chi::log.Log() << "Printing weights";
  double sum = 0.0;
  double weightsum = 0.0;
  size_t sizew = an_quad->weights.size();
  for (auto& u : an_quad->weights)
  {
    weightsum += u;
    chi::log.Log() << u;
  }
  chi::log.Log() << "Weight sum " << weightsum;
  chi::log.Log() << "Integrating Zi";
  for (size_t i =0;i<sizew;++i)
  {
    double weightVal =an_quad->weights[i];
    double zival = an_quad->omegas[i].z;
    sum +=  weightVal * zival;
  }
  chi::log.Log() << sum;
  //PRINTING D2M
  an_quad->BuildDiscreteToMomentOperator();
  auto d2m = an_quad->GetDiscreteToMomentOperator();
  chi::log.Log() << "Now printing the D2M Matrix";
  for (const std::vector<double>& u : d2m)
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
  std::shared_ptr<chi_math::AngularQuadratureTriangle> an_quad;
  try
  {

    auto ang_quad = chi::angular_quadrature_stack.at(handle);
    an_quad = std::static_pointer_cast<chi_math::AngularQuadratureTriangle>(ang_quad);
//    auto& quad = chi::GetStackItem<chi_math::AngularQuadratureTriangle>(chi::angular_quadrature_stack,handle,fname);
//    quad.BuildMomentToDiscreteOperator();
    an_quad->BuildDiscreteToMomentOperator();
    auto m2d = an_quad->GetDiscreteToMomentOperator();
    chi::log.Log() << "Now printing the M2D Matrix";
    for (std::vector<double> u : m2d)
    {
      std::string vec = "";
      for (auto i: u)
        vec += std::to_string(i) + " ";
      chi::log.Log() << vec;
    }
  }
  catch (const std::out_of_range& o)
  {
    chi::log.LogAllError()
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    chi::Exit(EXIT_FAILURE);
  }

  return 0;
}