#include "chi_runtime.h"
#include "ChiLua/chi_lua.h"
#include "chi_log.h"
#include "ChiMath/Quadratures/angular_quadrature_triangle.h"

int chiGetTriangleQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetTriangleQuadrature",1,num_args);

  int handle = lua_tonumber(L,1);

  std::shared_ptr<chi_math::AngularQuadratureTriangle> quad;
  try{
    auto ang_quad = chi::angular_quadrature_stack.at(handle);
    quad = std::static_pointer_cast<chi_math::AngularQuadratureTriangle>(ang_quad);
  }
  catch (const std::out_of_range& o){
    chi::log.LogAllError()
      << "chiGetProductQuadrature: Invalid quadrature handle.";
    chi::Exit(EXIT_FAILURE);
  }

  lua_newtable(L);
  for (size_t n=0; n<quad->weights.size(); ++n)
  {
    lua_pushnumber(L,n+1);
    lua_newtable(L);

    lua_pushstring(L,"weight");
    lua_pushnumber(L,quad->weights[n]);
    lua_settable(L,-3);

    lua_pushstring(L,"polar");
    lua_pushnumber(L,quad->abscissae[n].theta);
    lua_settable(L,-3);

    lua_pushstring(L,"azimuthal");
    lua_pushnumber(L,quad->abscissae[n].phi);
    lua_settable(L,-3);

    lua_settable(L,-3);
  }

  return 1;
}

int chiPrintD2M(lua_State* L)
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
    auto d2m = an_quad->GetDiscreteToMomentOperator();
    chi::log.Log() << "Now printing the D2M Matrix";
    for (const std::vector<double>& u : d2m)
    {
      std::string vec = "";
      for (auto i: u)
        vec += std::to_string(i) + " ";
      chi::log.Log() << vec;
    }
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
  catch (std::logic_error& a){
    chi::log.LogAllError()
      << "D2M not built.";
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
    auto m2d = an_quad->GetMomentToDiscreteOperator();
    chi::log.Log() << "Now printing the M2D Matrix";
    for (std::vector<double> u : m2d)
    {
      std::string vec = "";
      for (auto i: u)
        vec += std::to_string(i) + " ";
      chi::log.Log() << vec;
    }
  }
  catch (std::logic_error& a){
    chi::log.LogAllError()
      << "M2D not built.";
    chi::Exit(EXIT_FAILURE);
  }

  return 0;
}

int chiCheckIdentity(lua_State* L)
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
    auto m2d = an_quad->GetMomentToDiscreteOperator();
    auto d2m = an_quad->GetDiscreteToMomentOperator();
    chi::log.Log() << "Now printing the Check with 1e-8 precision\n";
    auto check = chi_math::MatMul(m2d, d2m);
    size_t num_angle = an_quad->weights.size();
    for (size_t i = 0; i < num_angle; ++i)
      for (size_t j = 0; j < num_angle; ++j)
        if (fabs(check[i][j]) < 1e-8)
          check[i][j] = 0.0;
    chi_math::PrintMatrix(check);
  }
  catch (std::logic_error& a){
    chi::log.LogAllError()
      << "M2D or D2M not built.";
    chi::Exit(EXIT_FAILURE);
  }

  return 0;
}
