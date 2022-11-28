#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "ChiMath/Quadratures/angular_product_quadrature.h"

#include "chi_log.h"

//########################################################## Get product quadrature
/** Get the values of a product quadrature

\param QuadHandle int Handle to an existing product quadrature.

\return Table A lua table with each entry being another table with entries
       .weight .polar .azimuthal.

\ingroup LuaQuadrature
\author Jan*/
int chiGetProductQuadrature(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetProductQuadrature",1,num_args);

  int handle = lua_tonumber(L,1);

  std::shared_ptr<chi_math::ProductQuadrature> quad;
  try{
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