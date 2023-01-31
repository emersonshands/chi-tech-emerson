#ifndef CHI_PHYSICS_MATERIAL_PROPERTY_BASE_H
#define CHI_PHYSICS_MATERIAL_PROPERTY_BASE_H

#include <string>
#include <vector>

#include "ChiLua/chi_lua.h"

namespace chi_physics
{
enum class PropertyType
{
  SCALAR_VALUE = 1,
  TRANSPORT_XSECTIONS = 10,
  ISOTROPIC_MG_SOURCE = 11
};

//###################################################################
/** Base class for material properties.*/
class MaterialProperty
{
private:
  const PropertyType type;
public:
  std::string property_name;

  explicit MaterialProperty(PropertyType in_type) : type(in_type) {}

  virtual ~MaterialProperty() = default;

  PropertyType Type() { return type; }

  virtual double GetScalarValue() { return 0.0; }

  virtual void PushLuaTable(lua_State *L);
};

}

#endif //CHI_PHYSICS_MATERIAL_PROPERTY_H