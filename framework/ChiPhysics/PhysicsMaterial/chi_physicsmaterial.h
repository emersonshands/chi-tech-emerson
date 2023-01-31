#ifndef CHI_PHYSICS_MATERIAL_H
#define CHI_PHYSICS_MATERIAL_H

#include "ChiPhysics/chi_physics_namespace.h"
#include "material_property_base.h"

#include <vector>
#include <memory>

namespace chi_physics
{

//###################################################################
/** Base class for materials used in physics simulations.*/
class Material
{
public:
  std::vector<std::shared_ptr<MaterialProperty>> properties{};
  std::string name="Unnamed Material";

};

}//namespace chi_physics

#endif