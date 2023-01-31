#ifndef CHI_PHYSICS_NAMESPACE_H
#define CHI_PHYSICS_NAMESPACE_H
/**\defgroup LuaPhysics C Physics*/

#include <petscksp.h>

namespace chi_physics
{
  enum class OperationType
  {
    SINGLE_VALUE = 0,
    FROM_ARRAY   = 1,
    SIMPLEXS0    = 20,
    SIMPLEXS1    = 21,
    EXISTING     = 22,
    CHI_XSFILE   = 23
  };

  class FieldFunction;
  class FieldFunction;
  class Solver;

  //03 Utils
  std::string GetPETScConvergedReasonstring(KSPConvergedReason reason);
}


#endif