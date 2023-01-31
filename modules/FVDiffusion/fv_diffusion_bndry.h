#ifndef FV_DIFFUSION_BOUNDARY_H
#define FV_DIFFUSION_BOUNDARY_H

#include <array>

namespace fv_diffusion
{
  class Boundary;
  
  enum class BoundaryType : int
  {
    Reflecting = 1,
    Dirichlet  = 2,
    Neumann    = 3,
    Robin      = 4,
    Vacuum     = 5
  };
}

//###################################################################
/**Parent class for diffusion boundaries*/
class fv_diffusion::Boundary
{
  public :
  BoundaryType type = BoundaryType::Dirichlet;

  std::array<double, 3> values = {0.,0.,0.};
};

#endif //FV_DIFFUSION_BOUNDARY_H