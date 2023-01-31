#ifndef _surfmesher_predefined_h
#define _surfmesher_predefined_h

#include "../surfacemesher.h"

//######################################################### Class def
/**Surface mesher that will not modify the mesh.
Meant for loading 2D meshes and just connecting boundaries
to elements.*/
class chi_mesh::SurfaceMesherPredefined : public chi_mesh::SurfaceMesher
{
public:
  SurfaceMesherPredefined() : SurfaceMesher(SurfaceMesherType::Predefined)
  {}
  //02 Execute
  void Execute() override;
};

#endif