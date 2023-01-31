#ifndef _chi_surfacemesh_h
#define _chi_surfacemesh_h

#include<stdio.h>
#include <vector>

#include"../chi_mesh.h"

//###################################################################
/** Generic surface mesh class.
This class facilitates many functions within the mesh environment including
logically determining volumes.*/
class chi_mesh::SurfaceMesh
{
public:
  std::vector<chi_mesh::Vertex>   vertices;
  std::vector<chi_mesh::Vertex>   tex_vertices; ///< Texture vertices
  std::vector<chi_mesh::Normal>   normals;
  std::vector<chi_mesh::Face>     faces;
  std::vector<chi_mesh::Edge>     lines;
  /*
   * DO NOT PLACE COMMENTS HERE UNLESS REALLY NECESSARY
   */
  std::vector<chi_mesh::PolyFace*> poly_faces; ///<Polygonal faces

  std::vector<int> physical_region_map;

public:
  //constrdestr.cc
        SurfaceMesh();
       ~SurfaceMesh();
  friend std::ostream& operator<<(std::ostream& os,  SurfaceMesh& dt);
  //loadexport.cc
  int   ImportFromOBJFile(const std::string& fileName,bool as_poly=false,
                          const chi_mesh::Vector3& transform=Vector3(0,0,0));
  int   ImportFromTriangleFiles(const char* fileName, bool as_poly);
  int   ImportFromMshFiles(const char* fileName, bool as_poly);
  void  ExportToOBJFile(const char* fileName);
  void  ExportToPolyFile(const char* fileName);
  static SurfaceMesh*
  CreateFromDivisions(std::vector<double>& vertices_1d_x,
                      std::vector<double>& vertices_1d_y);

  //internalconn.cc
  void  UpdateInternalConnectivity();

  //checksense.cc
  bool  CheckNegativeSense(double x, double y, double z);

  //splitbypatch.cc
  void  SplitByPatch(std::vector<chi_mesh::SurfaceMesh*>& patches);

  //extractopenedges.cc
  void  ExtractOpenEdgesToObj(const char* fileName);

  //meshstats.cc
  void  CheckCyclicDependencies(int num_angles);
  void  GetMeshStats();
  void  ComputeLoadBalancing(std::vector<double>& x_cuts,
                             std::vector<double>& y_cuts);

};

#endif
