#include "pwl_polyhedron.h"

/**Precomputes the shape function values of a face-side pair
 * at a quadrature point*/
double chi_math::PolyhedronMappingFE_PWL::FaceSideShape(unsigned int face_index,
                                                        unsigned int side_index,
                                                        unsigned int i,
                                                        const chi_mesh::Vector3& qpoint,
                                                        bool on_surface/*=false*/) const
{
  double value = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  value += TetShape(index, qpoint, on_surface);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    value += betaf* TetShape(1, qpoint, on_surface);}
  value += alphac* TetShape(3, qpoint, on_surface);

  return value;
}

/**Precomputes the gradx-shape function values of a face-side pair
 * at a quadrature point*/
double chi_math::PolyhedronMappingFE_PWL::FaceSideGradShape_x(unsigned int face_index,
                                                              unsigned int side_index,
                                                              unsigned int i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 2) * tetdfdz;

  return value;
}

/**Precomputes the grady-shape function values of a face-side pair
 * at a quadrature point*/
double chi_math::PolyhedronMappingFE_PWL::FaceSideGradShape_y(unsigned int face_index,
                                                              unsigned int side_index,
                                                              unsigned int i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 2) * tetdfdz;

  return value;
}

/**Precomputes the gradz-shape function values of a face-side pair
 * at a quadrature point*/
double chi_math::PolyhedronMappingFE_PWL::FaceSideGradShape_z(unsigned int face_index,
                                                              unsigned int side_index,
                                                              unsigned int i) const
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 2) * tetdfdz;

  return value;
}