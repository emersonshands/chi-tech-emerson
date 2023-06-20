#include "RPPLogicalVolume.h"

#include "ChiObject/object_maker.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, RPPLogicalVolume);

chi_objects::InputParameters RPPLogicalVolume::GetInputParameters()
{
  chi_objects::InputParameters params = LogicalVolume::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
  "\\defgroup chi_mesh__RPPLogicalVolume RPPLogicalVolume\n"
  "\\ingroup LuaLogicVolumes");
  // clang-format on

  params.AddOptionalParameter("xmin", 0.0, "X-min of the volume");
  params.AddOptionalParameter("xmax", 1.0, "X-max of the volume");
  params.AddOptionalParameter("ymin", 0.0, "Y-min of the volume");
  params.AddOptionalParameter("ymax", 1.0, "Y-max of the volume");
  params.AddOptionalParameter("zmin", 0.0, "Z-min of the volume");
  params.AddOptionalParameter("zmax", 1.0, "Z-max of the volume");

  params.AddOptionalParameter(
    "infx", false, "Flag, when true, will ignore xmin and xmax.");
  params.AddOptionalParameter(
    "infy", false, "Flag, when true, will ignore ymin and ymax.");
  params.AddOptionalParameter(
    "infz", false, "Flag, when true, will ignore zmin and zmax.");

  return params;
}

RPPLogicalVolume::RPPLogicalVolume(const chi_objects::InputParameters& params)
  : LogicalVolume(params),
    xmin_(params.GetParamValue<double>("xmin")),
    xmax_(params.GetParamValue<double>("xmax")),
    ymin_(params.GetParamValue<double>("ymin")),
    ymax_(params.GetParamValue<double>("ymax")),
    zmin_(params.GetParamValue<double>("zmin")),
    zmax_(params.GetParamValue<double>("zmax")),
    infx_(params.GetParamValue<bool>("infx")),
    infy_(params.GetParamValue<bool>("infy")),
    infz_(params.GetParamValue<bool>("infz"))
{
}

#define XMAX 0
#define XMIN 1
#define YMAX 2
#define YMIN 3
#define ZMAX 4
#define ZMIN 5

bool RPPLogicalVolume::Inside(const chi_mesh::Vector3& point) const
{
  constexpr std::array<bool, 6> true_condition = {
    true, true, true, true, true, true};
  std::array<bool, 6> condition = {false, false, false, false, false, false};

  if (point.x <= xmax_ or infx_) condition[XMAX] = true;
  if (point.x >= xmin_ or infx_) condition[XMIN] = true;
  if (point.y <= ymax_ or infy_) condition[YMAX] = true;
  if (point.y >= ymin_ or infy_) condition[YMIN] = true;
  if (point.z <= zmax_ or infz_) condition[ZMAX] = true;
  if (point.z >= zmin_ or infz_) condition[ZMIN] = true;

  return condition == true_condition;
}

} // namespace chi_mesh