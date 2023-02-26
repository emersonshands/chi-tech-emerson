#include "chi_unpartitioned_mesh.h"

#include "ChiMesh/MeshContinuum/chi_grid_vtk_utils.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridReader.h>

#include <vtkInformation.h>

#include "chi_runtime.h"
#include "chi_log.h"

#define ErrorReadingFile(fname) \
std::runtime_error("Failed to open file: " + options.file_name + \
" in call to " + #fname + ".")

//###################################################################
/**Reads a VTK unstructured mesh. This reader will use the following
 * options:
 * - `file_name`, of course.
 * - `material_id_fieldname`, cell data for material_id.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromPVTU(const chi_mesh::UnpartitionedMesh::Options &options)
{
  chi::log.Log() << "Reading PVTU file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromVTU);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  //======================================== Separate the blocks
  // For vtu files this is very simple. The
  // output of the reader is an UnstructuredGrid.
  auto ugrid_main = vtkUGridPtr(reader->GetOutput());
  std::vector<vtkUGridPtrAndName> grid_blocks = {{ugrid_main,""}};

  //======================================== Process blocks
  const int max_dimension = chi_mesh::FindHighestDimension(grid_blocks);
  auto ugrid = chi_mesh::ConsolidateAndCleanBlocks(grid_blocks, max_dimension);

  //======================================== Copy Data
  CopyUGridCellsAndPoints(*ugrid, options.scale);

  //======================================== Set material ids
  const auto material_ids = chi_mesh::BuildCellMaterialIDsFromField(
    ugrid, options.material_id_fieldname, options.file_name);
  SetMaterialIDsFromList(material_ids);

  //======================================== Always do this
  chi_mesh::MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1: dimension = DIMENSION_1; break;
    case 2: dimension = DIMENSION_2; break;
    case 3: dimension = DIMENSION_3; break;
    default: break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  chi::log.Log() << "Done reading PVTU file: " << options.file_name << ".";
}

