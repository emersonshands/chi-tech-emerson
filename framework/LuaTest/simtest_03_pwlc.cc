#include "chi_lua.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

namespace chi_unit_sim_tests
{

/**This is a simple test of the Finite Volume spatial discretization applied
 * to Laplace's problem. */
int chiSimTest03_PWLC(lua_State* L)
{
    chi::log.Log() << "Coding Tutorial 3";

    //============================================= Get grid
    auto grid_ptr = chi_mesh::GetCurrentHandler().GetGrid();
    const auto& grid = *grid_ptr;

    chi::log.Log() << "Global num cells: " << grid.GetGlobalNumberOfCells();

    //============================================= Make SDM
    typedef std::shared_ptr<chi_math::SpatialDiscretization> SDMPtr;
    SDMPtr sdm_ptr = chi_math::SpatialDiscretization_PWLC::New(grid_ptr);
    const auto& sdm = *sdm_ptr;

    const auto& OneDofPerNode = sdm.UNITARY_UNKNOWN_MANAGER;

    const size_t num_local_dofs = sdm.GetNumLocalDOFs(OneDofPerNode);
    const size_t num_globl_dofs = sdm.GetNumGlobalDOFs(OneDofPerNode);

    chi::log.Log() << "Num local DOFs: " << num_local_dofs;
    chi::log.Log() << "Num globl DOFs: " << num_globl_dofs;

    //============================================= Initializes Mats and Vecs
    const auto n = static_cast<int64_t>(num_local_dofs);
    const auto N = static_cast<int64_t>(num_globl_dofs);
    Mat A;
    Vec x,b;

    A = chi_math::PETScUtils::CreateSquareMatrix(n,N);
    x = chi_math::PETScUtils::CreateVector(n,N);
    b = chi_math::PETScUtils::CreateVector(n,N);

    std::vector<int64_t> nodal_nnz_in_diag;
    std::vector<int64_t> nodal_nnz_off_diag;
    sdm.BuildSparsityPattern(nodal_nnz_in_diag,nodal_nnz_off_diag, OneDofPerNode);

    chi_math::PETScUtils::InitMatrixSparsity(A,
                                             nodal_nnz_in_diag,
                                             nodal_nnz_off_diag);

    //============================================= Assemble the system
    chi::log.Log() << "Assembling system: ";
    for (const auto& cell : grid.local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();

      const size_t num_nodes = cell_mapping.NumNodes();
      MatDbl Acell(num_nodes, VecDbl(num_nodes, 0.0));
      VecDbl cell_rhs(num_nodes, 0.0);

      for (size_t i=0; i<num_nodes; ++i)
      {
        for (size_t j=0; j<num_nodes; ++j)
        {
          double entry_aij = 0.0;
          for (size_t qp : qp_data.QuadraturePointIndices())
          {
            entry_aij +=
              qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp)) *
              qp_data.JxW(qp);
          }//for qp
          Acell[i][j] = entry_aij;
        }//for j
        for (size_t qp : qp_data.QuadraturePointIndices())
          cell_rhs[i] += 1.0 * qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      }//for i

      //======================= Flag nodes for being on dirichlet boundary
      std::vector<bool> node_boundary_flag(num_nodes, false);
      const size_t num_faces = cell.faces.size();
      for (size_t f=0; f<num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        if (face.has_neighbor) continue;

        const size_t num_face_nodes = face.vertex_ids.size();
        for (size_t fi=0; fi<num_face_nodes; ++fi)
        {
          const uint i = cell_mapping.MapFaceNode(f,fi);
          node_boundary_flag[i] = true;
        }//for fi
      }//for face f

      //======================= Develop node mapping
      std::vector<int64_t> imap(num_nodes, 0); //node-mapping
      for (size_t i=0; i<num_nodes; ++i)
        imap[i] = sdm.MapDOF(cell, i);

      //======================= Assembly into system
      for (size_t i=0; i<num_nodes; ++i)
      {
        if (node_boundary_flag[i]) //if dirichlet node
        {
          MatSetValue(A, imap[i], imap[i], 1.0, ADD_VALUES);
          VecSetValue(b, imap[i], 0.0, ADD_VALUES);
        }
        else
        {
          for (size_t j=0; j<num_nodes; ++j)
          {
            if (not node_boundary_flag[j])
              MatSetValue(A, imap[i], imap[j], Acell[i][j], ADD_VALUES);
          }//for j
          VecSetValue(b, imap[i], cell_rhs[i], ADD_VALUES);
        }
      }//for i
    }//for cell

    chi::log.Log() << "Global assembly";

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    chi::log.Log() << "Done global assembly";

    //============================================= Create Krylov Solver
    chi::log.Log() << "Solving: ";
    auto petsc_solver =
      chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
        A,               //Matrix
        "PWLCDiffSolver",  //Solver name
        KSPCG,           //Solver type
        PCGAMG,          //Preconditioner type
        1.0e-6,          //Relative residual tolerance
        1000);            //Max iterations

    //============================================= Solve
    KSPSolve(petsc_solver.ksp,b,x);

    chi::log.Log() << "Done solving";

    //============================================= Extract PETSc vector
    std::vector<double> field;
    sdm.LocalizePETScVector(x,field,OneDofPerNode);

    //============================================= Clean up
    KSPDestroy(&petsc_solver.ksp);

    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);

    chi::log.Log() << "Done cleanup";

    //============================================= Create Field Function
    auto ff = std::make_shared<chi_physics::FieldFunction>(
      "Phi",                                           //Text name
      sdm_ptr,                                         //Spatial Discr.
      chi_math::Unknown(chi_math::UnknownType::SCALAR) //Unknown
    );

    ff->UpdateFieldVector(field);

    ff->ExportToVTK("CodeTut3_PWLC");

  return 0;
}

}//namespace chi_unit_tests
