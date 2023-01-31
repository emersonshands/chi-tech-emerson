#include "mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Assemble_A_bext()
{
  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  //============================================= Assemble the system
  unsigned int i_two_grid = do_two_grid ? 1 : 0;

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();
    const size_t num_nodes   = cell_mapping.NumNodes();

    const auto& xs   = matid_to_xs_map.at(cell.material_id);
    const auto& qext = matid_to_src_map.at(cell.material_id);
    double collapsed_D=0.,collapsed_sig_a=0.;
    if (do_two_grid)
    {
      const auto &xstg = map_mat_id_2_tginfo.at(cell.material_id);
      collapsed_D = xstg.collapsed_D;
      collapsed_sig_a = xstg.collapsed_sig_a;
    }

    std::vector<VecDbl> rhs_cell;
    rhs_cell.resize(num_groups);
    for (uint g=0; g < num_groups; ++g)
      rhs_cell[g].resize(num_nodes, 0.0);

    std::vector<MatDbl> Acell;
    Acell.resize(num_groups+i_two_grid);
    for (uint g=0; g < num_groups+i_two_grid; ++g)
      Acell[g].resize(num_nodes, VecDbl(num_nodes, 0.0));

    for (size_t i=0; i<num_nodes; ++i)
    {
      for (size_t j=0; j<num_nodes; ++j)
      {
        double entry_mij = 0.0;
        double entry_kij = 0.0;
        for (size_t qp : qp_data.QuadraturePointIndices())
        {
          entry_mij +=  qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp)
                        * qp_data.JxW(qp);
          entry_kij +=  qp_data.ShapeGrad(i, qp).Dot(qp_data.ShapeGrad(j, qp))
                        * qp_data.JxW(qp);
        }//for qp
        for (uint g=0; g < num_groups; ++g)
        {
          const double Dg     = xs->diffusion_coeff[g];
          const double sigr_g = xs->sigma_removal[g];
          Acell[g][i][j] = entry_mij * sigr_g + entry_kij * Dg;
        }
        if (do_two_grid)
          Acell[num_groups][i][j] = entry_mij * collapsed_sig_a
                                  + entry_kij * collapsed_D;
      }//for j
      double entry_rhsi = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        entry_rhsi += qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      for (uint g=0; g < num_groups; ++g)
        rhs_cell[g][i] = entry_rhsi * (qext->source_value_g[g]);
    }//for i

    //======================= Deal with BC (all based on variations of Robin)
    const size_t num_faces = cell.faces.size();
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      // not a boundary face
	    if (face.has_neighbor) continue; 
	  
      auto& bndry = boundaries[face.neighbor_id];

      // Robin boundary
      //   for two-grid, it is homogenous Robin
      if (bndry.type == BoundaryType::Robin)
      {
        const auto  qp_face_data = cell_mapping.MakeFaceQuadraturePointData( f );
        const size_t num_face_nodes = face.vertex_ids.size();

        auto& aval = bndry.mg_values[0];
        auto& bval = bndry.mg_values[1];
        auto& fval = bndry.mg_values[2];
        if (do_two_grid)
        {
          aval.push_back(0.25);
          bval.push_back(0.5);
          fval.push_back(0.0);
        }


        // sanity check, Assert if b=0
        for (uint g=0; g < num_groups+i_two_grid; ++g)
        {
          if (std::fabs(bval[g]) < 1e-8)
            throw std::logic_error("if b=0, this is a Dirichlet BC, not a Robin BC");
        }

        // true Robin when a!=0, otherwise, it is a Neumann:
        // only do this part if true Robin (i.e., a!=0)
        for (uint g=0; g < num_groups+i_two_grid; ++g)
        {
          if (std::fabs(aval[g]) > 1.0e-8)
          {
            // loop over nodes of that face
            for (size_t fi=0; fi<num_face_nodes; ++fi)
            {
              const uint i = cell_mapping.MapFaceNode(f,fi);

              double entry_rhsi = 0.0;
              for (size_t qp : qp_face_data.QuadraturePointIndices() )
                entry_rhsi +=  qp_face_data.ShapeValue(i, qp) * qp_face_data.JxW(qp);
              if (g<num_groups) // check due to two-grid
                rhs_cell[g][i] +=  fval[g] / bval[g] * entry_rhsi;

              for (size_t fj=0; fj<num_face_nodes; ++fj)
              {
                const uint j = cell_mapping.MapFaceNode(f,fj);
                double entry_aij = 0.0;
                for (size_t qp : qp_face_data.QuadraturePointIndices())
                  entry_aij +=  qp_face_data.ShapeValue(i, qp)
                                * qp_face_data.ShapeValue(j, qp)
                                * qp_face_data.JxW(qp);
                Acell[g][i][j] += aval[g] / bval[g] * entry_aij;
              }//for fj
            }//for fi
          }//end true Robin
        }//for g
      }//if Robin
    }//for face f
 
    //======================= Develop node mapping
    std::vector<int64_t> imap(num_nodes, 0); //node-mapping
    for (size_t i=0; i<num_nodes; ++i)
      imap[i] = sdm.MapDOF(cell, i);

    //======================= Assembly into system
    for (uint g=0; g < num_groups; ++g)
      for (size_t i=0; i<num_nodes; ++i)
        VecSetValue(bext[g], imap[i], rhs_cell[g][i], ADD_VALUES);

    for (uint g=0; g < num_groups+i_two_grid; ++g)
      for (size_t i=0; i<num_nodes; ++i)
        for (size_t j=0; j<num_nodes; ++j)
          MatSetValue(A[g], imap[i], imap[j], Acell[g][i][j], ADD_VALUES);

  }//for cell
 
  chi::log.Log() << "Global assembly";

  for (uint g=0; g < num_groups; ++g)
  {
    VecAssemblyBegin(bext[g]);
    VecAssemblyEnd(bext[g]);
  }
  for (uint g=0; g < num_groups+i_two_grid; ++g)
  {
    MatAssemblyBegin(A[g], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A[g], MAT_FINAL_ASSEMBLY);
  }

//  PetscViewer viewer;
//  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A2_before_bc.m",&viewer);
//  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//  MatView(A[0],viewer);
//  PetscViewerPopFormat(viewer);
//  PetscViewerDestroy(&viewer);
//  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"bext2_before_bc.m",&viewer);
//  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
//  VecView(bext[0],viewer);
//  PetscViewerPopFormat(viewer);
//  PetscViewerDestroy(&viewer);

  chi::log.Log() << "Done global assembly";

}