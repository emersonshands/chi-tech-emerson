#include "diffusion_solver.h"

#include <chi_log.h>
extern ChiLog& chi_log;

/**\defgroup LuaDiffusionBasicOptions Basic Options
 * \ingroup LuaDiffusion
 *
Option name           | Type   | Default Value | Description
----------------------|------- |---------------|------------
discretization_method | string | 500           | Spatial discretization method.
max_iters             | int    | 500           | Maximum iterations for solver convergence.
residual_tolerance    | float  | 1.0e-8        | Residual convergence tolerance.
property_map_D        | int    | 0             | Material property index to use for diffusion coefficient
property_map_q        | int    | 1             | Material property index to use for source.
property_map_sigma    | int    | 2             | Material property index to use for interaction coefficient.
 * Very nice stuff.*/

chi_diffusion::Solver::Solver(const std::string& in_solver_name):
  chi_physics::Solver(in_solver_name, {{"discretization_method", std::string()},
                                       {"max_iters", int64_t(500)},
                                       {"residual_tolerance", 1.0e-8},
                                       {"property_map_D",int64_t(0)},
                                       {"property_map_q",int64_t(1)},
                                       {"property_map_sigma",int64_t(2)}})
{}

chi_diffusion::Solver::~Solver()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Cleaning up diffusion solver: " << TextName();

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
  KSPDestroy(&ksp);

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Done cleaning up diffusion solver: " << TextName();
}