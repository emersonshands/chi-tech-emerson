#include "lbs_adjoint.h"

#include "ChiMath/chi_math.h"
#include "../ChiMath/SerialNewtonIteration/serial_newton_iteration.h"

#include "chi_runtime.h"
#include "chi_runtime.h"
#include "chi_log.h"

void lbs_adjoint::TestFunction()
{
  std::cout << "Test Function!\n";
}

std::array<double,2> lbs_adjoint::
  MakeExpRepFromP1(const std::array<double, 4> &P1_moments,
                   bool verbose/*=false*/)
{
  //=============================================
  // Custom function to implement the non-linear equations
  // that make up the system to solve the a and b coefficients
  // of an exponential representation of the form exp(a + b*mu).
  // The coefficients must be such that the following equations
  // hold:
  // 2\pi \int_{-1}^{+1} e^{a + b\mu} d\mu = \phi
  // 2\pi \int_{-1}^{+1} e^{a + b\mu} \mu d\mu = ||J||_2
  //
  // These integrals are analytical in \mu and therefore results
  // in the equations for the F-method. The Jacobian matrix is returned
  // with the J-method.
  //
  // Note that the function only accept the current normalized with
  // 1/phi. So phi is always 1.0.
  class CustomF : public chi_math::NonLinearFunction
  {
  private:
    double J_x, J_y, J_z;

  public:
    /**Constructor. Parameters are the 3 components of the currents
     * normalized with 1/phi.*/
    explicit CustomF(const std::array<double,3>& input) :
    J_x(input[0]), J_y(input[1]), J_z(input[2]) {}

    /**Function evaluation at vector-x.*/
    VecDbl F(const VecDbl& x) const override
    {
      assert(x.size() == 2);
      const double a = x[0];
      const double b = x[1];
      const double FOUR_PI = 4.0*M_PI;

      double size_J = chi_math::Vec2Norm({J_x, J_y, J_z});

      return {(FOUR_PI/b)   * exp(a) * sinh(b) - 1.0,
              (FOUR_PI/b/b) * exp(a) * (b * cosh(b) - sinh(b)) - size_J};

    }
    /**Jacobian evaluation at vector-x.*/
    MatDbl J(const VecDbl& x) const override
    {
      assert(x.size() == 2);
      const double a = x[0];
      const double b = x[1];
      const double FOUR_PI = 4.0*M_PI;

      return {{(FOUR_PI/b)     * exp(a) * sinh(b),
               (FOUR_PI/b/b)   * exp(a) * (b * cosh(b) - sinh(b))},
              {(FOUR_PI/b/b)   * exp(a) * (b * cosh(b) - sinh(b)),
               (FOUR_PI/b/b/b) * exp(a) * ((b*b + 2) * sinh(b) - 2 * b * cosh(b))}};
    }
  };

  //======================================== Convert P1 moments to phi and J
  double phi      = P1_moments[0];
  double J_x      = P1_moments[1];
  double J_y      = P1_moments[2];
  double J_z      = P1_moments[3];

  //======================================== Compute initial ratio size_J/phi
  double size_J_i = chi_math::Vec2Norm({J_x, J_y, J_z});
  double ratio_i  = size_J_i/phi;

  if (phi < 1.0e-10 or ratio_i > 0.9)
  { J_x *= 0.9 / size_J_i; J_y *= 0.9 / size_J_i; J_z *= 0.9 / size_J_i; }
  else
  { J_x /= phi; J_y /= phi; J_z /= phi; }

  double size_J_f = chi_math::Vec2Norm({J_x, J_y, J_z});
  double ratio_f  = size_J_f;

  if (verbose)
  {
    std::stringstream outstr;
    outstr << "P1 moments initial: " <<P1_moments[0]<<" "
                                     <<P1_moments[1]<<" "
                                     <<P1_moments[2]<<" "
                                     <<P1_moments[3]<<" |J|=";
    outstr << size_J_i << " ratio=" << ratio_i << "\n";
    outstr << "P1 moments initial: " <<1.0<<" "<<J_x<<" "<<J_y<<" "<<J_z<<" |J|=";
    outstr << size_J_f << " ratio=" << ratio_f << "\n";

    chi::log.Log() << outstr.str();
  }

  if (size_J_f < 1.0e-10)
  {
    double a = log(phi/4.0/M_PI);
    double b = 0.0;

    if (verbose) { chi::log.Log() << "Solution: " << a << " " << b; }

    return {a, b};
  }
  else
  {
    CustomF custom_function({J_x, J_y, J_z});
    auto solution = chi_math::NewtonIteration(
      custom_function, //Non-linear function
      {1.0,0.1},       //Initial guess
      100,             //Max iterations
      1.0e-8,          //Tolerance
      verbose);        //Verbose output?

    double a = solution[0];
    double b = solution[1];

    if (verbose) { chi::log.Log() << "Solution: " << a << " " << b; }

    return {a, b};
  }
}
