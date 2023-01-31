#ifndef _chi_math_rng_h
#define _chi_math_rng_h

#include <random>

namespace chi_math
{
//#########################################################
/**Random number generator based on threefry.*/
class RandomNumberGenerator
{
private:
  std::mt19937_64    mt1993764_generator;
  std::uniform_real_distribution<double> distribution;

public:
  RandomNumberGenerator();
  RandomNumberGenerator(int seed);
  double Rand();
};
}

#endif
