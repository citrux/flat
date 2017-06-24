#pragma once
#include <cmath>
#include <limits>

/*
* ## Physical (CGS with eV) and mathematical constants
*/
const double e = 4.8e-10;
const double v_f = 1e8;
const double v_s = 2.6e6;
const double c = 3e10;
const double eV = 1.6e-12;
const double hbar = 6.56e-16; /* eV */
const double k = 8.6e-5;      /* eV */
const double pi = 3.14159265f;

/*
* ## Needed mathematical functions
*/
inline double sqr(double x) { return x * x; };
inline double dirac_delta(double x, double half_width) {
  if (half_width == 0)
    return (x == 0) ? std::numeric_limits<double>::infinity() : 0;

  double sigma = half_width;
  return 1.0 / sqrt(2 * pi) / sigma * exp(-sqr(x / sigma) / 2);
}
