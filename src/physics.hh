#pragma once
#include <limits>
#include <cmath>

/*
* ## Physical (CGS with eV) and mathematical constants
*/
const float e = 4.8e-10;
const float v_f = 1e8;
const float v_s = 2.6e6;
const float c = 3e10;
const float eV = 1.6e-12;
const float hbar = 6.56e-16; /* eV */
const float k = 8.6e-5;   /* eV */
const float pi = 3.14159265f;


/*
* ## Needed mathematical functions
*/
inline float sqr(float x) { return x * x; };
inline float dirac_delta(float x, float width) {
    if (width == 0)
        return (x == 0) ? std::numeric_limits<float>::infinity() : 0;

    float sigma = width / 2;
    return 1.0 / sqrt(2 * pi) / sigma * exp(- sqr(x / sigma) / 2);
}
