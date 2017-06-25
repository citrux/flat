#pragma once

#include <climits>
/*
* ## Random number generator
*
* Simple and fast xorshift for use in simulation
*
*/
struct Rng {
  unsigned int x, y, z, w;

  Rng(unsigned int seed) : x(seed), y(362436069), z(521288629), w(88675123){};
  Rng() : x(0), y(362436069), z(521288629), w(88675123){};
  double uniform() {
    unsigned int t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = w ^ (w >> 19) ^ t ^ (t >> 8);
    return ((double)w) / UINT_MAX;
  }
};
