#pragma once
/*
* ## Random number generator
*
* Simple and fast xorshift for use in simulation
*
*/
struct Rng {
    unsigned int x, y, z, w;

    Rng(unsigned int seed) : x(seed),
        y(362436069),
        z(521288629),
        w(88675123) {};

    float uniform() {
        unsigned int t = x ^ (x << 1);
        x = y;
        y = z;
        z = w;
        w = w ^ (w >> 19) ^ t ^ (t >> 8);
        return ((float)w) / UINT_MAX;
    }
};