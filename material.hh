#pragma once

#include "vec2.hh"
#include "rng.hh"
#include <list>
#include <string>

struct Particle {
    Vec2 p;
    float r;
    Rng rng;
    void reset_r();
    void anisotropic_scattering(float momentum);
    Particle(unsigned int seed);
}

struct BandScatteringIntegral {
    float energy;
    float momentum;
    float integral;
}

struct Band {
    float gamma;
    float delta;

    float acoustic_phonon_constant;
    float optical_phonon_constant;
    float optical_phonon_energy;

    list<BandScatteringEntry> *table;
    int energy_samples;
    int momentum_samples;
    float momentum_precision;

    float energy(float momentum);
    float energy(Vec2 momentum);
    float velocity(float momentum);
    Vec2 velocity(Vec2 momentum);

    bool acoustic_phonon_scattering(Particle & p);
    bool optical_phonon_scattering(Particle & p);

    Band();
    ~Band();
}