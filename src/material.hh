#pragma once

#include "vec2.hh"
#include "rng.hh"
#include <list>
#include <vector>
#include <string>

struct Band;

struct Particle {
    Vec2 p;
    float r;
    Rng rng;
    Band *band;
    void reset_r();
    void isotropic_scattering(float momentum);
    Particle(unsigned int seed);
};

struct BandScatteringIntegral {
    float momentum;
    float integral;
};

struct BandScatteringEntry {
    float energy;
    std::list<BandScatteringIntegral> integrals;
};

struct Band {
    const float gamma = 0.35;
    const float delta = 0.1;

    float acoustic_phonon_constant;
    float optical_phonon_constant;
    float optical_phonon_energy;

    std::vector<BandScatteringEntry> table;
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
    ~Band() {};
};