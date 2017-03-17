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
    const float delta =  0.1;

    const float optical_phonon_energy = 0.196;
    float acoustic_phonon_constant;
    float optical_phonon_constant;
    

    std::vector<BandScatteringEntry> table;
    int energy_samples;
    int momentum_samples;
    float momentum_precision;

    inline float min_energy() const {
        return delta * gamma / std::sqrt(gamma * gamma + 4 * delta * delta);
        }

    inline float energy(float momentum) const {
        const float delta2 = delta * delta;
        const float gamma2 = gamma * gamma;
        const float gamma4 = gamma2 * gamma2;
        const float p2 = momentum * momentum;
        return std::sqrt(delta2 + gamma2 / 2 + p2 - std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2));
    }

    inline float energy(Vec2 const & momentum) const {
        return energy(momentum.len());
    }

    inline float velocity(float momentum) const {
        return velocity(Vec2(momentum, 0)).x;
    }

    inline Vec2 velocity(Vec2 const & momentum) const {
        const float delta2 = delta * delta;
        const float gamma2 = gamma * gamma;
        const float gamma4 = gamma2 * gamma2;
        const float p2 = momentum.dot(momentum);
        return momentum * (1 - 0.5 * (gamma2 + 4 * delta2) / std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) / energy(momentum);
    }


    bool acoustic_phonon_scattering(Particle & p, float dt);
    bool optical_phonon_scattering(Particle & p, float dt);

    Band(float temperature);
    ~Band() {};
};