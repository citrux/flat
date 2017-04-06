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

struct ScatteringResult {
    float momentum;
    float rate;
};

class Band {
public:
    virtual float min_energy() const = 0;
    virtual float energy(float momentum) const = 0;
    virtual float energy(Vec2 const & momentum) const = 0;
    virtual float velocity(float momentum) const = 0;
    virtual Vec2 velocity(Vec2 const & momentum) const = 0;
    virtual std::list<ScatteringResult> acoustic_phonon_scattering(Particle & p) = 0;
    virtual std::list<ScatteringResult> optical_phonon_scattering(Particle & p) = 0;
    Vec2 momentum_scattering(float momentum, Particle & p);
};

namespace Bigraphene {
class Lower : public Band {
public:
    const float gamma = 0.35;
    const float delta = 0;

    const float optical_phonon_energy = 0.196;
    float acoustic_phonon_constant;
    float optical_phonon_constant;

    std::vector<BandScatteringEntry> table;
    int energy_samples;
    int momentum_samples;
    float momentum_precision;

    float min_energy() const;
    float energy(float momentum) const;
    float energy(Vec2 const & momentum) const;
    float velocity(float momentum) const;
    Vec2 velocity(Vec2 const & momentum) const;
    std::list<ScatteringResult> acoustic_phonon_scattering(Particle & p);
    std::list<ScatteringResult> optical_phonon_scattering(Particle & p);

    Lower(float temperature);
    ~Lower() {};
};
class Upper : public Band {
public:
    const float gamma = 0.35;
    const float delta = 0;

    const float optical_phonon_energy = 0.196;
    float acoustic_phonon_constant;
    float optical_phonon_constant;

    std::vector<BandScatteringEntry> table;
    int energy_samples;
    int momentum_samples;
    float momentum_precision;

    float min_energy() const;
    float energy(float momentum) const;
    float energy(Vec2 const & momentum) const;
    float velocity(float momentum) const;
    Vec2 velocity(Vec2 const & momentum) const;
    std::list<ScatteringResult> acoustic_phonon_scattering(Particle & p);
    std::list<ScatteringResult> optical_phonon_scattering(Particle & p);

    Upper(float temperature);
    ~Upper() {};
};
};

namespace Graphene {
class Bnd : public Band {
public:
    const float delta = 0.13;

    const float optical_phonon_energy = 0.196;
    float acoustic_phonon_constant;
    float optical_phonon_constant;

    float min_energy() const;
    float energy(float momentum) const;
    float energy(Vec2 const & momentum) const;
    float velocity(float momentum) const;
    Vec2 velocity(Vec2 const & momentum) const;
    std::list<ScatteringResult> acoustic_phonon_scattering(Particle & p);
    std::list<ScatteringResult> optical_phonon_scattering(Particle & p);

    Bnd(float temperature);
    ~Bnd() {};
};
}
