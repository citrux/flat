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
    Vec2 p;
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
};

class BigrapheneLower : public Band {
public:
    const float gamma = 0.35;
    const float delta =  0.1;

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

    BigrapheneLower(float temperature);
    ~BigrapheneLower() {};
};