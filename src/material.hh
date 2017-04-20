#pragma once

#include "vec2.hh"
#include "rng.hh"
#include <list>
#include <vector>
#include <string>

struct Wave {
    Vec2 E;
    float H;
    float omega;
    float phi;
};

class Band;

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

class Material {
public:
    std::vector<Band*> bands;
    //virtual float vertical_transition(int from, int to, Wave const & wave) = 0;
};

class Bigraphene : public Material {
public:
    Bigraphene(float temperature, float delta, float number_of_bands);
};

class Graphene : public Material {
public:
    Graphene(float temperature, float delta);
};

template <typename T>
inline float bisection(T f, float xmin, float xmax, float precision) {
    while (xmax - xmin > precision) {
        float x = (xmax + xmin) / 2;
        if (f(x) * f(xmin) <= 0) {
            xmax = x;
        } else {
            xmin = x;
        }
    }
    return (xmax + xmin) / 2;
}


