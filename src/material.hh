#pragma once

#include "vec2.hh"
#include "rng.hh"
#include <list>
#include <vector>
#include <string>

struct Wave {
    Vec2 E;
    double H;
    double omega;
    double phi;
    double photon_energy;
};

class Band;

class Particle {
public:
    Vec2 p;
    double r;
    Rng rng;
    Band *band;
    void reset_r();
    Particle(unsigned int seed);
};

struct BandScatteringIntegral {
    double momentum;
    double integral;
};

struct BandScatteringEntry {
    double energy;
    std::list<BandScatteringIntegral> integrals;
};

struct ScatteringResult {
    double momentum;
    double rate;
};

class Band {
public:
    virtual double min_energy() const = 0;
    virtual double energy(double momentum) const = 0;
    virtual double energy(Vec2 const & momentum) const = 0;
    virtual double velocity(double momentum) const = 0;
    virtual Vec2 velocity(Vec2 const & momentum) const = 0;
    virtual std::list<ScatteringResult> acoustic_phonon_scattering(Particle * p) = 0;
    virtual std::list<ScatteringResult> optical_phonon_scattering(Particle * p) = 0;
    Vec2 momentum_scattering(double momentum, Particle * p);
};

class Material {
public:
    std::vector<Band*> bands;
    virtual Particle* create_particle(int seed) = 0;
    virtual double vertical_transition(Particle * p, Band *from, Band *to, Wave const & wave, double de) = 0;
};

namespace materials {

    class Bigraphene : public Material {
    public:
        double delta;
        Bigraphene(double temperature, double delta, double number_of_bands);
        Particle* create_particle(int seed);
        double vertical_transition(Particle * p, Band *from, Band *to, Wave const & wave, double de);
    };

    class Graphene : public Material {
    public:
        Graphene(double temperature, double delta);
        Particle* create_particle(int seed);
        double vertical_transition(Particle * p, Band *from, Band *to, Wave const & wave, double de) { return 0; };
    };
}

template <typename T>
inline double bisection(T f, double xmin, double xmax, double precision) {
    while (xmax - xmin > precision) {
        double x = (xmax + xmin) / 2;
        if (f(x) * f(xmin) <= 0) {
            xmax = x;
        } else {
            xmin = x;
        }
    }
    return (xmax + xmin) / 2;
}


