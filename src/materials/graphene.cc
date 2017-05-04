#include "../material.hh"
#include "../physics.hh"

namespace materials {
const float rho = 7.7e-8;
const float Dak = 18;
const float Dopt = 1.4e9;

class Bnd : public Band {
public:
    const float delta;

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

    Bnd(float temperature, float delta=0.13);
    ~Bnd() {};
};

Graphene::Graphene(float temperature, float delta) {
    Bnd *band = new Bnd(temperature, delta);
    bands = {band};
}


Bnd::Bnd(float temperature, float _delta) : delta(_delta) {
    const float T = temperature;
    acoustic_phonon_constant = k * T * sqr(Dak) * eV / (2 * hbar *
            rho * sqr(v_s * v_f) * hbar * hbar) / 2 / pi;
    optical_phonon_constant = sqr(Dopt) * eV / (4 * hbar *
            optical_phonon_energy * rho * sqr(v_f)) / 2 / pi;
}

float Bnd::min_energy() const {
    return delta;
}

float Bnd::energy(float momentum) const {
    const float delta2 = delta * delta;
    const float p2 = momentum * momentum;
    return std::sqrt(delta2 + p2);
}

float Bnd::energy(Vec2 const & momentum) const {
    return energy(momentum.len());
}

float Bnd::velocity(float momentum) const {
    return velocity(Vec2(momentum, 0)).x;
}

Vec2 Bnd::velocity(Vec2 const & momentum) const {
    const float nrg = energy(momentum);
    return momentum / nrg;
}

std::list<ScatteringResult> Bnd::acoustic_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p);
    std::list<ScatteringResult> result;
    if (e > delta) {
        float momentum = std::sqrt(e * e - delta * delta);
        result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum / velocity(momentum)});
    }
    return result;
}

std::list<ScatteringResult> Bnd::optical_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p) - optical_phonon_energy;
    std::list<ScatteringResult> result;
    if (e > delta) {
        float momentum = std::sqrt(e * e - delta * delta);
        result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum / velocity(momentum)});
    }
    return result;
}
}
