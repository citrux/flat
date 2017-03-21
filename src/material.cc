#include <cmath>
#include "material.hh"
#include "physics.hh"
#include <cstdio>

Vec2 momentum_scattering(float momentum, Particle & particle) {
    float theta = 2 * pi * particle.rng.uniform();
    float prob = particle.rng.uniform();
    while ((std::cos(theta) + 1) / 2 < prob) {
        theta = 2 * pi * particle.rng.uniform();
        prob = particle.rng.uniform();
    }
    float psi = atan2(particle.p.y, particle.p.x); 
    return {momentum * std::cos(psi + theta), momentum * std::sin(psi + theta)};
}

template <typename T>
float bisection(T f, float xmin, float xmax, float precision) {
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

/*
 * Particle
 */
Particle::Particle(unsigned int seed) {
    rng = Rng(seed);
    reset_r();
}

void Particle::reset_r() {
    r = -log(rng.uniform());
}

/*
 * BigrapheneLower
 */
BigrapheneLower::BigrapheneLower(float temperature) {
    const float rho = 2 * 7.7e-8;
    const float Dak = 18;
    const float Dopt = 1.4e9;

    const float T = temperature;
    
    acoustic_phonon_constant = sqr(k * T * Dak) * eV / (2 * hbar *
            rho * sqr(v_s * v_f) * hbar * hbar);
    optical_phonon_constant = k * T * sqr(Dopt) * eV / (4 * hbar *
            optical_phonon_energy * rho * sqr(v_f));

    const float max_momentum = 2;
    const float crit_momentum = delta * std::sqrt(2 - 4 * delta * delta / (gamma * gamma + 4 * delta * delta));
    const float max_energy = energy(max_momentum);
    const float min_energy = gamma * delta / std::sqrt(gamma * gamma + 4 * delta * delta);
    
    energy_samples = 10000;
    momentum_precision = 1e-5;
    
    table = std::vector<BandScatteringEntry>(energy_samples);
    for (int i = 0; i < energy_samples; ++i) {
        float e = (min_energy * (energy_samples - 1 - i) + max_energy * i) / (energy_samples - 1);
        table[i].energy = e;
        if (e < delta) {
            float p = bisection([this,e](float p){return energy(p) - e;}, 0, crit_momentum, momentum_precision);
            table[i].integrals.push_back({p, 2 * pi * p / std::abs(velocity(p))});
        }
        if (e < max_energy) {
            float p = bisection([this,e](float p){return energy(p) - e;}, crit_momentum, max_momentum, momentum_precision);
            table[i].integrals.push_back({p, 2 * pi * p / std::abs(velocity(p))});
        }
    }
}

float BigrapheneLower::min_energy() const {
    return delta * gamma / std::sqrt(gamma * gamma + 4 * delta * delta);
}

float BigrapheneLower::energy(float momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum * momentum;
    return std::sqrt(delta2 + gamma2 / 2 + p2 - std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2));
}

float BigrapheneLower::energy(Vec2 const & momentum) const {
    return energy(momentum.len());
}

float BigrapheneLower::velocity(float momentum) const {
    return velocity(Vec2(momentum, 0)).x;
}

Vec2 BigrapheneLower::velocity(Vec2 const & momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum.dot(momentum);
    return momentum * (1 - 0.5 * (gamma2 + 4 * delta2) / std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) / energy(momentum);
}

std::list<ScatteringResult> BigrapheneLower::acoustic_phonon_scattering(Particle & p) {
    float e = energy(p.p);
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({momentum_scattering(x.momentum, p), acoustic_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum_scattering(momentum, p), acoustic_phonon_constant * 2 * pi * momentum});
    }
    return result;
}

std::list<ScatteringResult> BigrapheneLower::optical_phonon_scattering(Particle & p) {
    float e = energy(p.p) - optical_phonon_energy;
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({momentum_scattering(x.momentum, p), optical_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum_scattering(momentum, p), optical_phonon_constant * 2 * pi * momentum});
    }
    return result;
}