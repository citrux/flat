#include <cmath>
#include "material.hh"
#include "physics.hh"
#include <cstdio>

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

void Particle::isotropic_scattering(float momentum) {
    float theta = 2 * M_PI * rng.uniform();
    float prob = rng.uniform();
    while ((std::cos(theta) + 1) / 2 < prob) {
        theta = 2 * M_PI * rng.uniform();
        prob = rng.uniform();
    }
    float psi = atan2(p.y, p.x); 
    p = {momentum * std::cos(psi + theta), momentum * std::sin(psi + theta)};
}
/*
 * Band
 */
Band::Band(float temperature) {
    const float rho = 2 * 7.7e-8;
    const float Dak = 18;
    const float Dopt = 1.4e9;

    const float T = temperature;
    
    acoustic_phonon_constant = sqr(k * T * Dak) * eV / (2 * hbar *
            rho * sqr(v_s * v_f) * hbar * hbar);
    optical_phonon_constant = k * T * sqr(Dopt) * eV / (4 * hbar *
            optical_phonon_energy * rho * sqr(v_f));

    const float max_momentum = 5;
    const float max_energy = energy(max_momentum);
    const float min_energy = gamma * delta / std::sqrt(gamma * gamma + 4 * delta * delta);
    
    energy_samples = 1000;
    momentum_samples = 1000;
    momentum_precision = 1e-5;
    
    table = std::vector<BandScatteringEntry>(energy_samples);
    for (int i = 0; i < energy_samples; ++i) {
        float e = (min_energy * (energy_samples - 1 - i) + max_energy * i) / (energy_samples - 1);
        table[i].energy = e;
        for (int j = 0; j < momentum_samples - 1; ++j) {
            float p1 = j * max_momentum / (momentum_samples - 1);
            float p2 = (j + 1) * max_momentum / (momentum_samples - 1);
            float p = p = (p1 + p2) / 2;
            if ((energy(p1) - e) * (energy(p2) - e) <= 0) {
                while (p2 - p1 > momentum_precision) {
                    p = (p1 + p2) / 2;
                    if ((energy(p1) - e) * (energy(p) - e) <= 0) {
                        p2 = p;
                    } else {
                        p1 = p;
                    }
                } 
                table[i].integrals.push_back({p, 2 * M_PI * p / std::abs(velocity(p))});
            }
        }
    }
}

bool Band::acoustic_phonon_scattering(Particle & p, float dt) {
    float e = energy(p.p);
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    if (i >= 0 && i < energy_samples) {
        for (auto x: table[i].integrals) {
            p.r -= acoustic_phonon_constant * x.integral * dt;
            if (p.r < 0) {
                p.isotropic_scattering(x.momentum);
                p.reset_r();
                p.band = this;
                return true;
            }
        }
    }
    return false;
}

bool Band::optical_phonon_scattering(Particle & p, float dt) {
    float e = energy(p.p) - optical_phonon_energy;
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    if (i >= 0 && i < energy_samples) {
        for (auto x: table[i].integrals) {
            p.r -= optical_phonon_constant * x.integral * dt;
            if (p.r < 0) {
                p.isotropic_scattering(x.momentum);
                p.reset_r();
                p.band = this;
                return true;
            }
        }
    }
    return false;
}