#include <cmath>
#include "material.hh"

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

void Particle::anisotropic_scattering(float momentum) {
    float theta = 2 * M_PI * rng.uniform();
    p = {momentum * cos(theta), momentum * sin(theta)};
}
/*
 * Band
 */
Band::Band() {
    gamma = 0.35;
    delta = 1e-3;

    float max_momentum = 10;
    float max_energy = energy(max_momentum);
    float min_energy = gamma * delta / sqrt(gamma * gamma + 4 * delta * delta);
    
    energy_samples = 1000;
    momentum_samples = 1000;
    momentum_precision = 1e-5;
    
    table = new list<BandScatteringIntegral>[energy_samples];
    for (int i = 0; i < energy_samples; ++i) {
        float e = (min_energy * (energy_samples - 1 - i) + max_energy * i) / (energy_samples - 1);
        for (int j = 0; j < momentum_samples - 1; ++j) {
            float p1 = j * max_momentum / (momentum_samples - 1);
            float p2 = (j + 1) * max_momentum / (momentum_samples - 1);
            if ((energy(p1) - e) * (energy(p2) - e) <= 0) {
                while (p2 - p1 > momentum_precision) {
                    p = (p1 + p2) / 2;
                    if ((energy(p1) - e) * (energy(p) - e) <= 0) {
                        p2 = p;
                    } else {
                        p1 = p;
                    }
                } 
                table[i].push_back({e, p, 2 * M_PI * p / velocity(p)});
            }
        }
    }
}

float Band::energy(float momentum) {
    float delta2 = delta * delta;
    float gamma2 = gamma * gamma;
    float gamma4 = gamma2 * gamma2;
    float p2 = momentum * momentum;
    return sqrt(delta2 + gamma2 / 2 + p2 - sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2));
}

float Band::energy(Vec2 momentum) {
    return energy(momentum.len())
}

float Band::velocity(float momentum) {
    return velocity(Vec2(momentum, 0)).len();
}

Vec2 Band::velocity(Vec2 momentum) {
    float delta2 = delta * delta;
    float gamma2 = gamma * gamma;
    float gamma4 = gamma2 * gamma2;
    float p2 = momentum * momentum;
    return momentum * (1 - 0.5 * (gamma2 + 4 * delta2) / sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) / energy(momentum);
}

bool acoustic_phonon_scattering(Particle & p) {
    float e = energy(p.p);
    if (e > table[0].at(0).energy && e < table[energy_samples-1].at(0).energy) {
        int i = 0;
        int j = energy_samples - 1;
        while (j - i > 1) {
            int k = (i + j) / 2;
            if (e >= table[k].at(0).energy) {
                i = k;
            } else {
                j = k;
            }
        }
        for (auto x: table[i]) {
            r -= acoustic_phonon_constant * x.integral;
            if (r < 0) {
                p.anisotropic_scattering(x.momentum);
                r.reset_r();
                return true;
            }
        }
    }
    return false;
}
bool optical_phonon_scattering(Particle & p) {
    float e = energy(p.p) - optical_phonon_energy;
    if (e > table[0].at(0).energy && e < table[energy_samples-1].at(0).energy) {
        int i = 0;
        int j = energy_samples - 1;
        while (j - i > 1) {
            int k = (i + j) / 2;
            if (e >= table[k].at(0).energy) {
                i = k;
            } else {
                j = k;
            }
        }
        for (auto x: table[i]) {
            r -= optical_phonon_constant * x.integral;
            if (r < 0) {
                p.anisotropic_scattering(x.momentum);
                r.reset_r();
                return true;
            }
        }
    }
    return false;
}