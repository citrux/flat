#include "../material.hh"
#include "../physics.hh"

namespace materials {
const float gamma = 0.35;
const float optical_phonon_energy = 0.196;

class Lower : public Band {
public:
    const float delta;

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

    Lower(float temperature, float delta=0);
    ~Lower() {};
};

class Upper : public Band {
public:
    const float delta;

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

    Upper(float temperature, float delta=0);
    ~Upper() {};
};

Bigraphene::Bigraphene(float temperature, float delta, float number_of_bands) : delta(delta) {
    Lower *lower = new Lower(temperature, delta);
    bands = {lower};
    if (number_of_bands > 1) {
        Upper *upper = new Upper(temperature, delta);
        bands.push_back(upper);
    }
}


float Bigraphene::vertical_transition(Particle & p, Band *from, Band *to, Wave const & wave) {
    if (from == to) { return 0; }
    float e1 = from->energy(p.p);
    float e2 = to->energy(p.p);
    float d = delta;
    float p2 = p.p.dot(p.p);
    float g = 0.35;
    float g2 = g * g;
    float xi = 1; // fix it!!!
    float theta = std::atan2(p.p.y, p.p.x);
    float lambda = (std::pow(d + e1, 2) - p2) * (std::pow(d + e2, 2) - p2) / g2;
    float numerator = std::pow(e1 * e1 - d * d + lambda, 2) * (wave.E.x * wave.E.x + wave.E.y * wave.E.y - 2 * xi * wave.E.x * wave.E.y * std::sin(wave.phi)) +
        std::pow(e1 - d, 2) / (e2 - d) * std::pow(e2 * e2 - d * d + lambda, 2) * (wave.E.x * wave.E.x + wave.E.y * wave.E.y + 2 * xi * wave.E.x * wave.E.y * std::sin(wave.phi)) +
        2 * (e1 - d) / (e2 - d) * (e1 * e1 - d * d + lambda) * (e2 * e2 - d * d + lambda) * ((wave.E.x * wave.E.x - wave.E.y * wave.E.y) * std::cos(2 * theta) + 2 * wave.E.x * wave.E.y * std::sin(2 * theta));
    float denomenator = std::pow(wave.omega, 2) * std::pow(g, -4) * (g * g * p2 * std::pow(p2 - 4 * d * d, 2) * std::pow(e1 + d, -2) * std::pow(e2 * e2 - d * d, -2) * (std::pow(d + e1, 2) + p2) + (1 + p2 * std::pow(p2 - 4 * d * d, 2) * std::pow(e1 + d, -2) * std::pow(e2 * e2 - d * d, -2)) * std::pow(std::pow(e1 + d,2) - p2, 2)) * (g * g * (std::pow(e2 + d, 2) + p2) + (1 + p2 * std::pow(e2 - d, -2)) * std::pow(std::pow(e2 + d, 2) - p2, 2));
    return pi / 8 / hbar * numerator / denomenator * dirac_delta(e2 - e1 - hbar * wave.omega, 0);
}

const float rho = 2 * 7.7e-8;
const float Dak = 18;
const float Dopt = 1.4e9;

/*
 * Lower
 */
Lower::Lower(float temperature, float _delta) : delta(_delta) {
    const float T = temperature;

    acoustic_phonon_constant = k * T * sqr(Dak) * eV / (2 * hbar *
            rho * sqr(v_s * v_f) * hbar * hbar) / 2 / pi;
    optical_phonon_constant = sqr(Dopt) * eV / (4 * hbar *
            optical_phonon_energy * rho * sqr(v_f)) / 2 / pi;

    const float max_momentum = 2;
    const float crit_momentum = delta * std::sqrt(2 - 4 * delta * delta / (gamma * gamma + 4 * delta * delta));
    const float max_energy = energy(max_momentum);

    energy_samples = 20000;
    momentum_precision = 1e-5;

    table = std::vector<BandScatteringEntry>(energy_samples);
    table[0].energy = min_energy();
    table[0].integrals.push_back({crit_momentum,
                    2.0f * pi * crit_momentum / std::sqrt(0 + 1e-3f)});
    for (int i = 1; i < energy_samples; ++i) {
        float e = (min_energy() * (energy_samples - 1 - i) + max_energy * i) / (energy_samples - 1);
        table[i].energy = e;
        if (e < delta) {
            float p = bisection([this,e](float p){return energy(p) - e;}, 0, crit_momentum, momentum_precision);
            table[i].integrals.push_back({p,
                    2.0f * pi * p / std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
        }
        if (e < max_energy) {
            float p = bisection([this,e](float p){return energy(p) - e;}, crit_momentum, max_momentum, momentum_precision);
            table[i].integrals.push_back({p,
                    2.0f * pi * p / std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
        }
    }
}

float Lower::min_energy() const {
    return delta * gamma / std::sqrt(gamma * gamma + 4 * delta * delta);
}

float Lower::energy(float momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum * momentum;
    return std::sqrt(((p2 - delta2) * (p2 - delta2) + gamma2 * delta2) / (delta2 + gamma2 / 2 + p2 + std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)));
}

float Lower::energy(Vec2 const & momentum) const {
    return energy(momentum.len());
}

float Lower::velocity(float momentum) const {
    return velocity(Vec2(momentum, 0)).x;
}

Vec2 Lower::velocity(Vec2 const & momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum.dot(momentum);
    const float nrg = energy(momentum);
    if (nrg == 0) {
        return {0, 0};
    }
    return momentum * (1 - 0.5 * (gamma2 + 4 * delta2) / std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) / nrg;
}

std::list<ScatteringResult> Lower::acoustic_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p);
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({x.momentum, acoustic_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum});
    }
    return result;
}

std::list<ScatteringResult> Lower::optical_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p) - optical_phonon_energy;
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({x.momentum, optical_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum});
    }
    return result;
}
/*
 * Upper
 */
Upper::Upper(float temperature, float _delta) : delta(_delta) {
    const float T = temperature;

    acoustic_phonon_constant = k * T * sqr(Dak) * eV / (2 * hbar *
            rho * sqr(v_s * v_f) * hbar * hbar) / 2 / pi;
    optical_phonon_constant = sqr(Dopt) * eV / (4 * hbar *
            optical_phonon_energy * rho * sqr(v_f)) / 2 / pi;

    const float max_momentum = 2;
    const float max_energy = energy(max_momentum);

    energy_samples = 10000;
    momentum_precision = 1e-5;

    table = std::vector<BandScatteringEntry>(energy_samples);
    for (int i = 0; i < energy_samples; ++i) {
        float e = (min_energy() * (energy_samples - 1 - i) + max_energy * i) / (energy_samples - 1);
        table[i].energy = e;
        float p = bisection([this,e](float p){return energy(p) - e;}, 0, max_momentum, momentum_precision);
        table[i].integrals.push_back({p, 2 * pi * p / std::abs(velocity(p))});
    }
}

float Upper::min_energy() const {
    return std::sqrt(delta*delta + gamma*gamma);
}

float Upper::energy(float momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum * momentum;
    return std::sqrt(delta2 + gamma2 / 2 + p2 + std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2));
}

float Upper::energy(Vec2 const & momentum) const {
    return energy(momentum.len());
}

float Upper::velocity(float momentum) const {
    return velocity(Vec2(momentum, 0)).x;
}

Vec2 Upper::velocity(Vec2 const & momentum) const {
    const float delta2 = delta * delta;
    const float gamma2 = gamma * gamma;
    const float gamma4 = gamma2 * gamma2;
    const float p2 = momentum.dot(momentum);
    const float nrg = energy(momentum);
    return momentum * (1 + 0.5 * (gamma2 + 4 * delta2) / std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) / nrg;
}

std::list<ScatteringResult> Upper::acoustic_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p);
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({x.momentum, acoustic_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum});
    }
    return result;
}

std::list<ScatteringResult> Upper::optical_phonon_scattering(Particle & p) {
    float e = p.band->energy(p.p) - optical_phonon_energy;
    int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
    std::list<ScatteringResult> result;
    if (i >= 0 && i < energy_samples) {
        for (auto const & x: table[i].integrals) {
            result.push_back({x.momentum, optical_phonon_constant * x.integral});
        }
    }
    if (i >= energy_samples-1) {
        float momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
        result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum});
    }
    return result;
}
}
