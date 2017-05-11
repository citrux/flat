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
    std::list<ScatteringResult> acoustic_phonon_scattering(Particle * p);
    std::list<ScatteringResult> optical_phonon_scattering(Particle * p);

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
    std::list<ScatteringResult> acoustic_phonon_scattering(Particle * p);
    std::list<ScatteringResult> optical_phonon_scattering(Particle * p);

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

class BigrapheneParticle : public Particle {
public:
    int valley;
    BigrapheneParticle(int seed) : Particle(seed) { valley = (rng.uniform() > 0.5) ? 1 : -1;}
};

Particle* Bigraphene::create_particle(int seed) {return new BigrapheneParticle(seed);}

float Bigraphene::vertical_transition(Particle * p, Band *from, Band *to, Wave const & wave, float de) {
    if (from == to) { return 0; }
    float e1 = from->energy(p->p);
    float e2 = to->energy(p->p);
    float d = delta;
    float p2 = p->p.dot(p->p);
    float g = gamma;
    float g2 = g * g;
    float xi = ((BigrapheneParticle *)p)->valley; // fix it!!!
    float theta = std::atan2(p->p.y, p->p.x);
    float a = g2 * (e1 - d) * (e2 - d) / wave.omega;
    float b = (std::pow(e1 + d, 2) - p2) * (std::pow(e2 + d, 2) - p2) / wave.omega;
    Vec2 c1 = wave.E.x * Vec2((a+b)*(e1+e2)+2*(a-b)*d, (b-a)*(e2-e1)) * Vec2(std::cos(xi * theta), std::sin(xi * theta));
    Vec2 c2 = wave.E.y * Vec2((a+b)*(e1+e2)+2*(a-b)*d, (b-a)*(e2-e1)) * Vec2(std::sin(xi * theta), -std::cos(xi * theta));
    Vec2 c = c1 + c2.rotate(-wave.phi);
    float wtf = p2 * std::pow((p2 - 4 * d * d) / (e1 + d) / (e2 * e2 - d * d), 2);
    float norm1 = g2 * std::pow(e1 - d, 2) + g2 * std::pow(e1 + d, 2) * wtf + (wtf + 1) * std::pow(std::pow(e1 + d, 2) - p2, 2);
    float norm2 = g2 * p2 * std::pow(e2 - d, 2) + g2 * std::pow(e2 * e2 - d * d, 2) + (std::pow(e2-d,2) + p2) * std::pow(std::pow(e2 + d, 2) - p2, 2);
    float norm = norm1 * norm2;
    return pi / 8 / hbar * c.dot(c) / norm * dirac_delta(std::abs(e2 - e1) - wave.photon_energy, de);
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
                    2.0f * pi * p / (float) std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
        }
        if (e < max_energy) {
            float p = bisection([this,e](float p){return energy(p) - e;}, crit_momentum, max_momentum, momentum_precision);
            table[i].integrals.push_back({p,
                    2.0f * pi * p / (float)std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
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

std::list<ScatteringResult> Lower::acoustic_phonon_scattering(Particle * p) {
    float e = p->band->energy(p->p);
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

std::list<ScatteringResult> Lower::optical_phonon_scattering(Particle * p) {
    float e = p->band->energy(p->p) - optical_phonon_energy;
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

std::list<ScatteringResult> Upper::acoustic_phonon_scattering(Particle * p) {
    float e = p->band->energy(p->p);
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

std::list<ScatteringResult> Upper::optical_phonon_scattering(Particle * p) {
    float e = p->band->energy(p->p) - optical_phonon_energy;
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
