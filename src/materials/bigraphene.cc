#include "../material.hh"
#include "../physics.hh"

namespace materials {
const double gamma = 0.35;
const double optical_phonon_energy = 0.196;

class Lower : public Band {
public:
  const double delta;

  double acoustic_phonon_constant;
  double optical_phonon_constant;

  std::vector<BandScatteringEntry> table;
  int energy_samples;
  int momentum_samples;
  double momentum_precision;

  double min_energy() const;
  double energy(double momentum) const;
  double energy(Vec2 const &momentum) const;
  double velocity(double momentum) const;
  Vec2 velocity(Vec2 const &momentum) const;
  std::list<ScatteringResult> acoustic_phonon_scattering(Particle *p);
  std::list<ScatteringResult> optical_phonon_scattering(Particle *p);

  Lower(double temperature, double delta = 0);
  ~Lower(){};
};

class Upper : public Band {
public:
  const double delta;

  double acoustic_phonon_constant;
  double optical_phonon_constant;

  std::vector<BandScatteringEntry> table;
  int energy_samples;
  int momentum_samples;
  double momentum_precision;

  double min_energy() const;
  double energy(double momentum) const;
  double energy(Vec2 const &momentum) const;
  double velocity(double momentum) const;
  Vec2 velocity(Vec2 const &momentum) const;
  std::list<ScatteringResult> acoustic_phonon_scattering(Particle *p);
  std::list<ScatteringResult> optical_phonon_scattering(Particle *p);

  Upper(double temperature, double delta = 0);
  ~Upper(){};
};

Bigraphene::Bigraphene(double temperature, double delta, double number_of_bands)
    : delta(delta) {
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
  BigrapheneParticle(int seed) : Particle(seed) {
    valley = (rng.uniform() > 0.5) ? 1 : -1;
  }
};

Particle *Bigraphene::create_particle(int seed) {
  return new BigrapheneParticle(seed);
}

double Bigraphene::vertical_transition(Particle *p, Band *dest,
                                       Wave const &wave, double de) {
  if (p->band == dest) {
    return 0;
  }
  double e1 = p->band->energy(p->p);
  double e2 = dest->energy(p->p);
  double d = delta;
  double p2 = p->p.dot(p->p);
  double g = gamma;
  double g2 = g * g;
  double xi = ((BigrapheneParticle *)p)->valley; // fix it!!!
  double theta = std::atan2(p->p.y, p->p.x);
  double a = g2 * (e1 - d) * (e2 - d) / wave.omega;
  double b =
      (std::pow(e1 + d, 2) - p2) * (std::pow(e2 + d, 2) - p2) / wave.omega;
  Vec2 c1 = wave.E.x *
            Vec2((a + b) * (e1 + e2) + 2 * (a - b) * d, (b - a) * (e2 - e1)) *
            Vec2(std::cos(xi * theta), std::sin(xi * theta));
  Vec2 c2 = wave.E.y *
            Vec2((a + b) * (e1 + e2) + 2 * (a - b) * d, (b - a) * (e2 - e1)) *
            Vec2(std::sin(xi * theta), -std::cos(xi * theta));
  Vec2 c = c1 + c2.rotate(-wave.phi);
  double wtf =
      p2 * std::pow((p2 - 4 * d * d) / (e1 + d) / (e2 * e2 - d * d), 2);
  double norm1 = g2 * std::pow(e1 - d, 2) + g2 * std::pow(e1 + d, 2) * wtf +
                 (wtf + 1) * std::pow(std::pow(e1 + d, 2) - p2, 2);
  double norm2 =
      g2 * p2 * std::pow(e2 - d, 2) + g2 * std::pow(e2 * e2 - d * d, 2) +
      (std::pow(e2 - d, 2) + p2) * std::pow(std::pow(e2 + d, 2) - p2, 2);
  double norm = norm1 * norm2;
  return pi / 8 / hbar * c.dot(c) / norm *
         dirac_delta(std::abs(e2 - e1) - wave.photon_energy, de);
}

const double rho = 2 * 7.7e-8;
const double Dak = 18;
const double Dopt = 1.4e9;

/*
 * Lower
 */
Lower::Lower(double temperature, double _delta) : delta(_delta) {
  const double T = temperature;

  acoustic_phonon_constant = k * T * sqr(Dak) * eV /
                             (2 * hbar * rho * sqr(v_s * v_f) * hbar * hbar) /
                             2 / pi;
  optical_phonon_constant =
      sqr(Dopt) * eV / (4 * hbar * optical_phonon_energy * rho * sqr(v_f)) / 2 /
      pi;

  const double max_momentum = 2;
  const double crit_momentum =
      delta *
      std::sqrt(2 - 4 * delta * delta / (gamma * gamma + 4 * delta * delta));
  const double max_energy = energy(max_momentum);

  energy_samples = 20000;
  momentum_precision = 1e-5;

  table = std::vector<BandScatteringEntry>(energy_samples);
  table[0].energy = min_energy();
  table[0].integrals.push_back(
      {crit_momentum, 2.0f * pi * crit_momentum / std::sqrt(0 + 1e-3f)});
  for (int i = 1; i < energy_samples; ++i) {
    double e = (min_energy() * (energy_samples - 1 - i) + max_energy * i) /
               (energy_samples - 1);
    table[i].energy = e;
    if (e < delta) {
      double p = bisection([this, e](double p) { return energy(p) - e; }, 0,
                           crit_momentum, momentum_precision);
      table[i].integrals.push_back(
          {p, 2.0f * pi * p /
                  (double)std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
    }
    if (e < max_energy) {
      double p = bisection([this, e](double p) { return energy(p) - e; },
                           crit_momentum, max_momentum, momentum_precision);
      table[i].integrals.push_back(
          {p, 2.0f * pi * p /
                  (double)std::sqrt(std::pow(velocity(p), 2) + 1e-3f)});
    }
  }
}

double Lower::min_energy() const {
  return delta * gamma / std::sqrt(gamma * gamma + 4 * delta * delta);
}

double Lower::energy(double momentum) const {
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double gamma4 = gamma2 * gamma2;
  const double p2 = momentum * momentum;
  return std::sqrt(((p2 - delta2) * (p2 - delta2) + gamma2 * delta2) /
                   (delta2 + gamma2 / 2 + p2 +
                    std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)));
}

double Lower::energy(Vec2 const &momentum) const {
  return energy(momentum.len());
}

double Lower::velocity(double momentum) const {
  return velocity(Vec2(momentum, 0)).x;
}

Vec2 Lower::velocity(Vec2 const &momentum) const {
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double gamma4 = gamma2 * gamma2;
  const double p2 = momentum.dot(momentum);
  const double nrg = energy(momentum);
  if (nrg == 0) {
    return {0, 0};
  }
  return momentum * (1 -
                     0.5 * (gamma2 + 4 * delta2) /
                         std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) /
         nrg;
}

std::list<ScatteringResult> Lower::acoustic_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p);
  int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
  std::list<ScatteringResult> result;
  if (i >= 0 && i < energy_samples) {
    for (auto const &x : table[i].integrals) {
      result.push_back({x.momentum, acoustic_phonon_constant * x.integral});
    }
  }
  if (i >= energy_samples - 1) {
    double momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
    result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum});
  }
  return result;
}

std::list<ScatteringResult> Lower::optical_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p) - optical_phonon_energy;
  int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
  std::list<ScatteringResult> result;
  if (i >= 0 && i < energy_samples) {
    for (auto const &x : table[i].integrals) {
      result.push_back({x.momentum, optical_phonon_constant * x.integral});
    }
  }
  if (i >= energy_samples - 1) {
    double momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
    result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum});
  }
  return result;
}
/*
 * Upper
 */
Upper::Upper(double temperature, double _delta) : delta(_delta) {
  const double T = temperature;

  acoustic_phonon_constant = k * T * sqr(Dak) * eV /
                             (2 * hbar * rho * sqr(v_s * v_f) * hbar * hbar) /
                             2 / pi;
  optical_phonon_constant =
      sqr(Dopt) * eV / (4 * hbar * optical_phonon_energy * rho * sqr(v_f)) / 2 /
      pi;

  const double max_momentum = 2;
  const double max_energy = energy(max_momentum);

  energy_samples = 10000;
  momentum_precision = 1e-5;

  table = std::vector<BandScatteringEntry>(energy_samples);
  for (int i = 0; i < energy_samples; ++i) {
    double e = (min_energy() * (energy_samples - 1 - i) + max_energy * i) /
               (energy_samples - 1);
    table[i].energy = e;
    double p = bisection([this, e](double p) { return energy(p) - e; }, 0,
                         max_momentum, momentum_precision);
    table[i].integrals.push_back({p, 2 * pi * p / std::abs(velocity(p))});
  }
}

double Upper::min_energy() const {
  return std::sqrt(delta * delta + gamma * gamma);
}

double Upper::energy(double momentum) const {
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double gamma4 = gamma2 * gamma2;
  const double p2 = momentum * momentum;
  return std::sqrt(delta2 + gamma2 / 2 + p2 +
                   std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2));
}

double Upper::energy(Vec2 const &momentum) const {
  return energy(momentum.len());
}

double Upper::velocity(double momentum) const {
  return velocity(Vec2(momentum, 0)).x;
}

Vec2 Upper::velocity(Vec2 const &momentum) const {
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double gamma4 = gamma2 * gamma2;
  const double p2 = momentum.dot(momentum);
  const double nrg = energy(momentum);
  return momentum * (1 +
                     0.5 * (gamma2 + 4 * delta2) /
                         std::sqrt(gamma4 / 4 + (gamma2 + 4 * delta2) * p2)) /
         nrg;
}

std::list<ScatteringResult> Upper::acoustic_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p);
  int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
  std::list<ScatteringResult> result;
  if (i >= 0 && i < energy_samples) {
    for (auto const &x : table[i].integrals) {
      result.push_back({x.momentum, acoustic_phonon_constant * x.integral});
    }
  }
  if (i >= energy_samples - 1) {
    double momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
    result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum});
  }
  return result;
}

std::list<ScatteringResult> Upper::optical_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p) - optical_phonon_energy;
  int i = (e - table[0].energy) / (table[1].energy - table[0].energy);
  std::list<ScatteringResult> result;
  if (i >= 0 && i < energy_samples) {
    for (auto const &x : table[i].integrals) {
      result.push_back({x.momentum, optical_phonon_constant * x.integral});
    }
  }
  if (i >= energy_samples - 1) {
    double momentum = e + .5 * std::sqrt(gamma * gamma + 4 * delta * delta);
    result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum});
  }
  return result;
}
}
