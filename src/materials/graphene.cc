#include "../material.hh"
#include "../physics.hh"

namespace materials {
const double rho = 7.7e-8;
const double Dak = 18;
const double Dopt = 1.4e9;

class Bnd : public Band {
public:
  const double delta;

  const double optical_phonon_energy = 0.196;
  double acoustic_phonon_constant;
  double optical_phonon_constant;

  double min_energy() const;
  double energy(double momentum) const;
  double energy(Vec2 const &momentum) const;
  double velocity(double momentum) const;
  Vec2 velocity(Vec2 const &momentum) const;
  std::list<ScatteringResult> acoustic_phonon_scattering(Particle *p);
  std::list<ScatteringResult> optical_phonon_scattering(Particle *p);

  Bnd(double temperature, double delta = 0.13);
  ~Bnd(){};
};

Graphene::Graphene(double temperature, double delta) {
  Bnd *band = new Bnd(temperature, delta);
  bands = {band};
}

Particle *Graphene::create_particle(int seed) { return new Particle(this, seed); }

Bnd::Bnd(double temperature, double _delta) : delta(_delta) {
  const double T = temperature;
  acoustic_phonon_constant = k * T * sqr(Dak) * eV /
                             (2 * hbar * rho * sqr(v_s * v_f) * hbar * hbar) /
                             2 / pi;
  optical_phonon_constant =
      sqr(Dopt) * eV / (4 * hbar * optical_phonon_energy * rho * sqr(v_f)) / 2 /
      pi;
}

double Bnd::min_energy() const { return delta; }

double Bnd::energy(double momentum) const {
  const double delta2 = delta * delta;
  const double p2 = momentum * momentum;
  return std::sqrt(delta2 + p2);
}

double Bnd::energy(Vec2 const &momentum) const {
  return energy(momentum.len());
}

double Bnd::velocity(double momentum) const {
  return velocity(Vec2(momentum, 0)).x;
}

Vec2 Bnd::velocity(Vec2 const &momentum) const {
  const double nrg = energy(momentum);
  return momentum / nrg;
}

std::list<ScatteringResult> Bnd::acoustic_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p);
  std::list<ScatteringResult> result;
  if (e > delta) {
    double momentum = std::sqrt(e * e - delta * delta);
    result.push_back({momentum, acoustic_phonon_constant * 2 * pi * momentum /
                                    velocity(momentum)});
  }
  return result;
}

std::list<ScatteringResult> Bnd::optical_phonon_scattering(Particle *p) {
  double e = p->band->energy(p->p) - optical_phonon_energy;
  std::list<ScatteringResult> result;
  if (e > delta) {
    double momentum = std::sqrt(e * e - delta * delta);
    result.push_back({momentum, optical_phonon_constant * 2 * pi * momentum /
                                    velocity(momentum)});
  }
  return result;
}
}
