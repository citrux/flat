#include <cstdio>

#include "material.hh"
#include "physics.hh"

/*
 * Particle
 */
Particle::Particle(Material* m, unsigned int seed) {
  material = m;
  rng = Rng(seed);
  reset_r();
}

void Particle::reset_r() { r = -log(rng.uniform()); }

Vec2 Band::momentum_scattering(double momentum, Particle *particle) {
  FILE *fd = fopen("theta.dat", "a");
  double theta = 2 * pi * particle->rng.uniform();
  double prob = particle->rng.uniform();
  while (0.5 * (std::cos(theta) + 1) < prob) {
    theta = 2 * pi * particle->rng.uniform();
    prob = particle->rng.uniform();
  }
  fprintf(fd, "%f\n", theta);
  fclose(fd);
  double psi = atan2(particle->p.y, particle->p.x);
  return {momentum * std::cos(psi + theta), momentum * std::sin(psi + theta)};
}
