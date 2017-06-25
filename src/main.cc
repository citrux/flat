#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

#include "material.hh"
#include "physics.hh"
#include "stats.hh"

template <typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

int main(int argc, char const *argv[]) {
  bool dumping = false;
  bool verbose = false;
  char dumpname[256] = "dump.bin";
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-d")) {
      dumping = true;
      if (i + 1 < argc)
        strcpy(dumpname, argv[i + 1]);
    }
    if (!strcmp(argv[i], "-v")) {
      verbose = true;
    }
  }

  std::string material, material_name;
  Vec2 Ec;
  double Hc, T;
  int n;
  int number_of_waves;
  double dt, all_time;

  /* input */
  std::cin >> material;
  std::cin >> Ec.x >> Ec.y;
  std::cin >> Hc;
  std::cin >> number_of_waves;
  std::vector<Wave> waves(number_of_waves);
  for (int i = 0; i < number_of_waves; i++) {
    std::cin >> waves[i].E.x >> waves[i].E.y;
    std::cin >> waves[i].H;
    std::cin >> waves[i].omega;
    std::cin >> waves[i].phi;
  }
  std::cin >> T;
  std::cin >> n;
  std::cin >> dt;
  std::cin >> all_time;

  auto material_params = split(material, '_');
  material_name = material_params[0];

  if (material_name != "bigraphene" && material_name != "graphene") {
    puts("incorrect material");
    exit(1);
  }
  if (verbose) {
    printf("threads: %d\n", omp_get_num_threads());
    printf("Field configuration:\n");
    printf("Ec: {%e, %e}\n", Ec.x, Ec.y);
    printf("Hc: %e\n", Hc);
  }

  const double field_dimensionless_factor = e * v_f * dt / eV;
  Ec *= field_dimensionless_factor;
  Hc *= v_f / c * field_dimensionless_factor;
  for (auto &w : waves) {
    w.E *= field_dimensionless_factor;
    w.H *= v_f / c * field_dimensionless_factor;
    w.photon_energy = w.omega * hbar;
    w.omega *= dt;
  }

  int s = all_time / dt;

  std::vector<unsigned int> seeds(n);
  double *dump;
  if (dumping) {
    dump = new double[2 * n * s];
  }
  srand(time(nullptr));
  for (int i = 0; i < n; ++i) {
    seeds[i] = rand();
  }
  Material *mat;
  if (material_name == "bigraphene") {
    double delta = 0;
    int number_of_bands = 1;
    if (material_params.size() > 1) {
      delta = atof(material_params[1].c_str());
    }
    if (material_params.size() > 2) {
      number_of_bands = atoi(material_params[2].c_str());
    }
    mat = new materials::Bigraphene(T, delta, number_of_bands);
  } else {
    double delta = 0;
    if (material_params.size() > 1) {
      delta = atof(material_params[1].c_str());
    }
    mat = new materials::Graphene(T, delta);
  }
  double de = 0;
  std::vector<Data> datas(n, Data(waves.size(), mat->bands.size()));
  if (verbose) {
    puts("start calculation");
  }
/* first run: just calculate tau for calculating de */
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    Particle *particle = mat->create_particle(seeds[i]);
    particle->band = mat->bands[0];
    particle->p = {0, 0};
    /* Boltzmann-distributed initial condition */
    double prob;
    do {
      double p = particle->rng.uniform();
      prob = particle->rng.uniform();
      double theta = 2 * pi * particle->rng.uniform();
      particle->p = {p * std::cos(theta), p * std::sin(theta)};
    } while (exp(-(particle->band->energy(particle->p) -
                   particle->band->min_energy()) /
                 k / T) < prob);
    datas[i].power.assign(number_of_waves + 1, 0);
    datas[i].population.assign(mat->bands.size(), 0);

    /* simulation */
    int t = 0;
    while (particle->r > 0 && t < s) {
      Vec2 v = particle->band->velocity(particle->p);
      Vec2 f = Ec + Vec2(v.y, -v.x) * Hc;
      for (auto w : waves) {
        f += w.E * Vec2(std::cos(w.omega * t), std::cos(w.omega * t + w.phi)) +
             Vec2(v.y, -v.x) * w.H * std::cos(w.omega * t);
      }
      int band_index =
          std::find(mat->bands.begin(), mat->bands.end(), particle->band) -
          mat->bands.begin();
      particle->p += f;
      double wsum = 0;
      for (auto const band : mat->bands) {
        for (auto const &result : band->acoustic_phonon_scattering(particle))
          wsum += result.rate;
        for (auto const &result : band->optical_phonon_scattering(particle))
          wsum += result.rate;
      }
      particle->r -= wsum * dt;
      ++t;
    }
    datas[i].tau = t * dt;
  }

  de = hbar / mean(datas).tau;
  if (verbose) {
    printf("first run: tau = %e s, de = %e eV\n", hbar / de, de);
  }

/* second run */
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    Particle *particle = mat->create_particle(seeds[i]);
    particle->band = mat->bands[0];
    particle->p = {0, 0};
    /* Boltzmann-distributed initial condition */
    double prob;
    do {
      double p = particle->rng.uniform();
      prob = particle->rng.uniform();
      double theta = 2 * pi * particle->rng.uniform();
      particle->p = {p * std::cos(theta), p * std::sin(theta)};
    } while (exp(-(particle->band->energy(particle->p) -
                   particle->band->min_energy()) /
                 k / T) < prob);
    datas[i].power.assign(number_of_waves + 1, 0);
    datas[i].population.assign(mat->bands.size(), 0);

    /* simulation */
    for (int t = 0; t < s; ++t) {
      if (dumping) {
        dump[2 * t * n + 2 * i] = particle->p.x;
        dump[2 * t * n + 2 * i + 1] = particle->p.y;
      }
      Vec2 v = particle->band->velocity(particle->p);
      Vec2 f = Ec + Vec2(v.y, -v.x) * Hc;
      for (auto w : waves) {
        f += w.E * Vec2(std::cos(w.omega * t), std::cos(w.omega * t + w.phi)) +
             Vec2(v.y, -v.x) * w.H * std::cos(w.omega * t);
      }
      if (std::isnan(f.x) || std::isnan(f.y)) {
        puts("force is nan!");
        printf("%e %e\n", f.x, f.y);
        printf("%e %e\n", v.x, v.y);
        printf("%e %e\n", particle->p.x, particle->p.y);
        printf("%e\n", particle->band->energy(particle->p));
        exit(1);
      }
      if (std::isnan(v.x) || std::isnan(v.y)) {
        puts("wtf");
        printf("%e %e\n", v.x, v.y);
        printf("%e %e\n", particle->p.x, particle->p.y);
        exit(1);
      }
      int band_index =
          std::find(mat->bands.begin(), mat->bands.end(), particle->band) -
          mat->bands.begin();
      datas[i].population[band_index] += 1;
      datas[i].v += v;
      datas[i].power[0] += v.dot(Ec);
      for (int j = 0; j < number_of_waves; j++) {
        const double Ex = waves[j].E.x;
        const double Ey = waves[j].E.y;
        const double omega = waves[j].omega;
        const double phi = waves[j].phi;
        Vec2 E = {Ex * std::cos(omega * t), Ey * std::cos(omega * t + phi)};
        datas[i].power[j + 1] += v.dot(E);
      }
      particle->p += f;
      double wsum = 0;
      for (auto const band : mat->bands) {
        for (auto const &result : band->acoustic_phonon_scattering(particle))
          wsum += result.rate;
        for (auto const &result : band->optical_phonon_scattering(particle))
          wsum += result.rate;
        for (auto const &wave : waves)
          wsum += mat->vertical_transition(particle, band, wave, de);
      }
      particle->r -= wsum * dt;
      if (particle->r < 0) {
        double w = particle->rng.uniform() * wsum;
        for (auto const band : mat->bands) {
          for (auto const &result :
               band->acoustic_phonon_scattering(particle)) {
            w -= result.rate;
            if (w < 0) {
              ++datas[i].acoustic_phonon_scattering_count;
              particle->reset_r();
              particle->band = band;
              particle->p =
                  band->momentum_scattering(result.momentum, particle);
              goto end;
            }
          }
          for (auto const &result : band->optical_phonon_scattering(particle)) {
            w -= result.rate;
            if (w < 0) {
              ++datas[i].optical_phonon_scattering_count;
              particle->reset_r();
              particle->band = band;
              particle->p =
                  band->momentum_scattering(result.momentum, particle);
              goto end;
            }
          }
          for (auto const &wave : waves) {
            w -= mat->vertical_transition(particle, band, wave, de);
            if (w < 0) {
              ++datas[i].vertical_transitions_count;
              int wave_index = 0;
              for (int i = 0; i < waves.size(); ++i) {
                if (&wave == &waves[i]) {
                  wave_index = i;
                }
              }
              datas[i].power[wave_index + 1] +=
                  band->energy(particle->p) -
                  particle->band->energy(particle->p);
              particle->reset_r();
              particle->band = band;
              goto end;
            }
          }
        }
      }
    end:
      wsum = 0;
    }
    datas[i].v /= s;
    datas[i].power /= s;
    datas[i].population /= s;
    datas[i].tau = all_time / (datas[i].acoustic_phonon_scattering_count +
                               datas[i].optical_phonon_scattering_count + 1);
  }
  if (dumping) {
    FILE *fd = fopen(dumpname, "wb");
    fwrite(&n, sizeof(double), 1, fd);
    fwrite(&s, sizeof(double), 1, fd);
    fwrite(dump, sizeof(double), 2 * n * s, fd);
    fclose(fd);
  }
  /* statistics */
  double student_coeff = 3;
  Data m = mean(datas);
  Data sd = stdev(datas) / std::sqrt(n) * student_coeff;

  /* output (stdout) */
  printf("{\n");
  printf("\"v_x\" : [%e, %e],\n", m.v.x, sd.v.x);
  printf("\"v_y\" : [%e, %e],\n", m.v.y, sd.v.y);
  printf("\"power\" : [\n");
  for (int i = 0; i <= number_of_waves; i++) {
    printf("\t[%e, %e]", m.power[i], sd.power[i]);
    putchar(",]"[i == number_of_waves]);
    putchar("\n,"[i == number_of_waves]);
  }
  printf("\n");
  printf("\"tau\" : [%e, %e],\n", m.tau, sd.tau);
  printf("\"acoustic\" : [%d, %d],\n", m.acoustic_phonon_scattering_count,
         sd.acoustic_phonon_scattering_count);
  printf("\"optical\" : [%d, %d],\n", m.optical_phonon_scattering_count,
         sd.optical_phonon_scattering_count);
  printf("\"vertical\" : [%d, %d],\n", m.vertical_transitions_count,
         sd.vertical_transitions_count);
  printf("\"population\" : [\n");
  for (int i = 0; i < mat->bands.size(); i++) {
    printf("\t[%e, %e]", m.population[i], sd.population[i]);
    putchar(",]"[i == mat->bands.size() - 1]);
    printf("\n");
  }
  printf("}\n");

  /* debug info */
  /*
  FILE *f = fopen("tau.dat", "w");
  for (auto & d: datas) {
      fprintf(f, "%e\n", d.tau);
  }
  fclose(f);

  f = fopen("acoustic.dat", "w");
  for (auto & d: datas) {
      fprintf(f, "%d\n", d.acoustic_phonon_scattering_count);
  }
  fclose(f);

  f = fopen("optical.dat", "w");
  for (auto & d: datas) {
      fprintf(f, "%d\n", d.optical_phonon_scattering_count);
  }
  fclose(f);
  */
  return 0;
}
