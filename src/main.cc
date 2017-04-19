#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include "physics.hh"
#include "material.hh"
#include "stats.hh"

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

template<typename Out>
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

struct Wave {
    Vec2 E;
    float H;
    float omega;
    float phi;
};

int main(int argc, char const *argv[])
{
    bool dumping = false;
    char dumpname[256] = "dump.bin";
    for(int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-d")) {
            dumping = true;
            if (i + 1 < argc)
                strcpy(dumpname, argv[i+1]);
        }
    }
    std::string material, material_name;
    Vec2 E, Ec;
    float H, Hc, omega, phi, T;
    int n;
    int number_of_waves;
    float dt, all_time;

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

    printf("threads: %d\n", omp_get_num_threads());
    if (material_name != "bigraphene" && material_name != "graphene") {
        puts("incorrect material");
        exit(1);
    }
    printf("Field configuration:\n");
    printf("Ec: {%e, %e}\n", Ec.x, Ec.y);
    printf("Hc: %e\n", Hc);
    printf("E: {%e, %e}\n", E.x, E.y);
    printf("H: %e\n", H);

    const float field_dimensionless_factor = e * v_f * dt / eV;
    Ec *= field_dimensionless_factor;
    Hc *= v_f / c * field_dimensionless_factor;
    for (auto & w: waves) {
        w.E *= field_dimensionless_factor;
        w.H *= v_f / c * field_dimensionless_factor;
        w.omega *= dt;
    }

    int s = all_time / dt;

    std::vector<Data> datas(n);
    std::vector<unsigned int> seeds(n);
    float* dump;
    if (dumping) { dump = new float[2 * n * s]; }
    srand(time(nullptr));
    for (int i = 0; i < n; ++i) {
        seeds[i] = rand();
    }
    std::vector<Band*> bands;
    if (material_name == "bigraphene") {
        float delta = 0;
        if (material_params.size() > 1) {
            delta = atof(material_params[1].c_str());
        }
        Bigraphene::Lower *lower = new Bigraphene::Lower(T, delta);
        bands = {lower};
        if (material_params.size() > 2) {
            Bigraphene::Upper *upper = new Bigraphene::Upper(T);
            bands.push_back(upper);
        }
    } else {
        float delta = 0;
        if (material_params.size() > 1) {
            delta = atof(material_params[1].c_str());
        }
        Graphene::Bnd *band = new Graphene::Bnd(T, delta);
        bands = {band};
    }
    puts("start calculation");
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        Particle particle(seeds[i]);
        particle.band = bands[0];
        particle.p = {0, 0};
        /* Boltzmann-distributed initial condition */
        float prob;
        do {
            float p = particle.rng.uniform();
            prob = particle.rng.uniform();
            float theta = 2 * pi * particle.rng.uniform();
            particle.p = {p * std::cos(theta), p * std::sin(theta)};
        } while (exp(-(particle.band->energy(particle.p) - particle.band->min_energy()) / k / T) < prob);
        datas[i].power.assign(number_of_waves + 1, 0);
        datas[i].population.assign(bands.size(), 0);

        /* simulation */
        for (int t = 0; t < s; ++t) {
            if (dumping) {
                dump[2 * t * n + 2 * i] = particle.p.x;
                dump[2 * t * n + 2 * i + 1] = particle.p.y;
            }
            Vec2 v = particle.band->velocity(particle.p);
            Vec2 f = Ec + Vec2(v.y, -v.x) * Hc;
            for (auto w: waves) {
                f += w.E * Vec2(std::cos(w.omega * t),
                                std::cos(w.omega * t + w.phi)) +
                     Vec2(v.y, -v.x) * w.H * std::cos(w.omega * t);
            }
            if (std::isnan(f.x) || std::isnan(f.y)) {
                puts("force is nan!");
                printf("%e %e\n", f.x, f.y);
                printf("%e %e\n", v.x, v.y);
                printf("%e %e\n", particle.p.x, particle.p.y);
                printf("%e\n", particle.band->energy(particle.p));
                exit(1);
            }
            if (std::isnan(v.x) || std::isnan(v.y)) {
                puts("wtf");
                printf("%e %e\n", v.x, v.y);
                printf("%e %e\n", particle.p.x, particle.p.y);
                exit(1);
            }
            int band_index = std::find(bands.begin(), bands.end(), particle.band) - bands.begin();
            datas[i].population[band_index] += 1;
            datas[i].v += v;
            datas[i].power[0] += v.dot(Ec);
            for (int j = 0; j < number_of_waves; j++) {
                const float Ex = waves[j].E.x;
                const float Ey = waves[j].E.y;
                const float omega = waves[j].omega;
                const float phi = waves[j].phi;
                Vec2 E = { Ex*std::cos(omega*t), Ey*std::cos(omega*t+phi) };
                datas[i].power[j+1] += v.dot(E);
            }
            particle.p += f;
            float wsum = 0;
            for (auto const band: bands) {
                for (auto const & result: band->acoustic_phonon_scattering(particle))
                    wsum += result.rate;
                for (auto const & result: band->optical_phonon_scattering(particle))
                    wsum += result.rate;
            }
            particle.r -= wsum * dt;
            if (particle.r < 0) {
                float w = particle.rng.uniform() * wsum;
                for (auto const band: bands) {
                    for (auto const & result: band->acoustic_phonon_scattering(particle)) {
                        w -= result.rate;
                        if (w < 0) {
                            ++datas[i].acoustic_phonon_scattering_count;
                            particle.reset_r();
                            particle.band = band;
                            particle.p = band->momentum_scattering(result.momentum, particle);
                            goto end;
                        }
                    }
                    for (auto const & result: band->optical_phonon_scattering(particle)) {
                        w -= result.rate;
                        if (w < 0) {
                            ++datas[i].optical_phonon_scattering_count;
                            particle.reset_r();
                            particle.band = band;
                            particle.p = band->momentum_scattering(result.momentum, particle);
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
        datas[i].tau = all_time / (datas[i].acoustic_phonon_scattering_count + datas[i].optical_phonon_scattering_count + 1);
    }
    if (dumping) {
        FILE *fd = fopen(dumpname, "wb");
        fwrite(&n, sizeof(float), 1, fd);
        fwrite(&s, sizeof(float), 1, fd);
        fwrite(dump, sizeof(float), 2 * n * s, fd);
        fclose(fd);
    }
    /* statistics */
    float student_coeff = 3;
    Data m = mean(datas);
    Data sd = stdev(datas) / (float)std::sqrt(n) * student_coeff;

    /* output (stdout) */
    printf("v_x = %e +/- %e\n", m.v.x, sd.v.x);
    printf("v_y = %e +/- %e\n", m.v.y, sd.v.y);
    printf("power:\n\tstatic = %e +/- %e\n", m.power[0], sd.power[0]);
    for (int i = 0; i < number_of_waves; i++) {
        printf("\twave %d = %e +/- %e\n", i, m.power[i+1], sd.power[i+1]);
    }
    printf("tau = %e +/- %e\n", m.tau, sd.tau);
    printf("acoustic = %d +/- %d\n", m.acoustic_phonon_scattering_count, sd.acoustic_phonon_scattering_count);
    printf("optical = %d +/- %d\n", m.optical_phonon_scattering_count, sd.optical_phonon_scattering_count);
    return 0;
}
