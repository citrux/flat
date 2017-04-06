#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <omp.h>
#include "physics.hh"
#include "material.hh"
#include "stats.hh"

int main(int argc, char const *argv[])
{
    bool dumping = false;
    for(int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-d")) {
            dumping = true;
        }
    }
    std::string material;
    Vec2 E, Ec;
    float H, Hc, omega, phi, T;
    int n;
    float dt, all_time;

    /* input */
    std::cin >> material;
    std::cin >> Ec.x >> Ec.y;
    std::cin >> Hc;
    std::cin >> E.x >> E.y;
    std::cin >> H;
    std::cin >> omega;
    std::cin >> phi;
    std::cin >> T;
    std::cin >> n;
    std::cin >> dt;
    std::cin >> all_time;

    if (material != "bigraphene" && material != "graphene") {
        puts("incorrect material");
        exit(1);
    }
    printf("Field configuration:\n");
    printf("Ec: {%e, %e}\n", Ec.x, Ec.y);
    printf("H: %e\n", H);
    printf("E: {%e, %e}\n", E.x, E.y);

    const float field_dimensionless_factor = e * v_f * dt / eV;
    Ec *= field_dimensionless_factor;
    E *= field_dimensionless_factor;
    Hc *= v_f / c * field_dimensionless_factor;
    H *= v_f / c * field_dimensionless_factor;
    omega *= dt;

    int s = all_time / dt;

    std::vector<Data> datas(n);
    std::vector<unsigned int> seeds(n);
    /*
    float* dump = new float[2 * n * s];
    */
    srand(time(nullptr));
    for (int i = 0; i < n; ++i) {
        seeds[i] = rand();
    }
    std::vector<Band*> bands;
    if (material == "bigraphene") {
        Bigraphene::Lower *lower = new Bigraphene::Lower(T);
        //Bigraphene::Upper *upper = new Bigraphene::Upper(T);
        bands = {lower/*, upper*/};
    } else {
        Graphene::Bnd *band = new Graphene::Bnd(T);
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

        /* simulation */
        for (int t = 0; t < s; ++t) {
            /*
            dump[2 * t * n + 2 * i] = particle.p.x;
            dump[2 * t * n + 2 * i + 1] = particle.p.y;
            */
            Vec2 v = particle.band->velocity(particle.p);
            Vec2 f = Ec + E * Vec2(std::cos(omega * t), std::cos(omega * t + phi)) +
                Vec2(v.y, -v.x) * (Hc + H * std::cos(omega * t));
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
            datas[i].v += v;
            datas[i].power += v.dot(E * Vec2(cos(omega * t), cos(omega * t + phi)));
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
    /*
    if (dumping) {
        FILE *fd = fopen("dump.bin", "wb");
        fwrite(&n, sizeof(float), 1, fd);
        fwrite(&s, sizeof(float), 1, fd);
        fwrite(dump, sizeof(float), 2 * n * s, fd);
        fclose(fd);
    }
    */
    /* statistics */
    float student_coeff = 3;
    Data m = mean(datas);
    Data sd = stdev(datas) / (float)std::sqrt(n) * student_coeff;

    /* output (stdout) */
    printf("v_x = %e +/- %e\n", m.v.x, sd.v.x);
    printf("v_y = %e +/- %e\n", m.v.y, sd.v.y);
    printf("power = %e +/- %e\n", m.power, sd.power);
    printf("tau = %e +/- %e\n", m.tau, sd.tau);
    printf("acoustic = %d +/- %d\n", m.acoustic_phonon_scattering_count, sd.acoustic_phonon_scattering_count);
    printf("optical = %d +/- %d\n", m.optical_phonon_scattering_count, sd.optical_phonon_scattering_count);
    return 0;
}
