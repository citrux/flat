#include <cstdio>
#include <ctime>
#include <iostream>
#include <omp.h>
#include "physics.hh"
#include "material.hh"
#include "stats.hh"

int test_prob() {
    auto band = new BigrapheneLower(300);
    const auto table = band->table;
    for (auto entry: table) {
        float e = entry.energy;
        float e1 = e - band->optical_phonon_energy;
        float ac_prob = 0;
        for (auto i: entry.integrals) {
            ac_prob += band->acoustic_phonon_constant * i.integral;
        }
        float op_prob = 0;
        int i = (e1 - table[0].energy) / (table[1].energy - table[0].energy);
        if (i >= 0) {
            auto entry1 = table[i];
            for (auto j: entry1.integrals) {
                op_prob += band->optical_phonon_constant * j.integral;
            }
        }
        printf("%e %e %e %e\n", e, ac_prob, op_prob, ac_prob + op_prob);
    }
    return 0;
}

int main(int argc, char const *argv[])
{
    Vec2 E, Ec;
    float H, omega, phi, T;
    int n;
    float dt, all_time;
    std::cin >> Ec.x >> Ec.y;
    std::cin >> H;
    std::cin >> E.x >> E.y;
    std::cin >> omega;
    std::cin >> phi;
    std::cin >> T;
    std::cin >> n;
    std::cin >> dt;
    std::cin >> all_time;

    printf("Field configuration:\n");
    printf("Ec: {%e, %e}\n", Ec.x, Ec.y);
    printf("H: %e\n", H);
    printf("E: {%e, %e}\n", E.x, E.y);

    const float field_dimensionless_factor = e * v_f * dt / eV;
    Ec *= field_dimensionless_factor;
    E *= field_dimensionless_factor;
    H *= v_f / c * field_dimensionless_factor;
    omega *= dt;

    int s = all_time / dt + 1;

    std::vector<Data> datas(n);
    std::vector<unsigned int> seeds(n);
    srand(time(nullptr));
    for (int i = 0; i < n; ++i) {
        seeds[i] = rand();
    }
    BigrapheneLower *lower = new BigrapheneLower(T);
    std::vector<Band*> bands = {lower};
    puts("start calculation");
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        Particle particle(seeds[i]);
        particle.band = bands[0];
        particle.p = {0, 0};
        float p = 0;
        float r1, r2;
        do {
            r1 = particle.rng.uniform();
            r2 = particle.rng.uniform();
            p = -log(r1);
        } while (exp(- (particle.band->energy(p) - particle.band->min_energy()) / k / T) < r2);
        float theta = 2 * M_PI * particle.rng.uniform();
        particle.p = {p * std::cos(theta), p * std::sin(theta)};
        for (int t = 0; t < s; ++t) {
            Vec2 v = particle.band->velocity(particle.p);
            Vec2 f = Ec + E * Vec2(cos(omega * t), cos(omega * t + phi)) + Vec2(v.y, -v.x) * H;
            datas[i].v += v;
            datas[i].power += v.dot(E * Vec2(cos(omega * t), cos(omega * t + phi)));
            particle.p += f;
            float wsum = 0;
            for (auto & band: bands) {
                for (auto & result: band->acoustic_phonon_scattering(particle))
                    wsum += result.rate;
                for (auto & result: band->optical_phonon_scattering(particle))
                    wsum += result.rate;
            }
            particle.r -= wsum * dt;
            if (particle.r < 0) {
                float w = particle.rng.uniform() * wsum;
                for (auto & band: bands) {       
                    for (auto & result: band->acoustic_phonon_scattering(particle)) {
                        w -= result.rate;
                        if (w < 0) {
                            ++datas[i].acoustic_phonon_scattering_count;
                            particle.reset_r();
                            particle.band = band;
                            goto end;
                        }
                    }
                    for (auto & result: band->optical_phonon_scattering(particle)) {
                        w -= result.rate;
                        if (w < 0) {
                            ++datas[i].optical_phonon_scattering_count;
                            particle.reset_r();
                            particle.band = band;
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
    float student_coeff = 3;
    Data m = mean(datas);
    Data sd = stdev(datas);
    printf("v_x = %e +/- %e\n", m.v.x, sd.v.x);
    printf("v_y = %e +/- %e\n", m.v.y, sd.v.y);
    printf("power = %e +/- %e\n", m.power, sd.power);
    printf("tau = %e +/- %e\n", m.tau, sd.tau);
    printf("acoustic = %d +/- %d\n", m.acoustic_phonon_scattering_count, sd.acoustic_phonon_scattering_count);
    printf("optical = %d +/- %d\n", m.optical_phonon_scattering_count, sd.optical_phonon_scattering_count);
    return 0;
}