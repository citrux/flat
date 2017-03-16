#include <cstdio>
#include <ctime>
#include <iostream>
#include <omp.h>
#include "physics.hh"
#include "material.hh"
#include "stats.hh"

int test_prob() {
    auto band = new Band();
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
    const Band lower = Band();
    std::vector<Band> bands = {lower};
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        Particle particle(seeds[i]);
        particle.band = &bands[0];
        particle.p = {0, 0};
        for (int t = 0; t < s; ++t) {
            Vec2 v = particle.band->velocity(particle.p);
            Vec2 f = Ec + E * Vec2(cos(omega * t), cos(omega * t + phi)) + Vec2(v.y, -v.x) * H;
            datas[i].v += v;
            datas[i].power += v.dot(E * Vec2(cos(omega * t), cos(omega * t + phi)));
            particle.p += f;
            for (auto band: bands) {
                if (band.acoustic_phonon_scattering(particle, dt)) {
                    ++datas[i].acoustic_phonon_scattering_count;
                }
                if (band.optical_phonon_scattering(particle, dt)) {
                    ++datas[i].optical_phonon_scattering_count;
                }
            }
        }
        datas[i].v /= s;
        datas[i].power /= s;
        datas[i].tau = all_time / (datas[i].acoustic_phonon_scattering_count + datas[i].optical_phonon_scattering_count + 1);
    }
    Data m = mean(datas);
    Data sd = stdev(datas);
    printf("v_x = %e +/- %e\n", m.v.x, sd.v.x);
    printf("v_y = %e +/- %e\n", m.v.y, sd.v.y);
    printf("power = %e +/- %e\n", m.power, sd.power);
    printf("tau = %e +/- %e\n", m.tau, sd.tau);
    return 0;
}