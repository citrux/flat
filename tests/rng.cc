#include <cstdio>
#include <cmath>
#include "rng.hh"
#include "physics.hh"

float phi(Rng & rng) {
    float x = 2 * pi * rng.uniform();
    float prob = rng.uniform();
    while((std::cos(x) + 1) / 2 < prob) {
        x = 2 * pi * rng.uniform();
        prob = rng.uniform();
    }
    return x;
}

int main(int argc, const char *argv[])
{
    int n = 100000;
    for (int i = 0; i < 8; ++i) {
        double mean = 0;
        double stdev = 0;
        Rng rng(i);
        for (int j = 0; j < n; j++) {
            auto x = rng.uniform();
            mean += x;
            stdev += x * x;
        }
        mean /= n;
        stdev /= n;
        stdev -= mean * mean;
        stdev = std::sqrt(stdev);
        printf("seed=%d: mean=%f stdev=%f\n", i, mean, stdev);
    }

    Rng rng(0);
    const int bins = 20;
    float phi_min = 0;
    float phi_max = 1;
    int dist[bins] = {0};
    float mean = 0;
    float stdev = 0;
    for (int i = 0; i < n; i++) {
        auto x = rng.uniform();
        ++dist[(int)(x / (phi_max - phi_min) * bins)];
        mean += x;
        stdev += x * x;
    }
    mean /= n;
    stdev /= n;
    stdev -= mean * mean;
    stdev = std::sqrt(stdev);
    printf("phi: mean=%f stdev=%f\n", mean, stdev);
    FILE *fd = fopen("phi.dat", "w");
    for (int i = 0; i < bins; i++) {
        fprintf(fd, "%f %d\n", phi_min + (phi_max - phi_min) / bins * (i + 0.5), dist[i]);
    }
    fclose(fd);
    return 0;
}
