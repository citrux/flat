#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>
#include "rng.hh"
#include "physics.hh"

float phi(Rng & rng) {
    float x = 2 * pi * rng.uniform();
    float prob = rng.uniform();
    while(0.5 * (std::cos(x) + 1) < prob) {
        x = 2 * pi * rng.uniform();
        prob = rng.uniform();
    }
    return x;
}

TEST(RngTest, Uniform)
{
    const int n = 100000;
    for (int i = 0; i < 100; ++i) {
        const int m = 100;
        int dist[m] = {0};
        Rng rng(i);
        for (int j = 0; j < n; j++) {
            auto x = rng.uniform();
            ++dist[static_cast<int>(x * m)];
        }
        for (int j = 0; j < m; j++) {
            float p = 1.0 / m;
            float q = 1 - p;
            EXPECT_LE(std::abs(dist[j] - n * p), 5 * std::sqrt(n * p * q)) << i << " " << j << " " << dist[j];
        }
    }
}

TEST(RngTest, Cos) {
    const int n = 100000;
    for (int i = 0; i < 100; ++i) {
        Rng rng(i);
        const int m = 100;
        int dist[m] = {0};
        for (int j = 0; j < n; j++) {
            auto x = phi(rng);
            ++dist[static_cast<int>(x / (2 * pi) * m)];
        }
        for (int j = 0; j < m; j++) {
            float theta = 2 * pi * (j + 0.5) / m;
            float p = (1 + std::cos(theta)) / m;
            float q = 1 - p;
            EXPECT_LE(std::abs(dist[j] - n * p), 5 * std::sqrt(n * p * q)) << i << " " << j << " " << dist[j];
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
