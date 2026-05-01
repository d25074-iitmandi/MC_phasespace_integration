#pragma once
#include <random>

struct RNG {
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;

    RNG(int seed);
    double uniform();
};
