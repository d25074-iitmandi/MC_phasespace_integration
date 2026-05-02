#include "rng.hpp"

RNG::RNG(int seed) : gen(seed), dist(0.0, 1.0) {}
//Distribution of particles
double RNG::uniform() {
    return dist(gen);
}
