#include "rng.h"

RNG::RNG(int seed) : gen(seed), dist(0.0, 1.0) {}

double RNG::uniform() {
    return dist(gen);
}
