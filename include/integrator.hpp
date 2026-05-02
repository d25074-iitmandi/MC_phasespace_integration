#pragma once

struct Result {
    double sigma;
    double error;
};

Result monte_carlo_serial(int N, double s);
Result monte_carlo_openmp(int N, double s, int nthreads);
