#include "integrator.h"
#include "physics.h"
#include "rng.h"
#include <cmath>
#include <omp.h>

Result monte_carlo_serial(int N, double s) {
    RNG rng(42);

    double sum_w = 0.0, sum_w2 = 0.0; //Weights

    for (int i = 0; i < N; i++) {
        double costh = 2.0 * rng.uniform() - 1.0;
        double M2 = matrix_element(costh);
        double w = M2 / (32.0 * M_PI * M_PI * s);

        sum_w += w;
        sum_w2 += w * w;
    }

    double volume = 4.0 * M_PI;
    double sigma = volume * sum_w / N;
    double var = (sum_w2 / N - pow(sum_w / N, 2)) / N;

    return {sigma, volume * sqrt(var)};
}

Result monte_carlo_openmp(int N, double s, int nthreads) {
    double sum_w = 0.0, sum_w2 = 0.0;

    omp_set_num_threads(nthreads);

    #pragma omp parallel
    {
        RNG rng(42 + omp_get_thread_num());

        #pragma omp for reduction(+:sum_w,sum_w2)
        for (int i = 0; i < N; i++) {
            double costh = 2.0 * rng.uniform() - 1.0;
            double M2 = matrix_element(costh);
            double w = M2 / (32.0 * M_PI * M_PI * s);

            sum_w += w;
            sum_w2 += w * w;
        }
    }

    double volume = 4.0 * M_PI;
    double sigma = volume * sum_w / N;
    double var = (sum_w2 / N - pow(sum_w / N, 2)) / N;

    return {sigma, volume * sqrt(var)};
}
