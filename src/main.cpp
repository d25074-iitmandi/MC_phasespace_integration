#include <iostream>
#include <chrono>
#include "integrator.hpp"
#include "physics.hpp"
#include "utils.hpp"

int main() {
    const int N = 1e6;
    const double sqrt_s = 3.2;
    const double s = sqrt_s * sqrt_s;

    // Serial run
    auto t1 = std::chrono::high_resolution_clock::now();
    Result serial = monte_carlo_serial(N, sqrt_s);
    auto t2 = std::chrono::high_resolution_clock::now();

    double serial_time = std::chrono::duration<double>(t2 - t1).count();

    std::cout << "Serial sigma: " << serial.sigma << "\n";
    std::cout << "Analytic: " << analytic_cross_section(sqrt_s) << "\n";
    std::cout << "Time: " << serial_time << " s\n\n";

    // OpenMP runs
    for (int t = 1; t <= 8; t *= 2) {
        auto start = std::chrono::high_resolution_clock::now();
        Result res = monte_carlo_openmp(N, sqrt_s, t);
        auto end = std::chrono::high_resolution_clock::now();

        double time = std::chrono::duration<double>(end - start).count();
        double speedup = serial_time / time;

        std::cout << "Threads: " << t
                  << " Time: " << time
                  << " Speedup: " << speedup << "\n";

        save_results(t, time, speedup);
    }

    return 0;
}
