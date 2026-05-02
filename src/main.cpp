#include <iostream>
#include <chrono>
#include "integrator.hpp"
#include "physics.hpp"
#include "utils.hpp"
#include <cmath>
#include "asymmetry.hpp"
#include <fstream>

int main() {
    //const int N = 1e8;
    const double sqrt_s = 0.22;
    const double s = sqrt_s * sqrt_s;
    double analytic = analytic_cross_section(sqrt_s);\
    double prev_error = 0.0;
    //Performance Analysis
    std::ofstream perf_file("results/performance.csv");
    perf_file << "logN,N,threads,time,speedup,efficiency,sigma,error\n";
    // Serial implementation
    for (int logN =3; logN <= 7; logN++){
    auto t1 = std::chrono::high_resolution_clock::now();
    int N = pow(10,logN);
    Result serial = monte_carlo_serial(N, sqrt_s);
    auto t2 = std::chrono::high_resolution_clock::now();

    double serial_time = std::chrono::duration<double>(t2 - t1).count();

    double abs_err = std::abs(serial.sigma - analytic);
    double rel_err = abs_err/analytic;
    double ratio    = (logN > 3) ? prev_error / serial.error : 0.0; // should → sqrt(10) ≈ 3.162

std::cout << "\n=== Serial  N=10^" << logN << " ===\n";
        std::cout << "Sigma_MC : " << serial.sigma  << " GeV^-2\n";
        std::cout << "Analytic : " << analytic       << " GeV^-2\n";
        std::cout << "Abs err  : " << abs_err        << " GeV^-2\n";
        std::cout << "MC error : " << serial.error   << "\n";
        std::cout << "Rel err  : " << rel_err        << " %\n";
        std::cout << "Ratio    : " << ratio          << "  (expect 3.162)\n";
        std::cout << "Time     : " << serial_time    << " s\n";


    // Serial counts as thread=1 in perf file
        perf_file << logN        << ","
                  << N           << ","
                  << 1           << ","
                  << serial_time << ","
                  << 1.0         << ","   // speedup = 1 for serial
                  << 100.0       << ","   // efficiency = 100%
                  << serial.sigma << ","
                  << serial.error << "\n";
    save_convergence_results(logN, N, serial.sigma, analytic, serial.error);

    // OpenMP Implementation
     std::cout << "\n--- OpenMP  N=10^" << logN << " ---\n";
        for (int t = 2; t <= 8; t *= 2) {
            auto start = std::chrono::high_resolution_clock::now();
            Result res  = monte_carlo_openmp(N, sqrt_s, t);
            auto end    = std::chrono::high_resolution_clock::now();
            double time_omp = std::chrono::duration<double>(end - start).count();

            double speedup    = serial_time / time_omp;
            double efficiency = speedup / t * 100.0;

            std::cout << "Threads: " << t
                      << "  Time: "      << time_omp
                      << "  Speedup: "   << speedup
                      << "  Efficiency: "<< efficiency << "%\n";

            perf_file << logN        << ","
                      << N           << ","
                      << t           << ","
                      << time_omp    << ","
                      << speedup     << ","
                      << efficiency  << ","
                      << res.sigma   << ","
                      << res.error   << "\n";
        }

        prev_error = serial.error;
    
    

    
}


std::cout << "\n=================================================\n";
std::cout << "   SECTION: Forward-Backward Asymmetry Check\n";
std::cout << "=================================================\n";

// 1. Detailed single-energy check at your primary working point
HemisphereResult afb = forward_backward_asymmetry(1000000, sqrt_s);
print_asymmetry_report(afb, sqrt_s);

// 2. Scan across energies — use the same grid as the energy scan
//    A_FB must be zero at every single energy point
std::vector<double> scan_energies = {
    0.215, 0.230, 0.250, 0.300, 0.500,
    1.0,   2.0,   3.2,   5.0,  10.0
};
auto afb_scan = asymmetry_energy_scan(200000, scan_energies);
print_asymmetry_scan_report(afb_scan);
    return 0;
}
