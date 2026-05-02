#include <iostream>
#include <chrono>
#include "integrator.hpp"
#include "physics.hpp"
#include "utils.hpp"
#include <cmath>
#include "validation.hpp"
int main() {
    //const int N = 1e8;
    const double sqrt_s = 3.2;
    const double s = sqrt_s * sqrt_s;
    double analytic = analytic_cross_section(sqrt_s);\
    double prev_error = 0.0;
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

    std::cout<<"\n=============Serial Implementation begins=============\n";
    std::cout<< "For N= 10^"<< logN <<"\n";     
    std::cout << "Sigma_MC: " << serial.sigma << " GeV^-2\n";
    std::cout << "Analytic: " << analytic << " GeV^-2\n";
    std::cout << "Absolute Error: " << abs_err << " GeV^-2\n";
    std::cout << "MC Error "<< serial.error <<"\n";
    std::cout << "Relative error = " << rel_err * 100.0 << " %\n";
    std::cout << "Error ratio    = " << ratio           << "\n";
    std::cout << "Time           = " << serial_time         << " s\n\n";

    prev_error = serial.error;
    save_convergence_results(logN, serial.sigma, analytic, ratio, serial.error );

    // OpenMP Implementation
    std::cout<< "\n=======OpenMp Implementation begins=======\n";
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

    
}
// --- In main(), after your convergence test ---

std::cout << "\n=== Validation: N=1e6, sqrt_s=" << sqrt_s << " GeV ===\n";

const int N_val  = 1e6;
const int N_BINS = 20;

auto bins    = angular_distribution(N_val, sqrt_s, N_BINS);
auto chi2res = chi_square_test(bins);
print_angular_report(bins, chi2res);
    return 0;
}
