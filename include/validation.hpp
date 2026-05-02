#pragma once
#include <vector>
#include <string>

struct BinResult {
    double costheta_center;   // bin centre
    double mc_dsigma;         // MC dσ/d(cosθ) for this bin
    double mc_error;          // statistical uncertainty on mc_dsigma
    double analytic_dsigma;   // analytic dσ/d(cosθ) at bin centre
};

struct ChiSquareResult {
    double chi2;              // raw χ²
    int    dof;               // degrees of freedom (n_bins - 1)
    double chi2_per_dof;      // χ²/dof — should be ~1.0
    double p_value;           // probability of observing this χ² or worse
};

// Returns per-bin results for the angular distribution
std::vector<BinResult> angular_distribution(int N, double sqrt_s,
                                            int n_bins = 20);

// Computes χ² from the angular distribution
ChiSquareResult chi_square_test(const std::vector<BinResult>& bins);

// Prints a formatted report of both results
void print_angular_report(const std::vector<BinResult>& bins,
                          const ChiSquareResult& chi2);