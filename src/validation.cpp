#include "validation.hpp"
#include "physics.hpp"
#include "rng.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>

// ---------------------------------------------------------------------------
// Angular distribution: dσ/d(cosθ)
// ---------------------------------------------------------------------------
// Strategy:
//   Sample cosθ uniformly over [-1,1] as usual, but NOW record which bin
//   each event falls into and accumulate sum_w and sum_w2 per bin.
//   After N events, the MC estimate of dσ/d(cosθ) in bin i is:
//
//       (dσ/d(cosθ))_i = (4π / Δ(cosθ)) * sum_w[i] / N
//
//   where 4π is the full solid angle volume and Δ(cosθ) = 2/n_bins is the
//   bin width. The 1/Δ(cosθ) converts from total-σ units to differential.
// ---------------------------------------------------------------------------

std::vector<BinResult> angular_distribution(int N, double sqrt_s, int n_bins) {
    double s        = sqrt_s * sqrt_s;
    double mu_mass = 0.105;
    double alpha = 1.0/137.0;
    double beta     = std::sqrt(1.0 - 4.0 * mu_mass * mu_mass / s);
    double bin_width = 2.0 / n_bins;          // cosθ range = [-1,1], width = 2
    double volume   = 4.0 * M_PI;             // azimuthal integral already done

    // Per-bin accumulators
    std::vector<double> sum_w (n_bins, 0.0);
    std::vector<double> sum_w2(n_bins, 0.0);

    RNG rng(42);
    for (int i = 0; i < N; i++) {
        double costh = 2.0 * rng.uniform() - 1.0;   // uniform in [-1,1]

        // Bin index: costh=-1 → bin 0, costh=+1 → bin n_bins-1
        int bin = static_cast<int>((costh + 1.0) / bin_width);
        if (bin >= n_bins) bin = n_bins - 1;         // guard upper edge

        double M2 = matrix_element(costh, sqrt_s);
        // Phase-space weight including β factor (2-body LIPS)
        double w  = beta * M2 / (32.0 * M_PI * M_PI * s);

        sum_w [bin] += w;
        sum_w2[bin] += w * w;
    }

    // Convert accumulators to dσ/d(cosθ) per bin
    std::vector<BinResult> results(n_bins);
    for (int b = 0; b < n_bins; b++) {
        double costh_center = -1.0 + (b + 0.5) * bin_width;

        // MC differential cross-section for this bin:
        //   dσ/d(cosθ) = volume * <w>_bin / bin_width
        //   where <w>_bin = sum_w[b] / N   (average over ALL events, not
        //   just those in the bin — uniform sampling means each bin
        //   contributes fraction bin_width/2 of all events on average)
        double dsigma_mc  = volume * sum_w[b] / N / bin_width;

        // Statistical variance of the mean, then propagated through the
        // same volume/N/bin_width prefactor
        double var_mean   = (sum_w2[b] / N - std::pow(sum_w[b] / N, 2)) / N;
        double error_mc   = volume * std::sqrt(var_mean) / bin_width;

        // Analytic dσ/d(cosθ) = (α²β/4s) * [(1+4m²/s) + (1-4m²/s)cos²θ]
        // This is the full QED result without the 2π from φ integration
        // (we keep it differential in cosθ only, integrating out φ gives 2π,
        //  so the total is 2π × this expression, and ∫ this dcosθ = analytic σ)
        double r          = 4.0 * mu_mass * mu_mass / s;
        double dsigma_an  = (alpha * alpha * beta / (4.0 * s))
                          * ((1.0 + r) + (1.0 - r) * costh_center * costh_center);
        // Multiply by 2π to match the full azimuthal-integrated differential
        dsigma_an *= 2.0 * M_PI;

        results[b] = { costh_center, dsigma_mc, error_mc, dsigma_an };
    }
    return results;
}

// ---------------------------------------------------------------------------
// χ² goodness-of-fit
// ---------------------------------------------------------------------------
// χ² = Σ_i  (MC_i - analytic_i)² / σ²_i
//
// where σ_i = mc_error[i] is the MC statistical uncertainty on bin i.
// dof = n_bins - 1  (one constraint: total normalisation is shared)
//
// p-value from the chi2 CDF using the regularised incomplete gamma function
// p = 1 - P(dof/2, χ²/2) = Q(dof/2, χ²/2)
// Computed via the continued-fraction approximation (no external library).
// ---------------------------------------------------------------------------

// Regularised upper incomplete gamma Q(a,x) via continued fractions
// Sufficient precision for p-value reporting; avoids needing <boost> or GSL
static double gammq(double a, double x) {
    if (x <= 0.0) return 1.0;

    double lnGammaA = std::lgamma(a);   // drop the custom implementation entirely

    const double FPMIN = 1.0e-300;
    const int    ITMAX = 200;
    const double EPS   = 3.0e-7;

    double b    = x + 1.0 - a;
    double c_cf = 1.0 / FPMIN;
    double d    = 1.0 / b;
    double h    = d;
    for (int i = 1; i <= ITMAX; i++) {
        double an = -i * (i - a);
        b  += 2.0;
        d   = an * d + b;
        if (std::abs(d)    < FPMIN) d    = FPMIN;
        c_cf = b + an / c_cf;
        if (std::abs(c_cf) < FPMIN) c_cf = FPMIN;
        d = 1.0 / d;
        double del = d * c_cf;
        h *= del;
        if (std::abs(del - 1.0) < EPS) break;
    }
    return std::exp(-x + a * std::log(x) - lnGammaA) * h;
}

ChiSquareResult chi_square_test(const std::vector<BinResult>& bins) {
    int n = static_cast<int>(bins.size());
    double chi2 = 0.0;

    for (const auto& b : bins) {
        if (b.mc_error <= 0.0) continue;   // skip empty bins
        double residual = b.mc_dsigma - b.analytic_dsigma;
        chi2 += (residual * residual) / (b.mc_error * b.mc_error);
    }

    int dof     = n - 1;      // -1 for shared normalisation constraint
    double p    = gammq(0.5 * dof, 0.5 * chi2);   // p-value = Q(dof/2, χ²/2)

    return { chi2, dof, chi2 / dof, p };
}

// ---------------------------------------------------------------------------
// Formatted report
// ---------------------------------------------------------------------------

void print_angular_report(const std::vector<BinResult>& bins,
                          const ChiSquareResult& chi2) {
    using std::setw;
    using std::setprecision;

    std::cout << "\n=== Angular distribution: dσ/d(cosθ) ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << setw(12) << "cosθ_centre"
              << setw(16) << "MC dσ/dcosθ"
              << setw(16) << "MC error"
              << setw(16) << "Analytic"
              << setw(12) << "pull\n";
    std::cout << std::string(72, '-') << "\n";

    for (const auto& b : bins) {
        double pull = (b.mc_error > 0)
                    ? (b.mc_dsigma - b.analytic_dsigma) / b.mc_error
                    : 0.0;
        std::cout << std::setw(12) << std::setprecision(3) << b.costheta_center
                  << std::setw(16) << std::setprecision(5) << std::scientific << b.mc_dsigma
                  << std::setw(16) << b.mc_error
                  << std::setw(16) << b.analytic_dsigma
                  << std::setw(12) << std::fixed << std::setprecision(2) << pull
                  << "\n";
    }

    std::cout << "\n=== Chi-square goodness of fit ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  χ²          = " << chi2.chi2         << "\n";
    std::cout << "  dof         = " << chi2.dof           << "\n";
    std::cout << "  χ²/dof      = " << chi2.chi2_per_dof  << "  (ideal: ~1.0)\n";
    std::cout << "  p-value     = " << chi2.p_value       << "\n";

    // Interpret the result
    std::cout << "\n  Interpretation: ";
    if (chi2.p_value > 0.05)
        std::cout << "PASS — MC consistent with analytic at 95% CL\n";
    else if (chi2.p_value > 0.01)
        std::cout << "MARGINAL — tension at 95% CL, acceptable at 99%\n";
    else
        std::cout << "FAIL — significant disagreement (p < 0.01)\n";

    // Pull distribution summary
    double pull_mean = 0.0, pull_rms = 0.0;
    for (const auto& b : bins) {
        double p = (b.mc_error > 0)
                 ? (b.mc_dsigma - b.analytic_dsigma) / b.mc_error : 0.0;
        pull_mean += p;
        pull_rms  += p * p;
    }
    int n = static_cast<int>(bins.size());
    pull_mean /= n;
    pull_rms   = std::sqrt(pull_rms / n);
    std::cout << "  Pull mean   = " << std::setprecision(3) << pull_mean
              << "  (ideal: 0.0)\n";
    std::cout << "  Pull RMS    = " << pull_rms
              << "  (ideal: 1.0)\n";
}