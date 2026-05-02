#include "asymmetry.hpp"
#include "physics.hpp"
#include "rng.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utils.hpp>
// ---------------------------------------------------------------------------
// Core hemisphere integrator
// ---------------------------------------------------------------------------
// We run a SINGLE loop and split events into forward (cosθ > 0) and backward
// (cosθ < 0) accumulators simultaneously. This is important: using two
// separate MC runs would introduce independent statistical fluctuations in σ_F
// and σ_B that are uncorrelated, inflating A_FB errors. Splitting one run
// gives negatively correlated estimators and a tighter uncertainty on A_FB.
// ---------------------------------------------------------------------------

HemisphereResult forward_backward_asymmetry(int N, double sqrt_s) {
    double s    = sqrt_s * sqrt_s;
    double mu_mass = 0.1056;
    double beta = std::sqrt(1.0 - 4.0 * mu_mass * mu_mass / s);

    // Forward (F) and backward (B) accumulators
    double sum_wF = 0.0, sum_w2F = 0.0;
    double sum_wB = 0.0, sum_w2B = 0.0;

    RNG rng(42);
    for (int i = 0; i < N; i++) {
        double costh = 2.0 * rng.uniform() - 1.0;
        double M2    = matrix_element(costh, sqrt_s);
        double w     = beta * M2 / (32.0 * M_PI * M_PI * s);

        if (costh > 0.0) {
            sum_wF  += w;
            sum_w2F += w * w;
        } else {
            sum_wB  += w;
            sum_w2B += w * w;
        }
    }

    // Each hemisphere samples cosθ over half the range [-1,0] or [0,1]
    // but the sampling is still uniform over the full [-1,1] range, so
    // the correct estimator is:
    //
    //   σ_F = (4π) * sum_wF / N       (not sum_wF / (N/2))
    //
    // This is because sum_wF already has roughly N/2 nonzero contributions,
    // and dividing by N (not N/2) accounts for the fact that we sampled the
    // full interval but only accumulated F events. The volume factor 4π
    // covers the full solid angle; the hemisphere split emerges naturally.

    double volume = 4.0 * M_PI;

    double sigma_F = volume * sum_wF / N;
    double sigma_B = volume * sum_wB / N;

    // Variance on each hemisphere estimate
    // Var(σ_F) = V² * (⟨w²⟩_F - ⟨w⟩_F²) / N
    // where ⟨w⟩_F = sum_wF/N  (average over ALL N events, not just F ones)
    double var_F = (sum_w2F / N - std::pow(sum_wF / N, 2)) / N;
    double var_B = (sum_w2B / N - std::pow(sum_wB / N, 2)) / N;

    double err_F = volume * std::sqrt(var_F);
    double err_B = volume * std::sqrt(var_B);

    // A_FB = (σ_F - σ_B) / (σ_F + σ_B)
    double denom   = sigma_F + sigma_B;
    double A_FB    = (denom > 0) ? (sigma_F - sigma_B) / denom : 0.0;

    // Error propagation:
    // ∂A_FB/∂σ_F =  2σ_B / (σ_F+σ_B)²
    // ∂A_FB/∂σ_B = -2σ_F / (σ_F+σ_B)²
    // σ²(A_FB) = (∂A/∂σ_F)² σ²_F + (∂A/∂σ_B)² σ²_B
    double dAdF    = (denom > 0) ?  2.0 * sigma_B / (denom * denom) : 0.0;
    double dAdB    = (denom > 0) ? -2.0 * sigma_F / (denom * denom) : 0.0;
    double A_FB_err = std::sqrt(dAdF*dAdF * var_F * volume*volume
                              + dAdB*dAdB * var_B * volume*volume);

    return { sigma_F, sigma_B, err_F, err_B,
             A_FB, A_FB_err, sigma_F + sigma_B };
}

// ---------------------------------------------------------------------------
// Energy scan of A_FB
// ---------------------------------------------------------------------------

std::vector<AsymmetryScanPoint> asymmetry_energy_scan(
        int N_per_point,
        const std::vector<double>& energies) {

    std::vector<AsymmetryScanPoint> results;
    results.reserve(energies.size());
    for (double sq : energies)
        results.push_back({ sq, forward_backward_asymmetry(N_per_point, sq) });
    return results;
}

// ---------------------------------------------------------------------------
// Detailed single-energy report
// ---------------------------------------------------------------------------

void print_asymmetry_report(const HemisphereResult& r, double sqrt_s) {
    double s       = sqrt_s * sqrt_s;
    double mu_mass = 0.1056;
    double beta    = std::sqrt(1.0 - 4.0 * mu_mass * mu_mass / s);

    // Significance: how many σ is A_FB away from 0?
    double significance = (r.A_FB_error > 0)
                        ? std::abs(r.A_FB) / r.A_FB_error : 0.0;

    std::cout << "\n=== Forward-Backward Asymmetry ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  sqrt_s          = " << sqrt_s       << " GeV\n";
    std::cout << "  beta            = " << beta          << "\n\n";

    std::cout << std::scientific << std::setprecision(5);
    std::cout << "  sigma_F         = " << r.sigma_F
              << " ± " << r.error_F    << " GeV^-2\n";
    std::cout << "  sigma_B         = " << r.sigma_B
              << " ± " << r.error_B    << " GeV^-2\n";
    std::cout << "  sigma_F+sigma_B = " << r.sigma_total << " GeV^-2\n\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  A_FB            = " << r.A_FB
              << " ± " << r.A_FB_error << "\n";
    std::cout << "  Significance    = " << significance
              << " σ  (expect < 2σ)\n";
    std::cout << "  Expected (QED)  = 0.000000  (exact at tree level)\n\n";

    // Analytic check: σ_F and σ_B should each equal σ_total/2
    double sigma_half     = r.sigma_total / 2.0;
    double asym_F         = (r.sigma_F - sigma_half) / r.error_F;
    double asym_B         = (r.sigma_B - sigma_half) / r.error_B;
    std::cout << "  σ_F pull from σ/2 = " << std::setprecision(3)
              << asym_F << " σ\n";
    std::cout << "  σ_B pull from σ/2 = " << asym_B    << " σ\n\n";

    // Pass/fail
    bool pass = significance < 2.0;
    std::cout << "  Result: " << (pass ? "PASS" : "FAIL")
              << " — A_FB is "
              << (pass ? "consistent with 0" : "SIGNIFICANTLY NONZERO")
              << " at 2σ level\n";

    if (!pass) {
        std::cout << "\n  *** Possible causes of nonzero A_FB:\n"
                  << "      - RNG not uniform over [-1,1]\n"
                  << "      - Integration bounds asymmetric\n"
                  << "      - Bug in matrix_element() sign of costheta term\n"
                  << "      - Bin edge at cosθ=0 mishandled\n";
    }
}

// ---------------------------------------------------------------------------
// Energy scan summary table
// ---------------------------------------------------------------------------

void print_asymmetry_scan_report(
        const std::vector<AsymmetryScanPoint>& scan) {

    std::cout << "\n=== A_FB Energy Scan ===\n\n";
    std::cout << std::left
              << std::setw(10) << "sqrt_s"
              << std::setw(8)  << "beta"
              << std::setw(14) << "sigma_F"
              << std::setw(14) << "sigma_B"
              << std::setw(12) << "A_FB"
              << std::setw(12) << "A_FB_err"
              << std::setw(10) << "signif."
              << std::setw(8)  << "status"
              << "\n";
    std::cout << std::string(88, '-') << "\n";

    int pass_count = 0;
    for (const auto& sp : scan) {
        const auto& r  = sp.hemi;
        double mu_mass = 0.1056;
        double s       = sp.sqrt_s * sp.sqrt_s;
        double beta    = std::sqrt(1.0 - 4.0 * mu_mass * mu_mass / s);
        double signif  = (r.A_FB_error > 0)
                       ? std::abs(r.A_FB) / r.A_FB_error : 0.0;
        save_asymmetry(sp.sqrt_s,sp.hemi.A_FB, sp.hemi.A_FB_error);
        
                       bool   pass    = signif < 2.0;
        if (pass) pass_count++;

        std::cout << std::left  << std::fixed
                  << std::setw(10) << std::setprecision(4) << sp.sqrt_s
                  << std::setw(8)  << std::setprecision(4) << beta
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << r.sigma_F
                  << std::setw(14) << r.sigma_B
                  << std::setw(12) << r.A_FB
                  << std::setw(12) << r.A_FB_error
                  << std::fixed << std::setprecision(2)
                  << std::setw(10) << signif
                  << std::setw(8)  << (pass ? "OK" : "FAIL")
                  << "\n";
    }

    // Overall summary
    std::cout << "\n  Pass rate: " << pass_count
              << " / "            << scan.size()
              << " points within 2σ of A_FB = 0\n";

    // Statistical expectation: with ~N/2 events per hemisphere, roughly
    // 95% of points should pass the 2σ test by chance alone even for a
    // perfect integrator. Flag if pass rate is too low.
    double pass_rate = static_cast<double>(pass_count) / scan.size();
    if (pass_rate >= 0.90)
        std::cout << "  RESULT: PASS — no evidence of sampling asymmetry\n";
    else
        std::cout << "  RESULT: FAIL — sampling asymmetry detected\n";

    // Extra: check that A_FB values are consistent with a normal distribution
    // centred on 0 by computing their mean and RMS across the scan
    double mean_AFB = 0.0, rms_AFB = 0.0;
    for (const auto& sp : scan) {
        mean_AFB += sp.hemi.A_FB / sp.hemi.A_FB_error;    // normalised by error
        rms_AFB  += std::pow(sp.hemi.A_FB / sp.hemi.A_FB_error, 2);
    }
    mean_AFB /= scan.size();
    rms_AFB   = std::sqrt(rms_AFB / scan.size());
    std::cout << "\n  Normalised A_FB distribution across scan:\n";
    std::cout << "    mean / σ_A  = " << std::fixed << std::setprecision(4)
              << mean_AFB << "  (expect ~0)\n";
    std::cout << "    RMS  / σ_A  = " << rms_AFB
              << "  (expect ~1)\n";
}