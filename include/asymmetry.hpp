#pragma once
#include <vector>

struct HemisphereResult {
    double sigma_F;          // forward cross-section  (cosθ > 0)
    double sigma_B;          // backward cross-section (cosθ < 0)
    double error_F;          // MC statistical error on sigma_F
    double error_B;          // MC statistical error on sigma_B
    double A_FB;             // (sigma_F - sigma_B) / (sigma_F + sigma_B)
    double A_FB_error;       // propagated statistical error on A_FB
    double sigma_total;      // sigma_F + sigma_B (consistency check vs direct MC)
};

struct AsymmetryScanPoint {
    double sqrt_s;
    HemisphereResult hemi;
};

// Single energy asymmetry check
HemisphereResult forward_backward_asymmetry(int N, double sqrt_s);

// Asymmetry across an energy range — A_FB should be 0 everywhere
std::vector<AsymmetryScanPoint> asymmetry_energy_scan(
    int N_per_point,
    const std::vector<double>& energies);

// Detailed report
void print_asymmetry_report(
    const HemisphereResult& result,
    double sqrt_s);

// Summary table for energy scan
void print_asymmetry_scan_report(
    const std::vector<AsymmetryScanPoint>& scan);