#pragma once
#include <vector>

void save_results(int threads, double time, double speedup);
void save_convergence_results(int logN, int N, double mc_sigma, double a_sigma, double mc_error);
void save_asymmetry(double sqrt_s, double a_fb, double a_fb_err);
