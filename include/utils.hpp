#pragma once
#include <vector>

void save_results(int threads, double time, double speedup);
void save_convergence_results(int logN, double mc_sigma, double a_sigma, double abs_err, double mc_error);
