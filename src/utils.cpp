#include "utils.hpp"
#include <fstream>

void save_convergence_results(int logN, int N, double mc_sigma, double a_sigma, double mc_error){
    std::ofstream conv_file("results/convergence.csv",std::ios::app);
    //conv_file << "logN,N,sigma_mc,sigma_analytic,mc_stat_error,abs_diff,rel_diff_pct\n";
    conv_file << logN << ","<< N <<","<< mc_sigma << "," << a_sigma << ","<< mc_error << "\n";

}
void save_asymmetry(double sqrt_s, double a_fb, double a_fb_err){
    std::ofstream afb_file("results/afb.csv", std::ios::app);
    //afb_file << "sqrt_s,beta,sigma_F,sigma_B,sigma_total, A_FB,A_FB_error,significance\n";
    afb_file << sqrt_s << "," << a_fb << "," << a_fb_err << "\n";
}
