#include "utils.hpp"
#include <fstream>
//Writing the output for scaling
void save_results(int threads, double time, double speedup) {
    std::ofstream file("results/scaling.txt", std::ios::app);
    file << threads << " " << time << " " << speedup << "\n";
}
void save_convergence_results(int logN, int N, double mc_sigma, double a_sigma, double mc_error, double abs_err, double rel_diff){
    std::ofstream conv_file("results/convergence.csv",std::ios::app);
    //conv_file << "logN,N,sigma_mc,sigma_analytic,mc_stat_error,abs_diff,rel_diff_pct\n";
    conv_file << logN << ","<< N <<","<< mc_sigma << "," << a_sigma << ","<< mc_error << ","<< abs_err << "," << abs_err/a_sigma*100 << "\n";

}
void save_asymmetry(double sqrt_s, double beta, double sigma_f, double sigma_b, double sigma_total, double a_fb, double a_fb_err, double sig){
    std::ofstream afb_file("results/afb.csv", std::ios::app);
    //afb_file << "sqrt_s,beta,sigma_F,sigma_B,sigma_total, A_FB,A_FB_error,significance\n";
    afb_file << sqrt_s << "," << beta << ","<< sigma_f << "," << sigma_b << "," << sigma_total << "," << a_fb << "," << a_fb_err << "," << sig << "\n";
}
