#include "utils.hpp"
#include <fstream>
//Writing the output for scaling
void save_results(int threads, double time, double speedup) {
    std::ofstream file("results/scaling.txt", std::ios::app);
    file << threads << " " << time << " " << speedup << "\n";
}
void save_convergence_results(int logN, double mc_sigma, double a_sigma, double abs_err, double mc_error){
    std::ofstream file("results/convergence.txt",std::ios::app);
    file << logN << " "<< mc_sigma << " " << a_sigma << " "<< abs_err << " " << mc_error << "\n";

}
