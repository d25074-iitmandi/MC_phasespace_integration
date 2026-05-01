#include "utils.h"
#include <fstream>

void save_results(int threads, double time, double speedup) {
    std::ofstream file("results/scaling.txt", std::ios::app);
    file << threads << " " << time << " " << speedup << "\n";
}
