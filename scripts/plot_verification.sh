#!/bin/bash
# run_plots.sh
mkdir -p plots
echo "Plotting convergence..."
python3 scripts/plot_convergence.py
echo "Plotting asymmetry..."
python3 scripts/plot_asymmetry.py
echo "Plotting performance analysis..."
python3 scripts/plot_performance.py
echo "All plots saved to plots/"