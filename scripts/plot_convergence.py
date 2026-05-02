import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy.stats import linregress

# ── Load data ──────────────────────────────────────────────────────────────
df = pd.read_csv("results/convergence.csv", header=None)

N           = df.iloc[:,1]
sigma_mc    = df.iloc[:,2]
sigma_an    = df.iloc[:,3]
stat_err    = df.iloc[:,4]


analytic_val = sigma_an[0]   # constant across rows

# ── Fit 1/sqrt(N) scaling to the MC stat error ─────────────────────────────
log_N   = np.log10(N)
log_err = np.log10(stat_err)
slope, intercept, r2, *_ = linregress(log_N, log_err)
N_fit   = np.logspace(2.8, 7.2, 200)
err_fit = 10**intercept * N_fit**slope

# ── Layout ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 10))
fig.suptitle(
    r"Monte Carlo Convergence: $e^+e^- \to \mu^+\mu^-$"
    f"\n" + r"$\sqrt{s}$ = " + f"{df.iloc[:,2].iloc[0]:.3e} reference  |  "
    r"Analytic $\sigma$ = " + f"{analytic_val:.4e} GeV$^{{-2}}$",
    fontsize=13, y=0.98
)
gs = gridspec.GridSpec(1, 2, hspace=0.38, wspace=0.32)
ax1 = fig.add_subplot(gs[0, 0])   # sigma_MC vs N
ax2 = fig.add_subplot(gs[0, 1])   # stat error vs N (log-log)


# ── Panel 1: sigma_MC converging to analytic ───────────────────────────────
ax1.semilogx(N, sigma_mc * 1e9, 'o-', color='steelblue',
             linewidth=1.8, markersize=7, label=r'MC $\sigma$')
ax1.fill_between(N,
                 (sigma_mc - stat_err) * 1e9,
                 (sigma_mc + stat_err) * 1e9,
                 alpha=0.25, color='steelblue', label=r'MC $\pm 1\sigma$')
ax1.axhline(analytic_val * 1e9, color='crimson', linewidth=1.8,
            linestyle='--', label='Analytic')
ax1.set_xlabel("Number of events  $N$", fontsize=11)
ax1.set_ylabel(r"$\sigma$ [$\times 10^{-9}$ GeV$^{-2}$]", fontsize=11)
ax1.set_title("Cross-section convergence", fontsize=11)
ax1.legend(fontsize=9)
ax1.grid(True, which='both', alpha=0.3)

# ── Panel 2: stat error vs N (log-log, expect slope = -0.5) ────────────────
ax2.loglog(N, stat_err, 'o', color='darkorange',
           markersize=8, zorder=5, label='MC stat error')
ax2.loglog(N_fit, err_fit, '--', color='gray',
           linewidth=1.6,
           label=f'Fit: slope = {slope:.3f}\n(ideal = −0.500)')
ax2.set_xlabel("Number of events  $N$", fontsize=11)
ax2.set_ylabel(r"Statistical error [GeV$^{-2}$]", fontsize=11)
ax2.set_title(r"Error scaling: $1/\sqrt{N}$ test", fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, which='both', alpha=0.3)

# Annotate slope
ax2.annotate(f"slope = {slope:.4f}\n" + r"$R^2$ = " + f"{r2**2:.6f}",
             xy=(0.62, 0.72), xycoords='axes fraction',
             fontsize=9,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))



plt.savefig("plots/convergence.pdf")#, bbox_inches='tight')
plt.show()
