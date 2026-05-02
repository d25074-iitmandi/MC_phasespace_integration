import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy.stats import norm

# ── Load data ──────────────────────────────────────────────────────────────
df = pd.read_csv("results/afb.csv", header=None)

sqrt_s    = df.iloc[:,0]
A_FB      = df.iloc[:,1]
A_err     = df.iloc[:,2]

# Normalised pulls: A_FB / sigma_A  (should be N(0,1))
pulls_norm = A_FB / A_err

mu_mass   = 0.1056   # GeV
threshold = 2 * mu_mass

# ── Layout ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 11))
fig.suptitle(
    r"Forward-Backward Asymmetry: $e^+e^- \to \mu^+\mu^-$"
    "\n" + r"QED tree level predicts $A_{FB} = 0$ exactly",
    fontsize=13, y=0.98
)
gs = gridspec.GridSpec(1, 1, hspace=0.42, wspace=0.35)
ax1 = fig.add_subplot(gs[0, :2])   # A_FB vs sqrt_s  (main)

# ── Panel 1: A_FB vs sqrt_s ────────────────────────────────────────────────
ax1.axhline(0, color='crimson', linewidth=2.0, linestyle='--',
            label='QED prediction: $A_{FB} = 0$', zorder=3)
ax1.fill_between(sqrt_s, -2 * A_err, 2 * A_err,
                 alpha=0.20, color='steelblue',
                 label=r'MC $\pm 2\sigma$ band')
ax1.fill_between(sqrt_s, -A_err, A_err,
                 alpha=0.35, color='steelblue',
                 label=r'MC $\pm 1\sigma$ band')
ax1.errorbar(sqrt_s, A_FB, yerr=A_err,
             fmt='o', color='steelblue', markersize=6,
             capsize=4, linewidth=1.4, zorder=5, label='$A_{FB}$ (MC)')
ax1.axvline(threshold, color='gray', linewidth=1.2,
            linestyle=':', alpha=0.7,
            label=f'Threshold $2m_\\mu$ = {threshold:.4f} GeV')
ax1.set_xlabel(r"$\sqrt{s}$  [GeV]", fontsize=12)
ax1.set_ylabel(r"$A_{FB}$", fontsize=12)
ax1.set_title(r"Forward-Backward Asymmetry vs $\sqrt{s}$", fontsize=11)
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xscale('log')

# Annotate pass rate
n_pass = np.sum(np.abs(A_FB) < 2 * A_err)
ax1.annotate(
    f"Pass rate (2σ): {n_pass}/{len(A_FB)}",
    xy=(0.03, 0.90), xycoords='axes fraction', fontsize=10,
    bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
              edgecolor='gray', alpha=0.9)
)


plt.savefig("plots/asymmetry.pdf", bbox_inches='tight')
plt.show()
