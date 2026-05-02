import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import pandas as pd

# ── Load data ────────────────────────────────────────────────────────────────
# convergence.csv has no header (it's commented out in utils.cpp)
# columns written: logN, N, sigma_mc, sigma_analytic, mc_stat_error
conv = pd.read_csv(
    "results/convergence.csv",
    names=["logN", "N", "sigma_mc", "sigma_analytic", "mc_stat_error"]
)

# Derive missing columns from what utils.cpp wrote
conv["abs_diff"]     = (conv["sigma_mc"] - conv["sigma_analytic"]).abs()
conv["rel_diff_pct"] = conv["abs_diff"] / conv["sigma_analytic"] * 100.0

# performance.csv is written inline in main.cpp and has a proper header
df = pd.read_csv("results/performance.csv")

# serial_time is not in convergence.csv — pull it from performance.csv
# where threads==1, one row per logN
serial_times = (
    df[df["threads"] == 1]
    .set_index("logN")["time"]
)
conv["serial_time"] = conv["logN"].map(serial_times)

all_logN          = sorted(df["logN"].unique())
all_N             = sorted(df["N"].unique())
threads_available = sorted(df["threads"].unique())

cmap   = cm.get_cmap("plasma", len(all_logN))
colors = {logN: cmap(i) for i, logN in enumerate(all_logN)}

def get_N_data(logN):
    return df[df["logN"] == logN].sort_values("threads")

def amdahl_fit(sub):
    sp   = sub["speedup"].values
    th   = sub["threads"].values.astype(float)
    mask = th > 1
    if mask.sum() == 0:
        return np.nan
    fs = np.mean(
        (1.0/sp[mask] - 1.0/th[mask]) / (1.0 - 1.0/th[mask])
    )
    return fs

# ── Layout ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.suptitle(
    r"OpenMP Scaling Analysis: $e^+e^-\!\to\!\mu^+\mu^-$ Monte Carlo"
    "\nSpeedup, efficiency and physics consistency across N and thread count",
    fontsize=13, y=0.99
)
gs = gridspec.GridSpec(3, 3, hspace=0.48, wspace=0.33)
ax1 = fig.add_subplot(gs[0, 0])   # speedup vs threads, coloured by N
ax2 = fig.add_subplot(gs[0, 1])   # efficiency vs threads, coloured by N
ax3 = fig.add_subplot(gs[0, 2])   # wall time vs threads, coloured by N
ax4 = fig.add_subplot(gs[1, 0])   # speedup vs N, coloured by thread count
ax5 = fig.add_subplot(gs[1, 1])   # efficiency vs N, coloured by thread count
ax6 = fig.add_subplot(gs[1, 2])   # Amdahl serial fraction vs N
ax7 = fig.add_subplot(gs[2, 0])   # time vs N (log-log scaling)
ax8 = fig.add_subplot(gs[2, 1])   # sigma consistency across threads, per N
ax9 = fig.add_subplot(gs[2, 2])   # efficiency heatmap: N × threads

p_smooth = np.linspace(1, max(threads_available) * 1.4, 300)

# ── Panel 1: Speedup vs threads (one line per N) ─────────────────────────────
ax1.plot(p_smooth, p_smooth, '--', color='lightgray',
         linewidth=1.4, label='Ideal', zorder=1)
for logN in all_logN:
    sub = get_N_data(logN)
    fs  = amdahl_fit(sub)
    amd = 1.0 / (fs + (1.0 - fs) / p_smooth) if not np.isnan(fs) else None
    if amd is not None:
        ax1.plot(p_smooth, amd, '-', color=colors[logN],
                 linewidth=1.0, alpha=0.4)
    ax1.plot(sub["threads"], sub["speedup"], 'o-',
             color=colors[logN], linewidth=1.8, markersize=7,
             label=f'$N=10^{{{logN}}}$')
ax1.set_xlabel("Threads", fontsize=11)
ax1.set_ylabel(r"Speedup  $S = T_1/T_p$", fontsize=11)
ax1.set_title("Speedup vs threads", fontsize=11)
ax1.legend(fontsize=8, ncol=2)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(left=0); ax1.set_ylim(bottom=0)
ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# ── Panel 2: Efficiency vs threads ───────────────────────────────────────────
ax2.axhline(100, color='lightgray', linewidth=1.4,
            linestyle='--', label='Ideal')
ax2.axhline(80,  color='green', linewidth=1.1,
            linestyle=':', alpha=0.7, label='80% target')
for logN in all_logN:
    sub = get_N_data(logN)
    ax2.plot(sub["threads"], sub["efficiency"], 's-',
             color=colors[logN], linewidth=1.8, markersize=7,
             label=f'$N=10^{{{logN}}}$')
ax2.set_xlabel("Threads", fontsize=11)
ax2.set_ylabel(r"Efficiency  $E = S/p \times 100\%$", fontsize=11)
ax2.set_title("Parallel efficiency vs threads", fontsize=11)
ax2.legend(fontsize=8, ncol=2)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(left=0); ax2.set_ylim(0, 115)
ax2.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# ── Panel 3: Wall time vs threads ────────────────────────────────────────────
for logN in all_logN:
    sub      = get_N_data(logN)
    t_serial = sub[sub["threads"] == 1]["time"].values[0]
    ax3.plot(sub["threads"], sub["time"] * 1000, 'D-',
             color=colors[logN], linewidth=1.8, markersize=7,
             label=f'$N=10^{{{logN}}}$')
    ax3.plot(p_smooth, t_serial * 1000 / p_smooth, ':',
             color=colors[logN], linewidth=0.9, alpha=0.45)
ax3.set_xlabel("Threads", fontsize=11)
ax3.set_ylabel("Wall time [ms]", fontsize=11)
ax3.set_title("Wall-clock time vs threads", fontsize=11)
ax3.set_yscale('log')
ax3.legend(fontsize=8, ncol=2)
ax3.grid(True, which='both', alpha=0.3)
ax3.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# ── Panel 4: Speedup vs N (one line per thread count) ────────────────────────
thread_colors = {t: cm.get_cmap("cool")(i / max(1, len(threads_available)-1))
                 for i, t in enumerate(threads_available)}
for t in threads_available:
    sub = df[df["threads"] == t].sort_values("N")
    ax4.semilogx(sub["N"], sub["speedup"], 'o-',
                 color=thread_colors[t], linewidth=1.8, markersize=7,
                 label=f'{int(t)} thread{"s" if t>1 else ""}')
ax4.set_xlabel("N (events)", fontsize=11)
ax4.set_ylabel(r"Speedup  $S$", fontsize=11)
ax4.set_title("Speedup vs N per thread count", fontsize=11)
ax4.legend(fontsize=8)
ax4.grid(True, which='both', alpha=0.3)
ax4.set_ylim(bottom=0)

# ── Panel 5: Efficiency vs N ──────────────────────────────────────────────────
ax5.axhline(80, color='green', linewidth=1.1,
            linestyle=':', alpha=0.7, label='80% target')
for t in threads_available:
    if t == 1:
        continue
    sub = df[df["threads"] == t].sort_values("N")
    ax5.semilogx(sub["N"], sub["efficiency"], 's-',
                 color=thread_colors[t], linewidth=1.8, markersize=7,
                 label=f'{int(t)} threads')
ax5.set_xlabel("N (events)", fontsize=11)
ax5.set_ylabel(r"Efficiency [%]", fontsize=11)
ax5.set_title("Efficiency vs N per thread count", fontsize=11)
ax5.legend(fontsize=8)
ax5.grid(True, which='both', alpha=0.3)
ax5.set_ylim(0, 115)

# ── Panel 6: Amdahl serial fraction vs N ─────────────────────────────────────
fs_vals, N_vals = [], []
for logN in all_logN:
    sub = get_N_data(logN)
    fs  = amdahl_fit(sub)
    if not np.isnan(fs):
        fs_vals.append(fs * 100)
        N_vals.append(10**logN)

ax6.semilogx(N_vals, fs_vals, 'D-', color='crimson',
             linewidth=2.0, markersize=9, label='Estimated $f_s$')
ax6.set_xlabel("N (events)", fontsize=11)
ax6.set_ylabel("Serial fraction $f_s$ [%]", fontsize=11)
ax6.set_title("Serial fraction vs N\n(Amdahl fit)", fontsize=11)
ax6.grid(True, which='both', alpha=0.3)
ax6.set_ylim(bottom=0)
# Annotate S_max for each point
for N_v, fs_v in zip(N_vals, fs_vals):
    s_max = 100.0 / fs_v if fs_v > 0 else float('inf')
    ax6.annotate(f"$S_{{max}}$={s_max:.1f}×",
                 xy=(N_v, fs_v), xytext=(4, 6),
                 textcoords='offset points', fontsize=7.5,
                 color='crimson')

# ── Panel 7: Time vs N (log-log, one line per thread count) ──────────────────
for t in threads_available:
    sub = df[df["threads"] == t].sort_values("N")
    ax7.loglog(sub["N"], sub["time"] * 1000, 'o-',
               color=thread_colors[t], linewidth=1.8, markersize=7,
               label=f'{int(t)} thread{"s" if t>1 else ""}')
# Reference O(N) line
N_ref = np.array([1e3, 1e7])
t_ref = df[(df["threads"]==1) & (df["N"]==1000)]["time"].values
if len(t_ref):
    ax7.loglog(N_ref, t_ref[0] * 1000 * N_ref / 1e3, 'k--',
               linewidth=1.2, alpha=0.5, label=r'$\mathcal{O}(N)$')
ax7.set_xlabel("N (events)", fontsize=11)
ax7.set_ylabel("Wall time [ms]", fontsize=11)
ax7.set_title("Time vs N (log-log)\nshould be O(N)", fontsize=11)
ax7.legend(fontsize=8, ncol=2)
ax7.grid(True, which='both', alpha=0.3)

# ── Panel 8: sigma consistency (deviation from serial, in units of error) ────
for logN in all_logN:
    sub      = get_N_data(logN)
    sig_ser  = sub[sub["threads"] == 1]["sigma"].values[0]
    err_ser  = sub[sub["threads"] == 1]["error"].values[0]
    dev      = (sub["sigma"] - sig_ser) / err_ser
    ax8.plot(sub["threads"], dev, 'o-',
             color=colors[logN], linewidth=1.8, markersize=7,
             label=f'$N=10^{{{logN}}}$')
ax8.axhline( 0, color='black',  linewidth=1.0)
ax8.axhline( 2, color='orange', linewidth=1.2,
             linestyle='--', label=r'$\pm 2\sigma$')
ax8.axhline(-2, color='orange', linewidth=1.2, linestyle='--')
ax8.set_xlabel("Threads", fontsize=11)
ax8.set_ylabel(r"$(\sigma_{OMP} - \sigma_{serial})/\sigma_{err}$", fontsize=11)
ax8.set_title("Physics consistency across threads", fontsize=11)
ax8.legend(fontsize=8, ncol=2)
ax8.grid(True, alpha=0.3)
ax8.set_ylim(-4, 4)
ax8.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

# ── Panel 9: Efficiency heatmap (N × threads) ─────────────────────────────────
pivot = df[df["threads"] > 1].pivot_table(
    index="logN", columns="threads", values="efficiency"
)
im = ax9.imshow(pivot.values, aspect='auto', cmap='RdYlGn',
                vmin=50, vmax=100, origin='lower')
ax9.set_xticks(range(len(pivot.columns)))
ax9.set_xticklabels([f'{int(t)}' for t in pivot.columns], fontsize=9)
ax9.set_yticks(range(len(pivot.index)))
ax9.set_yticklabels([f'$10^{{{int(n)}}}$' for n in pivot.index], fontsize=9)
ax9.set_xlabel("Threads", fontsize=11)
ax9.set_ylabel("N (events)", fontsize=11)
ax9.set_title("Efficiency heatmap [%]\n(green = efficient)", fontsize=11)
plt.colorbar(im, ax=ax9, label="Efficiency [%]")
# Annotate each cell
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.values[i, j]
        if not np.isnan(val):
            ax9.text(j, i, f"{val:.1f}",
                     ha='center', va='center', fontsize=8.5,
                     color='black' if val > 65 else 'white')

plt.savefig("plots/performance.pdf",bbox_inches='tight')
plt.show()

# ── Console summary ───────────────────────────────────────────────────────────
print("\n=== Scaling Summary ===")
print(f"{'N':<10} {'Threads':<10} {'Time(ms)':<12} "
      f"{'Speedup':<10} {'Efficiency':<12} {'Physics'}")
print("-" * 68)
for _, row in df.sort_values(["N","threads"]).iterrows():
    sub_ser = df[(df["N"]==row["N"]) & (df["threads"]==1)]
    if sub_ser.empty:
        continue
    sig_ser = sub_ser["sigma"].values[0]
    err_ser = sub_ser["error"].values[0]
    dev     = abs(row["sigma"] - sig_ser) / err_ser if err_ser > 0 else 0
    status  = "OK" if dev < 2.0 else "CHECK"
    print(f"10^{int(np.log10(row['N'])):<7} "
          f"{int(row['threads']):<10} "
          f"{row['time']*1000:<12.2f} "
          f"{row['speedup']:<10.3f} "
          f"{row['efficiency']:<12.2f} "
          f"{status}")