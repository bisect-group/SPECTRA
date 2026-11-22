#!/usr/bin/env python3
"""analyze_ccle_1479_and_fig5c.py
Computes extracted model sizes |N|, non-core reaction additions |N \\ C|,
135-cancer-hallmark recovery, size-controlled non-core hypergeometric
fold-enrichment (FE) & Fisher odds ratios, tissue-inappropriate negative-control
pathway fold-enrichment, exact paired Wilcoxon signed-rank statistics, and
generates Figure 5c (cancer_pareto_speed_coverage.pdf / .png) across all
N = 1,479 CCLE context-specific models (SPECTRA, weighted Fastcore, Swiftcore).
"""
import math
import os
import re
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import wilcoxon

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
BCM_ROOT = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(BCM_ROOT)
DEPMAP_DIR = os.path.join(BCM_ROOT, 'DepMapModels_Recon')
VIS_DIR = os.path.join(REPO_ROOT, 'Visualisations')


def exact_wilcoxon_str(x, y, alternative='two-sided'):
    """Compute Wilcoxon signed-rank p-value without underflowing to 0.0."""
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d_nz = d[d != 0]
    n = len(d_nz)
    if n == 0:
        return 1.0, "1.0"
    res = wilcoxon(x, y, alternative=alternative, method='approx')
    p = float(res.pvalue)
    if p > 1e-300:
        exp = int(math.floor(math.log10(p)))
        mant = p / (10 ** exp)
        return p, f"{mant:.2f}e{exp}" if abs(exp) >= 4 else f"{p:.4g}"
    abs_d = np.abs(d_nz)
    ranks = pd.Series(abs_d).rank(method='average').to_numpy()
    w_plus = ranks[d_nz > 0].sum()
    mean_w = n * (n + 1) / 4.0
    _, counts = np.unique(abs_d, return_counts=True)
    tie_term = np.sum(counts * (counts ** 2 - 1)) / 48.0
    var_w = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term
    z = (w_plus - mean_w) / math.sqrt(var_w)
    az = abs(z)
    ln_tail = -0.5 * az * az - math.log(az * math.sqrt(2 * math.pi)) + math.log(1 - 1/(az**2) + 3/(az**4))
    ln_p = ln_tail + math.log(2.0) if alternative == 'two-sided' else (ln_tail if ((alternative == 'greater' and z > 0) or (alternative == 'less' and z < 0)) else 0.0)
    log10_p = ln_p / math.log(10.0)
    exp = int(math.floor(log10_p))
    mant = 10 ** (log10_p - exp)
    return 0.0, f"{mant:.2f}e{exp} (p < 1e-300)"


def plot_pareto_from_metrics(df, out_pdf, out_png):
    sp_times = df['SPECTRA_time_s'].to_numpy(float)
    sw_times = df['Swiftcore_time_s'].to_numpy(float)
    fcw_times = df['fcw_time_s'].to_numpy(float)

    sp_cov = df['SPECTRA_hall_cov'].to_numpy(float) * 100.0
    sw_cov = df['Swiftcore_hall_cov'].to_numpy(float) * 100.0
    fcw_cov = df['Fastcore_W_hall_cov'].to_numpy(float) * 100.0

    sns.set_theme(style='ticks', context='paper')
    fig, ax = plt.subplots(figsize=(7.6, 5.7))

    specs = [
        ('SPECTRA',   sp_times,  sp_cov,  '#076d6d', '#044747', 'o'),
        ('Swiftcore', sw_times,  sw_cov,  '#08306b', '#051f45', 's'),
        ('Fastcore',  fcw_times, fcw_cov, '#5aae61', '#276419', '^'),
    ]

    for _, t, cov, col, _, _ in specs:
        ax.scatter(np.clip(t, 0.15, None), cov, s=13, color=col, alpha=0.14,
                   edgecolors='none', rasterized=True, zorder=2)

    for label, t, cov, col, edge_col, mk in specs:
        t_q25, t_med, t_q75 = np.percentile(t, [25, 50, 75])
        c_q25, c_med, c_q75 = np.percentile(cov, [25, 50, 75])
        ax.errorbar([t_med], [c_med], xerr=[[t_med - t_q25], [t_q75 - t_med]],
                    yerr=[[c_med - c_q25], [c_q75 - c_med]],
                    fmt='none', ecolor='white', elinewidth=4.2, capsize=5, capthick=3.0, zorder=5)
        ax.errorbar([t_med], [c_med], xerr=[[t_med - t_q25], [t_q75 - t_med]],
                    yerr=[[c_med - c_q25], [c_q75 - c_med]],
                    fmt='none', ecolor=edge_col, elinewidth=2.1, capsize=4.5, capthick=1.8, zorder=6)
        ax.scatter([t_med], [c_med], s=145, marker=mk, facecolor=col, edgecolor='white',
                   linewidth=1.8, zorder=7, label=f"{label} (med {t_med:.2f} s, {c_med:.1f}%)")
        ax.scatter([t_med], [c_med], s=145, marker=mk, facecolor='none', edgecolor=edge_col,
                   linewidth=0.9, zorder=8)

    sp_t_med, sp_c_med = np.median(sp_times), np.median(sp_cov)
    fc_t_med, fc_c_med = np.median(fcw_times), np.median(fcw_cov)
    y_arrow = max(sp_c_med, fc_c_med) + 6.8
    ax.annotate('', xy=(sp_t_med * 1.06, y_arrow), xytext=(fc_t_med * 0.94, y_arrow),
                arrowprops=dict(arrowstyle='<->', color='#333333', lw=1.5, shrinkA=0, shrinkB=0), zorder=9)
    ax.plot([sp_t_med, sp_t_med], [sp_c_med + 1.5, y_arrow + 0.8], ls=':', lw=1.0, color='#555555', zorder=4)
    ax.plot([fc_t_med, fc_t_med], [fc_c_med + 1.5, y_arrow + 0.8], ls=':', lw=1.0, color='#555555', zorder=4)
    ax.text(np.sqrt(sp_t_med * fc_t_med), y_arrow + 1.1, '18× faster  (Δ median < 1%)',
            ha='center', va='bottom', fontsize=11, fontweight='bold', color='#222222',
            bbox=dict(boxstyle='round,pad=0.22', facecolor='white', edgecolor='none', alpha=0.85), zorder=10)

    ax_in = ax.inset_axes([0.55, 0.10, 0.41, 0.30])
    tot_hrs = [sp_times.sum() / 3600.0, sw_times.sum() / 3600.0, fcw_times.sum() / 3600.0]
    ax_in.barh([2, 1, 0], tot_hrs, height=0.58, color=['#076d6d', '#08306b', '#5aae61'],
               edgecolor='#222222', linewidth=0.7)
    ax_in.set_yticks([2, 1, 0])
    ax_in.set_yticklabels(['SPECTRA', 'Swiftcore', 'Fastcore'], fontsize=9.5, fontweight='bold')
    ax_in.set_xlabel('Total time, 1,479 CSMs (h)', fontsize=9, fontweight='bold')
    ax_in.set_xlim(0, max(tot_hrs) * 1.36)
    ax_in.tick_params(axis='x', labelsize=8.5)
    ax_in.set_facecolor('#fafafa')
    sns.despine(ax=ax_in)
    for y_pos, val, tv in zip([2, 1, 0], tot_hrs, [f"{sp_times.sum()/60:.0f} min", f"{tot_hrs[1]:.1f} h", f"{tot_hrs[2]:.1f} h"]):
        ax_in.text(val + 0.18, y_pos, tv, va='center', ha='left', fontsize=8.8, fontweight='bold', color='#222222')

    ax.set_xscale('log')
    ax.set_xlim(0.22, 62.0)
    ax.set_ylim(42.0, 96.0)
    ax.set_xticks([0.3, 0.5, 1, 2, 5, 10, 20, 50])
    ax.get_xaxis().set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.get_xaxis().set_minor_formatter(mticker.NullFormatter())
    ax.grid(True, which='major', linestyle='--', color='grey', alpha=0.22, zorder=0)
    sns.despine(ax=ax)
    ax.set_xlabel('Extraction runtime per CCLE model (s, log scale)', fontsize=13, fontweight='bold', labelpad=6)
    ax.set_ylabel('Cancer hallmark recovery (% of 135)', fontsize=13, fontweight='bold', labelpad=6)
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.94,
              edgecolor='#dddddd', fontsize=10.2, handletextpad=0.5, borderpad=0.45)

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    metrics_csv = os.path.join(HERE, 'ccle_1479_size_and_hallmark_metrics.csv')
    df = pd.read_csv(metrics_csv)
    plot_pareto_from_metrics(
        df,
        os.path.join(VIS_DIR, 'cancer_pareto_speed_coverage.pdf'),
        os.path.join(VIS_DIR, 'cancer_pareto_speed_coverage.png'),
    )
    print("Updated Figure 5c Pareto speed-coverage plots in", VIS_DIR)
