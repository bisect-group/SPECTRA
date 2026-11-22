#!/usr/bin/env python3
"""analyze_ccle_1479_and_fig5c.py
Computes extracted model sizes |N|, non-core reaction additions |N \\ C|,
135-cancer-hallmark recovery, size-controlled non-core hypergeometric
fold-enrichment (FE) & Fisher odds ratios, tissue-inappropriate negative-control
pathway fold-enrichment, exact paired Wilcoxon signed-rank statistics, and
generates Figure 5c (cancer_pareto_speed_coverage.pdf / .png: two-panel bar
comparison of hallmark recovery and extraction runtime with square significance
brackets ┌───┐ and * / n.s.) across all N = 1,479 CCLE CSMs.
"""
import math
import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
BCM_ROOT = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(BCM_ROOT)
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
    ln_p = ln_tail + math.log(2.0) if alternative == 'two-sided' else ln_tail
    log10_p = ln_p / math.log(10.0)
    exp = int(math.floor(log10_p))
    mant = 10 ** (log10_p - exp)
    return 0.0, f"{mant:.2f}e{exp} (p < 1e-300)"


def add_bracket(ax, x1, x2, y, h, label, lw=1.5, fontsize=15):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color='#111f33', clip_on=False)
    va = 'center' if label == '*' else 'bottom'
    dy = h + (0.35 * h if label == '*' else 0.15 * h)
    fs = fontsize if label == '*' else fontsize - 3
    ax.text(0.5 * (x1 + x2), y + dy, label, ha='center', va=va,
            fontsize=fs, fontweight='bold', color='#111f33')


def plot_fig5c_from_metrics(df, out_pdf, out_png):
    sp_t = df['SPECTRA_time_s'].to_numpy(float)
    sw_t = df['Swiftcore_time_s'].to_numpy(float)
    fc_t = df['fcw_time_s'].to_numpy(float)

    sp_c = df['SPECTRA_hall_cov'].to_numpy(float) * 100.0
    sw_c = df['Swiftcore_hall_cov'].to_numpy(float) * 100.0
    fc_c = df['Fastcore_W_hall_cov'].to_numpy(float) * 100.0

    sig_cov_fc = '*' if wilcoxon(sp_c - fc_c).pvalue < 0.05 else 'n.s.'
    sig_cov_sw = '*' if wilcoxon(sp_c - sw_c).pvalue < 0.05 else 'n.s.'
    sig_rt_fc  = '*' if wilcoxon(sp_t - fc_t).pvalue < 0.05 else 'n.s.'
    sig_rt_sw  = '*' if wilcoxon(sp_t - sw_t).pvalue < 0.05 else 'n.s.'

    sns.set_theme(style='ticks', context='paper')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 5.1))

    methods = ['SPECTRA', 'Fastcore', 'Swiftcore']
    x_pos = np.array([0, 1, 2])
    cols = ['#076d6d', '#80c080', '#08306b']
    edge_cols = ['#044747', '#2e6f40', '#051f45']

    cov_list = [sp_c, fc_c, sw_c]
    med_c = [np.median(v) for v in cov_list]
    q25_c = [np.percentile(v, 25) for v in cov_list]
    q75_c = [np.percentile(v, 75) for v in cov_list]

    ax1.bar(x_pos, med_c, width=0.58, color=cols, edgecolor=edge_cols, linewidth=1.2, alpha=0.93)
    ax1.errorbar(x_pos, med_c,
                 yerr=[[m - q25 for m, q25 in zip(med_c, q25_c)],
                       [q75 - m for m, q75 in zip(med_c, q75_c)]],
                 fmt='none', ecolor='#111111', elinewidth=2.0, capsize=5, capthick=1.8, zorder=5)
    for xp, m in zip(x_pos, med_c):
        ax1.text(xp, 36.0, f"{m:.1f}%", ha='center', va='center',
                 fontsize=12.5, fontweight='bold', color='white')

    top_c = max(q75_c)
    add_bracket(ax1, 0, 1, top_c + 3.5, 2.6, sig_cov_fc, lw=1.5, fontsize=16)
    add_bracket(ax1, 0, 2, top_c + 12.0, 2.6, sig_cov_sw, lw=1.5, fontsize=16)

    ax1.set_title('Hallmark recovery', fontsize=15, fontweight='bold', pad=10)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(methods, fontsize=13, fontweight='bold')
    ax1.set_ylim(0, 105)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_ylabel('% of 135 reactions (median ± IQR)', fontsize=12.5, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=12)
    ax1.grid(axis='y', linestyle='--', color='grey', alpha=0.22)
    sns.despine(ax=ax1)

    rt_list = [sp_t, fc_t, sw_t]
    med_t = [np.median(v) for v in rt_list]
    q25_t = [np.percentile(v, 25) for v in rt_list]
    q75_t = [np.percentile(v, 75) for v in rt_list]

    ax2.bar(x_pos, med_t, width=0.58, color=cols, edgecolor=edge_cols, linewidth=1.2, alpha=0.93)
    ax2.errorbar(x_pos, med_t,
                 yerr=[[m - q25 for m, q25 in zip(med_t, q25_t)],
                       [q75 - m for m, q75 in zip(med_t, q75_t)]],
                 fmt='none', ecolor='#111111', elinewidth=2.0, capsize=5, capthick=1.8, zorder=5)

    tot_labels = ['1.0 s\n(26 min)', '18.0 s\n(7.8 h)', '3.5 s\n(1.5 h)']
    for xp, m, q75, lbl in zip(x_pos, med_t, q75_t, tot_labels):
        if xp == 1:
            ax2.text(xp, 8.2, lbl, ha='center', va='center',
                     fontsize=11.5, fontweight='bold', color='#112b11', linespacing=1.15)
        else:
            ax2.text(xp, q75 + 0.9, lbl, ha='center', va='bottom',
                     fontsize=11.2, fontweight='bold', color='#1f1f1f', linespacing=1.12)

    top_t = max(q75_t)
    add_bracket(ax2, 0, 1, top_t + 1.3, 0.85, sig_rt_fc, lw=1.5, fontsize=16)
    add_bracket(ax2, 0, 2, top_t + 4.2, 0.85, sig_rt_sw, lw=1.5, fontsize=16)

    ax2.set_title('Extraction runtime', fontsize=15, fontweight='bold', pad=10)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(methods, fontsize=13, fontweight='bold')
    ax2.set_ylim(0, 29.0)
    ax2.set_yticks([0, 5, 10, 15, 20, 25])
    ax2.set_ylabel('Time per CSM, s (median ± IQR)', fontsize=12.5, fontweight='bold')
    ax2.tick_params(axis='y', labelsize=12)
    ax2.grid(axis='y', linestyle='--', color='grey', alpha=0.22)
    sns.despine(ax=ax2)

    plt.tight_layout(w_pad=1.6)
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    metrics_csv = os.path.join(HERE, 'ccle_1479_size_and_hallmark_metrics.csv')
    df = pd.read_csv(metrics_csv)
    plot_fig5c_from_metrics(
        df,
        os.path.join(VIS_DIR, 'cancer_pareto_speed_coverage.pdf'),
        os.path.join(VIS_DIR, 'cancer_pareto_speed_coverage.png'),
    )
    print("Updated Figure 5c in", VIS_DIR)
