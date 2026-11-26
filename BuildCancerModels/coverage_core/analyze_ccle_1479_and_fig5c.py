#!/usr/bin/env python3
"""render_fig5c_3panels.py
Renders 3-panel Fig. 5c at the exact same canvas size figsize=(9.2, 5.1):
  Panel 1: Hallmark recovery (% of 135 reactions, median +/- IQR)
  Panel 2: Model size (Reactions per CSM |N|, median +/- IQR)
  Panel 3: Extraction runtime (Time per CSM in s + total cohort time, median +/- IQR)
With unpaired two-sided Mann-Whitney U square significance brackets (* / n.s.).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

RERUN_ROOT = '/scratch/pavan/q_4_4/rerun_69d2281'
OUT_DIR = os.path.join(RERUN_ROOT, 'results')
FIG_DIR = os.path.join(RERUN_ROOT, 'figures')


def add_bracket(ax, x1, x2, y, h, label, lw=1.4, fontsize=14.5):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color='#111f33', clip_on=False)
    va = 'center' if label == '*' else 'bottom'
    dy = h + (0.32 * h if label == '*' else 0.14 * h)
    fs = fontsize if label == '*' else fontsize - 3
    ax.text(0.5 * (x1 + x2), y + dy, label, ha='center', va=va,
            fontsize=fs, fontweight='bold', color='#111f33')


def main():
    df = pd.read_csv(os.path.join(OUT_DIR, 'ccle_1479_new_spectra_metrics.csv'))

    sp_c = df['SPECTRA_hall_cov'].to_numpy(float) * 100.0
    fc_c = df['Fastcore_W_hall_cov'].to_numpy(float) * 100.0
    sw_c = df['Swiftcore_hall_cov'].to_numpy(float) * 100.0

    sp_n = df['SPECTRA_N'].to_numpy(float)
    fc_n = df['Fastcore_W_N'].to_numpy(float)
    sw_n = df['Swiftcore_N'].to_numpy(float)

    sp_t = df['SPECTRA_time_4cpu_s'].to_numpy(float)
    fc_t = df['fcw_time_4cpu_s'].to_numpy(float)
    sw_t = df['Swiftcore_time_4cpu_s'].to_numpy(float)

    sig_cov_fc = '*' if mannwhitneyu(sp_c, fc_c, alternative='two-sided').pvalue < 0.05 else 'n.s.'
    sig_cov_sw = '*' if mannwhitneyu(sp_c, sw_c, alternative='two-sided').pvalue < 0.05 else 'n.s.'
    sig_n_fc   = '*' if mannwhitneyu(sp_n, fc_n, alternative='two-sided').pvalue < 0.05 else 'n.s.'
    sig_n_sw   = '*' if mannwhitneyu(sp_n, sw_n, alternative='two-sided').pvalue < 0.05 else 'n.s.'
    sig_rt_fc  = '*' if mannwhitneyu(sp_t, fc_t, alternative='two-sided').pvalue < 0.05 else 'n.s.'
    sig_rt_sw  = '*' if mannwhitneyu(sp_t, sw_t, alternative='two-sided').pvalue < 0.05 else 'n.s.'

    sns.set_theme(style='ticks', context='paper')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.2, 5.1))

    methods_lbl = ['SPECTRA', 'Fastcore', 'Swiftcore']
    x_pos = np.array([0, 1, 2])
    cols = ['#076d6d', '#80c080', '#08306b']
    edge_cols = ['#044747', '#2e6f40', '#051f45']
    inside_cols = ['white', '#112b11', 'white']
    bar_w = 0.65

    # --- Panel 1: Hallmark recovery ---
    cov_list = [sp_c, fc_c, sw_c]
    med_c = [np.median(v) for v in cov_list]
    q25_c = [np.percentile(v, 25) for v in cov_list]
    q75_c = [np.percentile(v, 75) for v in cov_list]

    ax1.bar(x_pos, med_c, width=bar_w, color=cols, edgecolor=edge_cols, linewidth=1.15, alpha=0.93)
    ax1.errorbar(x_pos, med_c,
                 yerr=[[m - q25 for m, q25 in zip(med_c, q25_c)],
                       [q75 - m for m, q75 in zip(med_c, q75_c)]],
                 fmt='none', ecolor='#111111', elinewidth=1.8, capsize=4.2, capthick=1.6, zorder=5)
    for i, (xp, m) in enumerate(zip(x_pos, med_c)):
        ax1.text(xp, 36.0, f"{m:.1f}%", ha='center', va='center',
                 fontsize=9.3, fontweight='bold', color=inside_cols[i])
    top_c = max(q75_c)
    add_bracket(ax1, 0, 1, top_c + 3.5, 2.5, sig_cov_fc, lw=1.4, fontsize=14.5)
    add_bracket(ax1, 0, 2, top_c + 12.0, 2.5, sig_cov_sw, lw=1.4, fontsize=14.5)
    ax1.set_title('Hallmark recovery', fontsize=12.6, fontweight='bold', pad=9)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(methods_lbl, fontsize=9.5, fontweight='bold')
    ax1.set_xlim(-0.56, 2.56)
    ax1.set_ylim(0, 105)
    ax1.set_yticks([0, 25, 50, 75, 100])
    ax1.set_ylabel('% of 135 reactions (median ± IQR)', fontsize=11.0, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=10.2)
    ax1.grid(axis='y', linestyle='--', color='grey', alpha=0.22)
    sns.despine(ax=ax1)

    # --- Panel 2: Extracted model size |N| (plotted in thousands for compact y-axis, exact integer inside bar) ---
    n_list = [sp_n / 1000.0, fc_n / 1000.0, sw_n / 1000.0]
    med_n = [np.median(v) for v in n_list]
    q25_n = [np.percentile(v, 25) for v in n_list]
    q75_n = [np.percentile(v, 75) for v in n_list]
    med_n_int = [int(round(np.median(v))) for v in [sp_n, fc_n, sw_n]]

    ax2.bar(x_pos, med_n, width=bar_w, color=cols, edgecolor=edge_cols, linewidth=1.15, alpha=0.93)
    ax2.errorbar(x_pos, med_n,
                 yerr=[[m - q25 for m, q25 in zip(med_n, q25_n)],
                       [q75 - m for m, q75 in zip(med_n, q75_n)]],
                 fmt='none', ecolor='#111111', elinewidth=1.8, capsize=4.2, capthick=1.6, zorder=5)
    for i, (xp, m_i) in enumerate(zip(x_pos, med_n_int)):
        ax2.text(xp, 3.25, f"{m_i:,}", ha='center', va='center',
                 fontsize=9.3, fontweight='bold', color=inside_cols[i])
    top_n = max(q75_n)
    add_bracket(ax2, 0, 1, top_n + 0.35, 0.24, sig_n_fc, lw=1.4, fontsize=14.5)
    add_bracket(ax2, 0, 2, top_n + 1.15, 0.24, sig_n_sw, lw=1.4, fontsize=14.5)
    ax2.set_title('Model size', fontsize=12.6, fontweight='bold', pad=9)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(methods_lbl, fontsize=9.5, fontweight='bold')
    ax2.set_xlim(-0.56, 2.56)
    ax2.set_ylim(0, 9.8)
    ax2.set_yticks([0, 2, 4, 6, 8])
    ax2.set_ylabel('Reactions, ×10³ (median ± IQR)', fontsize=11.0, fontweight='bold')
    ax2.tick_params(axis='y', labelsize=10.2)
    ax2.grid(axis='y', linestyle='--', color='grey', alpha=0.22)
    sns.despine(ax=ax2)

    # --- Panel 3: Extraction runtime ---
    rt_list = [sp_t, fc_t, sw_t]
    med_t = [np.median(v) for v in rt_list]
    q25_t = [np.percentile(v, 25) for v in rt_list]
    q75_t = [np.percentile(v, 75) for v in rt_list]

    ax3.bar(x_pos, med_t, width=bar_w, color=cols, edgecolor=edge_cols, linewidth=1.15, alpha=0.93)
    ax3.errorbar(x_pos, med_t,
                 yerr=[[m - q25 for m, q25 in zip(med_t, q25_t)],
                       [q75 - m for m, q75 in zip(med_t, q75_t)]],
                 fmt='none', ecolor='#111111', elinewidth=1.8, capsize=4.2, capthick=1.6, zorder=5)

    def fmt_tot(arr):
        tot_s = arr.sum()
        return f"{tot_s/60:.0f}m" if tot_s < 3600 else f"{tot_s/3600:.1f}h"

    top_t = max(q75_t)
    tot_labels = [f"{med_t[i]:.1f} s\n({fmt_tot(rt_list[i])})" for i in range(3)]
    for i, (xp, m, q25, q75, lbl) in enumerate(zip(x_pos, med_t, q25_t, q75_t, tot_labels)):
        if q25 > 0.36 * top_t:
            ax3.text(xp, 0.50 * q25, lbl, ha='center', va='center',
                     fontsize=9.1, fontweight='bold', color=inside_cols[i], linespacing=1.12)
        else:
            ax3.text(xp, q75 + 0.045 * top_t, lbl, ha='center', va='bottom',
                     fontsize=9.2, fontweight='bold', color='#1f1f1f', linespacing=1.10)

    h_b = 0.045 * top_t
    add_bracket(ax3, 0, 1, top_t + 0.08 * top_t, h_b, sig_rt_fc, lw=1.4, fontsize=14.5)
    add_bracket(ax3, 0, 2, top_t + 0.23 * top_t, h_b, sig_rt_sw, lw=1.4, fontsize=14.5)
    ax3.set_title('Extraction runtime', fontsize=12.6, fontweight='bold', pad=9)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(methods_lbl, fontsize=9.5, fontweight='bold')
    ax3.set_xlim(-0.56, 2.56)
    ax3.set_ylim(0, top_t * 1.42)
    ax3.set_ylabel('Time per CSM, s (median ± IQR)', fontsize=11.0, fontweight='bold')
    ax3.tick_params(axis='y', labelsize=10.2)
    ax3.grid(axis='y', linestyle='--', color='grey', alpha=0.22)
    sns.despine(ax=ax3)

    plt.tight_layout(w_pad=0.95)
    plt.savefig(os.path.join(FIG_DIR, 'cancer_5c_new_spectra_3bars.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'cancer_5c_new_spectra_3bars.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved 3-panel Fig. 5c cancer_5c_new_spectra_3bars.pdf/.png!", flush=True)


if __name__ == '__main__':
    main()
