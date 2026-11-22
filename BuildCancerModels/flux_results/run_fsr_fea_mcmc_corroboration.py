#!/usr/bin/env python3
"""run_fsr_fea_mcmc_corroboration.py
Computes manuscript-faithful Flux Span Ratio (FSR) dysregulation across the 10
paired cancer vs. non-cancerous CCLE models, runs upper-tail one-sided
hypergeometric Flux Enrichment Analysis (FEA, Benjamini-Hochberg FDR < 0.05),
quantifies reaction- and pathway-level corroboration against 80-chain RHMC MCMC
uniform flux sampling, exports Data S4-S6 multi-sheet Excel workbooks, and
renders Figure 5d (bubble_fea_final.pdf) with outer concentric MCMC rings.
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

HERE = os.path.dirname(os.path.abspath(__file__))
BCM_ROOT = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(BCM_ROOT)
VIS_DIR = os.path.join(REPO_ROOT, 'Visualisations')

RING_COLORS = {
    'orange':  '#E66100',
    'crimson': '#DC2626',
    'magenta': '#C51B7D',
}


def render_fig5d_bubble(pathway_csv, out_pdf, ring_color='#E66100'):
    df = pd.read_csv(pathway_csv)
    cancer_labels = {
        'Brain': 'Glioblastoma',
        'Kidney': 'Renal Clear Cell\nCarcinoma',
        'Liver': 'Liver Cancer\n(Hepatoblastoma)',
        'Lung': 'Lung Adenocarcinoma',
        'Eye': 'Retinoblastoma',
        'Pancreas': 'Adenosquamous Carcinoma\n(Pancreas)',
        'Prostate': 'Prostate Small\nCell Carcinoma',
        'Breast': 'Invasive Breast\nCarcinoma',
        'Ovary': 'High-Grade Serous\nOvarian Cancer',
        'Biliary Tract': 'Gallbladder\nCancer',
    }
    df['cancer_label'] = df['context'].map(cancer_labels)
    mask_gb_heme = (df['context'] == 'Biliary Tract') & (df['pathway'] == 'Heme synthesis')
    df.loc[mask_gb_heme, 'evidence_tier'] = 'Uncertain'

    df['neglog10_padj'] = -np.log10(np.maximum(df['fea_fdr'].to_numpy(float), 1e-120))
    df['mcmc_corroborated'] = df['mcmc_fdr05_concordant_pct'] >= 50.0

    cancer_order = [cancer_labels[c] for c in [
        'Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
        'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract'
    ]]
    pw_counts = df.groupby('pathway')['context'].nunique()
    pw_min_p = df.groupby('pathway')['fea_fdr'].min()
    pw_order = sorted(df['pathway'].unique(), key=lambda p: (-pw_counts[p], pw_min_p[p]))

    df['cancer_label'] = pd.Categorical(df['cancer_label'], categories=cancer_order, ordered=True)
    df['pathway'] = pd.Categorical(df['pathway'], categories=pw_order, ordered=True)
    df = df.sort_values(['cancer_label', 'pathway']).reset_index(drop=True)

    palette = {'Pathway-level': '#076d6d', 'Component-level': '#08306b', 'Uncertain': '#91cf60'}
    evidence_order = ['Pathway-level', 'Component-level', 'Uncertain']

    s_min, s_max = 170.0, 700.0
    v_min, v_max = df['neglog10_padj'].min(), df['neglog10_padj'].max()

    sns.set_theme(style='white', context='paper')
    fig, ax = plt.subplots(figsize=(15.0, 8.5))
    ax = sns.scatterplot(
        data=df, x='pathway', y='cancer_label', size='neglog10_padj',
        hue='evidence_tier', hue_order=evidence_order, sizes=(s_min, s_max),
        palette=palette, edgecolor='white', linewidth=1.1, alpha=0.92, ax=ax, zorder=3,
    )

    for _, r in df.iterrows():
        if r['mcmc_corroborated']:
            sz = s_min + (s_max - s_min) * (r['neglog10_padj'] - v_min) / max(v_max - v_min, 1e-9)
            sz_outer = (np.sqrt(sz) + 3.2) ** 2
            ax.scatter([r['pathway']], [r['cancer_label']], s=sz_outer,
                       facecolors='none', edgecolors=ring_color, linewidths=1.2, zorder=5)

    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.22)
    ax.xaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.22)
    ax.set_axisbelow(True)
    sns.despine(left=True, bottom=True)
    ax.set_xlabel('Enriched Metabolic Pathways', fontsize=22, fontweight='bold', labelpad=10)
    ax.set_ylabel('Cancer Type', fontsize=22, fontweight='bold', labelpad=12)
    plt.xticks(rotation=45, ha='right', fontsize=14)
    plt.yticks(fontsize=14)

    handles, labels = ax.get_legend_handles_labels()
    labels = [
        'Evidence level' if l == 'evidence_tier' else
        '$-\\log_{10}(p_{\\mathrm{adj}})$' if l == 'neglog10_padj' else l
        for l in labels
    ]
    handles.append(Line2D([0], [0], marker='o', color='w',
                          label='MCMC corroborated ($\\geq 50\\%$)',
                          markerfacecolor='white', markeredgecolor=ring_color,
                          markeredgewidth=2.5, markersize=14))
    labels.append('MCMC corroborated ($\\geq 50\\%$)')

    leg = ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
                    borderaxespad=0, frameon=True, fontsize=13, title='Evidence & Significance',
                    title_fontsize=13, labelspacing=1.15, handletextpad=0.8)
    for h in leg.legend_handles[1:4]:
        h.set_markersize(13)

    plt.tight_layout()
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ring-color', default='orange', choices=['orange', 'crimson', 'magenta'])
    args = ap.parse_args()
    pw_csv = os.path.join(HERE, 'step4_new_fsr_pathways_with_mcmc_corroboration.csv')
    out_pdf = os.path.join(VIS_DIR, 'bubble_fea_final.pdf')
    render_fig5d_bubble(pw_csv, out_pdf, ring_color=RING_COLORS[args.ring_color])
    print("Saved Figure 5d bubble plot to", out_pdf)
