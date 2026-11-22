"""s10_figures.py -- Reviewer-2 response figures.

  fig1_convergence      chain diagnostics (R-hat, ESS) -- evidence the MCMC mixed
  fig2_fsr_vs_delta     FSR against the sampled effect size, the core comparison
  fig3_failure_mode     worked examples of reactions FSR and sampling disagree on
  fig4_pathways         pathway enrichment, FSR-derived vs sampling-derived sets
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = '/home/pavan/scratch/SPECTRA_review_2'
OUT = os.path.join(ROOT, 'results')
FIG = os.path.join(ROOT, 'figures')
CHAINS = os.path.join(ROOT, 'sampling', 'chains')
CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
PAIRS = dict(zip(CONTEXTS, zip(
    ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
     'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961'],
    ['ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466',
     'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854'])))

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 150})


def fig_convergence(m):
    fig, ax = plt.subplots(1, 3, figsize=(9, 2.7))
    rh = np.r_[m.rhat_cancer.values, m.rhat_normal.values]
    rh = rh[np.isfinite(rh)]
    ax[0].hist(rh, bins=np.linspace(1, max(1.2, np.percentile(rh, 99.9)), 60),
               color='#4C72B0')
    ax[0].axvline(1.05, color='crimson', lw=1, ls='--', label='1.05')
    ax[0].set_xlabel(r'Gelman-Rubin $\hat{R}$'); ax[0].set_ylabel('reactions')
    ax[0].set_title(r'$\hat{R}$ across 4 chains'); ax[0].legend(frameon=False)

    ess = np.r_[m.ess_cancer.values, m.ess_normal.values]
    ax[1].hist(ess[np.isfinite(ess)], bins=60, color='#55A868')
    ax[1].set_xlabel('effective sample size'); ax[1].set_ylabel('reactions')
    ax[1].set_title('ESS (4 chains pooled)')

    by = m.groupby('context').agg(maxr=('rhat_cancer', 'max'),
                                  medess=('ess_cancer', 'median'))
    ax[2].barh(range(len(by)), by.medess.values, color='#C44E52')
    ax[2].set_yticks(range(len(by))); ax[2].set_yticklabels(by.index, fontsize=6)
    ax[2].set_xlabel('median ESS'); ax[2].set_title('per context')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig1_convergence.png'))
    plt.close(fig)


def fig_fsr_vs_delta(m):
    fig, axes = plt.subplots(2, 5, figsize=(13, 5.4), sharex=True, sharey=True)
    for axx, ctx in zip(axes.ravel(), CONTEXTS):
        s = m[(m.context == ctx) & np.isfinite(m.fsr) & (m.fsr > 0)]
        x = np.log2(s.fsr.values)
        y = s.cliffs_delta_abs.values
        axx.scatter(x, y, s=1.5, alpha=0.15, color='0.5', rasterized=True)
        for lab, col in (('up', '#C44E52'), ('down', '#4C72B0')):
            q = s[s.fsr_call == lab]
            axx.scatter(np.log2(q.fsr.values), q.cliffs_delta_abs.values,
                        s=4, alpha=0.7, color=col, label='FSR ' + lab, rasterized=True)
        axx.axhline(0.33, color='k', lw=0.6, ls=':')
        axx.axhline(-0.33, color='k', lw=0.6, ls=':')
        axx.axvline(1, color='k', lw=0.6, ls='--')   # FSR = 2
        axx.axvline(-1, color='k', lw=0.6, ls='--')  # FSR = 0.5
        axx.set_title(ctx, fontsize=8)
    axes[0, 0].legend(frameon=False, fontsize=6, markerscale=2)
    fig.supxlabel(r'$\log_2$ flux span ratio', fontsize=9)
    fig.supylabel(r"Cliff's $\delta$ on sampled $|v|$ (cancer vs normal)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig2_fsr_vs_delta.png'))
    plt.close(fig)


def fig_failure_mode(m):
    """Flux distributions for FSR up-calls that sampling does not support."""
    # The reviewer's case: FSR called these up because the span widened, but the
    # widening came from the minimum flux falling. Sampling shows what actually
    # happened -- the reaction runs in reverse, it is not carrying more forward flux.
    cand = m[m.min_driven & (m.fsr_call == 'up') &
             (m.median_signed_cancer < m.median_signed_normal)].copy()
    if len(cand) < 6:
        extra = m[(m.fsr_call == 'up') & (m.smp_call != 'up')].copy()
        extra['rank'] = extra.fsr / (extra.cliffs_delta_abs.abs() + 0.01)
        cand = pd.concat([cand, extra.sort_values('rank', ascending=False)])
    if cand.empty:
        print('  (nothing to show; skipping fig3)')
        return
    cand['score'] = (cand.median_signed_normal - cand.median_signed_cancer).abs()
    cand = cand.sort_values('score', ascending=False).head(6)

    fig, axes = plt.subplots(2, 3, figsize=(10, 5))
    for axx, (_, r) in zip(axes.ravel(), cand.iterrows()):
        ach_c, ach_n = PAIRS[r.context]
        vc = _load_rxn(ach_c, r.rxn)
        vn = _load_rxn(ach_n, r.rxn)
        if vc is None or vn is None:
            continue
        lo = min(np.percentile(vn, 0.5), np.percentile(vc, 0.5))
        hi = max(np.percentile(vn, 99.5), np.percentile(vc, 99.5))
        bins = np.linspace(lo, hi, 60) if hi > lo else 30
        axx.hist(vn, bins=bins, alpha=0.55, label='normal', color='#4C72B0', density=True)
        axx.hist(vc, bins=bins, alpha=0.55, label='cancer', color='#C44E52', density=True)
        axx.axvline(r.min_normal, color='#4C72B0', lw=0.8, ls='--')
        axx.axvline(r.max_normal, color='#4C72B0', lw=0.8, ls='--')
        axx.axvline(r.min_cancer, color='#C44E52', lw=0.8, ls=':')
        axx.axvline(r.max_cancer, color='#C44E52', lw=0.8, ls=':')
        axx.set_title('%s / %s\nFSR=%.1f  $\\delta_{|v|}$=%.2f   median $v$: %.1f$\\to$%.1f'
                      % (r.context, r.rxn, r.fsr, r.cliffs_delta_abs,
                         r.median_signed_normal, r.median_signed_cancer),
                      fontsize=7)
        axx.set_yticks([])
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle('FSR calls these upregulated because the flux span widened downward.\n'
                 'The sampled distributions show the reaction running in reverse, not carrying more flux.\n'
                 '(dashed = normal FVA min/max, dotted = cancer FVA min/max: all FSR sees)', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig3_failure_mode.png'))
    plt.close(fig)


_cache = {}


def _load_rxn(ach, rxn):
    if ach not in _cache:
        mats, rxns = [], None
        for c in range(1, 5):
            f = os.path.join(CHAINS, '%s_chain%d.mat' % (ach, c))
            if not os.path.exists(f):
                return None
            d = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
            mats.append(np.asarray(d['samples'], dtype=np.float32))
            rxns = [str(x) for x in np.atleast_1d(d['rxns'])]
        _cache[ach] = (rxns, np.hstack(mats))
    rxns, X = _cache[ach]
    if rxn not in rxns:
        return None
    return X[rxns.index(rxn), :]


def fig_pathways(fe):
    sig = fe[fe.p_adj < 0.05]
    groups = sorted(set(sig.group))
    if not groups:
        print('  (no significant pathways; skipping fig4)')
        return
    Mx = np.zeros((len(groups), len(CONTEXTS) * 2))
    for gi, g in enumerate(groups):
        for ci, ctx in enumerate(CONTEXTS):
            for mi, meth in enumerate(('FSR', 'sampling')):
                q = sig[(sig.group == g) & (sig.context == ctx) & (sig.method == meth)]
                if len(q):
                    Mx[gi, ci * 2 + mi] = -np.log10(max(q.p_adj.iloc[0], 1e-30))
    fig, ax = plt.subplots(figsize=(9, max(3, 0.17 * len(groups))))
    im = ax.imshow(Mx, aspect='auto', cmap='magma_r')
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([g[:42] for g in groups], fontsize=5.5)
    ax.set_xticks(np.arange(0.5, len(CONTEXTS) * 2, 2))
    ax.set_xticklabels(CONTEXTS, rotation=45, ha='right', fontsize=7)
    for ci in range(1, len(CONTEXTS)):
        ax.axvline(ci * 2 - 0.5, color='w', lw=0.8)
    ax.set_title('Enriched subsystems, FDR < 0.05 (left bar = FSR set, right bar = sampling set)',
                 fontsize=8)
    fig.colorbar(im, ax=ax, label=r'$-\log_{10}$ FDR', shrink=0.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'fig4_pathways.png'))
    plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    m = pd.read_csv(os.path.join(OUT, 'fsr_vs_sampling_per_reaction.csv'))
    print('fig1 convergence');  fig_convergence(m)
    print('fig2 FSR vs delta'); fig_fsr_vs_delta(m)
    print('fig3 failure mode'); fig_failure_mode(m)
    fp = os.path.join(OUT, 'fea_fsr_vs_sampling.csv')
    if os.path.exists(fp):
        print('fig4 pathways'); fig_pathways(pd.read_csv(fp))
    print('figures written to', FIG)


if __name__ == '__main__':
    main()
