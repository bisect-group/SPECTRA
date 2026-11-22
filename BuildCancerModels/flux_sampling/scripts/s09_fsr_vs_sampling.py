"""s09_fsr_vs_sampling.py -- Reviewer-2 response, Phases 4 and 5.

Puts the two methods side by side:

  * agreement between the FSR calls (reproduced exactly in s03) and the
    sampling calls (s08), per context and pooled;
  * what sampling says about the specific reactions the reviewer's criticism
    targets -- FSR "up" calls whose span widened because the minimum flux fell
    rather than because the maximum rose;
  * flux enrichment analysis (corrected one-sided hypergeometric, fea.py) on the
    sampling-derived reaction sets, compared with the same test applied to the
    FSR-derived sets, so the pathway-level conclusions can be compared without
    the p-value problem found in s07 confounding it.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fea import fea
from s07_validate_fea import model_subsystems

ROOT = '/home/pavan/scratch/SPECTRA_review_2'
OUT = os.path.join(ROOT, 'results')
CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
CANCER = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
          'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961']


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def main():
    fsr = pd.read_csv(os.path.join(OUT, 'fsr_reproduction_per_reaction.csv'))
    smp = pd.read_csv(os.path.join(OUT, 'sampling_per_reaction.csv'))

    m = fsr.merge(smp, on=['context', 'rxn'], how='inner',
                  suffixes=('_fsr', '_smp'))
    m['fsr_call'] = np.where(m.call_up_fsr, 'up',
                     np.where(m.call_down_fsr, 'down', 'none'))
    m['smp_call'] = np.where(m.call_up_smp, 'up',
                     np.where(m.call_down_smp, 'down', 'none'))
    m['dmax'] = m.max_cancer - m.max_normal
    m['dmin'] = m.min_cancer - m.min_normal
    m['min_driven'] = (-m.dmin) > m.dmax
    m.to_csv(os.path.join(OUT, 'fsr_vs_sampling_per_reaction.csv'), index=False)

    # ---------- per-context agreement ----------
    rows = []
    for ctx in CONTEXTS:
        s = m[m.context == ctx]
        fu = set(s.loc[s.fsr_call == 'up', 'rxn']);  su = set(s.loc[s.smp_call == 'up', 'rxn'])
        fd = set(s.loc[s.fsr_call == 'down', 'rxn']); sd = set(s.loc[s.smp_call == 'down', 'rxn'])
        rows.append({
            'context': ctx, 'n_common': len(s),
            'fsr_up': len(fu), 'smp_up': len(su),
            'up_overlap': len(fu & su), 'up_jaccard': jac(fu, su),
            'fsr_down': len(fd), 'smp_down': len(sd),
            'down_overlap': len(fd & sd), 'down_jaccard': jac(fd, sd),
            'fsr_up_called_down_by_sampling': len(fu & sd),
            'fsr_down_called_up_by_sampling': len(fd & su),
        })
    agree = pd.DataFrame(rows)
    agree.to_csv(os.path.join(OUT, 'fsr_vs_sampling_agreement.csv'), index=False)

    print('=== agreement between FSR calls and sampling calls ===')
    print(agree[['context', 'fsr_up', 'smp_up', 'up_overlap', 'up_jaccard',
                 'fsr_down', 'smp_down', 'down_overlap', 'down_jaccard']].to_string(index=False))
    print('\n=== direct contradictions ===')
    print(agree[['context', 'fsr_up_called_down_by_sampling',
                 'fsr_down_called_up_by_sampling']].to_string(index=False))

    # ---------- confusion matrix, pooled ----------
    print('\n=== pooled confusion (rows FSR, cols sampling) ===')
    print(pd.crosstab(m.fsr_call, m.smp_call).to_string())

    # ---------- what happens to the reviewer's cases ----------
    up = m[m.fsr_call == 'up']
    md = up[up.min_driven]
    print('\n=== FSR "up" reactions, verdict from sampling ===')
    print('all FSR up-calls (%d): %s' % (len(up), up.smp_call.value_counts().to_dict()))
    print('min-driven subset (%d): %s' % (len(md), md.smp_call.value_counts().to_dict()))
    if len(up):
        print('median |Cliff delta| : all up-calls %.3f | min-driven %.3f'
              % (up.cliffs_delta_abs.abs().median(),
                 md.cliffs_delta_abs.abs().median() if len(md) else float('nan')))

    # ---------- does FSR magnitude predict the real effect? ----------
    v = m[(m.fsr > 0) & np.isfinite(m.fsr)]
    if len(v) > 10:
        r = np.corrcoef(np.log2(v.fsr), v.cliffs_delta_abs)[0, 1]
        print('\nSpearman-like check: corr(log2 FSR, Cliff delta on |v|) = %.3f over %d reactions'
              % (r, len(v)))

    # ---------- enrichment on both sets, same corrected test ----------
    fea_rows = []
    for ctx, ach in zip(CONTEXTS, CANCER):
        s = m[m.context == ctx]
        rxns, ss = model_subsystems(ach)
        pos = {r: i for i, r in enumerate(rxns)}
        for label, sel in (('FSR', s.fsr_call != 'none'), ('sampling', s.smp_call != 'none')):
            q = [pos[r] for r in s.loc[sel, 'rxn'] if r in pos]
            if not q:
                continue
            for g in fea(ss, q):
                g.update({'context': ctx, 'method': label, 'query_size': len(q)})
                fea_rows.append(g)
    fe = pd.DataFrame(fea_rows)
    fe.to_csv(os.path.join(OUT, 'fea_fsr_vs_sampling.csv'), index=False)

    sig = fe[fe.p_adj < 0.05]
    print('\n=== enriched pathways (FDR<0.05, corrected one-sided test) ===')
    rows = []
    for ctx in CONTEXTS:
        a = set(sig[(sig.context == ctx) & (sig.method == 'FSR')].group)
        b = set(sig[(sig.context == ctx) & (sig.method == 'sampling')].group)
        rows.append({'context': ctx, 'fsr_pathways': len(a), 'sampling_pathways': len(b),
                     'shared': len(a & b), 'jaccard': jac(a, b),
                     'fsr_only': '; '.join(sorted(a - b))[:70],
                     'sampling_only': '; '.join(sorted(b - a))[:70]})
    pw = pd.DataFrame(rows)
    pw.to_csv(os.path.join(OUT, 'pathway_overlap.csv'), index=False)
    print(pw[['context', 'fsr_pathways', 'sampling_pathways', 'shared', 'jaccard']].to_string(index=False))

    # pan-cancer signature: pathways enriched in >=5 contexts, each method
    print('\n=== pan-cancer signature (enriched in >= 5 of 10 contexts) ===')
    for label in ('FSR', 'sampling'):
        cnt = sig[sig.method == label].groupby('group').context.nunique().sort_values(ascending=False)
        keep = cnt[cnt >= 5]
        print('%-9s %d pathways: %s' % (label, len(keep), ', '.join(keep.index[:12])))


if __name__ == '__main__':
    main()
