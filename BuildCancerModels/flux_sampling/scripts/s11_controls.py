"""s11_controls.py -- two controls the headline comparison needs.

CONTROL 1 -- DIRECTION vs MAGNITUDE.
s08 calls dysregulation on |v| (how much flux the reaction carries). For a
reversible reaction, |v| can rise because the reaction now runs hard in REVERSE,
which is a change of direction, not upregulation. This splits every call by
whether the signed flux moved the same way, so the two are not conflated -- and
in particular checks the FSR "min-driven" reactions, whose minimum flux fell,
since those are exactly the ones at risk of being reverse-flux cases.

CONTROL 2 -- SET-SIZE MATCHED ENRICHMENT.
Sampling calls 5-10x more reactions than FSR, and a hypergeometric test on a
larger query set has more power, so "sampling finds more pathways" is partly an
artefact of set size. This re-runs the enrichment on the top-N sampling
reactions by |Cliff's delta|, with N set to the FSR set size for that context,
so the two pathway lists are comparable.
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


def main():
    m = pd.read_csv(os.path.join(OUT, 'fsr_vs_sampling_per_reaction.csv'))

    # ---------------- Control 1 ----------------
    print('=== CONTROL 1: is an increase in |v| an increase, or a reversal? ===')
    up = m[m.smp_call == 'up'].copy()
    # signed median moved up as well -> genuine increase in the forward direction;
    # signed median moved down while |v| rose -> reaction is running in reverse
    up['signed_change'] = up.median_signed_cancer - up.median_signed_normal
    up['reversal'] = (up.median_signed_cancer * up.median_signed_normal < 0)
    print('sampling "up" calls: %d total' % len(up))
    print('  signed median also increased : %d (%.1f%%)'
          % ((up.signed_change > 0).sum(), 100 * (up.signed_change > 0).mean()))
    print('  signed median flipped sign   : %d (%.1f%%)'
          % (up.reversal.sum(), 100 * up.reversal.mean()))

    fu = m[m.fsr_call == 'up']
    md = fu[fu.min_driven]
    print('\nFSR up-calls (%d) and the "min-driven" subset (%d):' % (len(fu), len(md)))
    for lab, s in (('all FSR up', fu), ('min-driven', md)):
        conf = s[s.smp_call == 'up']
        rev = (conf.median_signed_cancer * conf.median_signed_normal < 0).sum() if len(conf) else 0
        neg = (conf.median_signed_cancer < conf.median_signed_normal).sum() if len(conf) else 0
        print('  %-11s n=%4d | sampling says up %4d | of those, signed median DECREASED %4d (%.0f%%),'
              ' sign flipped %3d'
              % (lab, len(s), len(conf), neg, 100 * neg / max(len(conf), 1), rev))
    print('\n  -> a "min-driven" FSR up-call is largely a reaction running harder in')
    print('     REVERSE, which |v| scores as more active but is a directional change.')

    # median |delta| by group, for the write-up
    print('\nmedian |Cliff delta| on |v| : all FSR up %.3f | min-driven %.3f | sampling-only up %.3f'
          % (fu.cliffs_delta_abs.abs().median(),
             md.cliffs_delta_abs.abs().median() if len(md) else np.nan,
             m[(m.smp_call == 'up') & (m.fsr_call == 'none')].cliffs_delta_abs.abs().median()))

    # ---------------- Control 2 ----------------
    print('\n=== CONTROL 2: enrichment with set sizes matched to FSR ===')
    rows, fea_rows = [], []
    for ctx, ach in zip(CONTEXTS, CANCER):
        s = m[m.context == ctx]
        n_fsr = int((s.fsr_call != 'none').sum())
        rxns, ss = model_subsystems(ach)
        pos = {r: i for i, r in enumerate(rxns)}

        sig = s[(s.smp_call != 'none')].copy()
        sig['absd'] = sig.cliffs_delta_abs.abs()
        top = sig.nlargest(n_fsr, 'absd')

        sets = {
            'FSR': list(s.loc[s.fsr_call != 'none', 'rxn']),
            'sampling_all': list(sig.rxn),
            'sampling_topN': list(top.rxn),
        }
        pw = {}
        for label, rl in sets.items():
            q = [pos[r] for r in rl if r in pos]
            res = fea(ss, q) if q else []
            pw[label] = {g['group'] for g in res if g['p_adj'] < 0.05}
            for g in res:
                g.update({'context': ctx, 'method': label, 'query_size': len(q)})
                fea_rows.append(g)
        a, b = pw['FSR'], pw['sampling_topN']
        rows.append({'context': ctx, 'n_query': n_fsr,
                     'fsr_pathways': len(a), 'sampling_topN_pathways': len(b),
                     'shared': len(a & b),
                     'jaccard': len(a & b) / max(len(a | b), 1),
                     'sampling_all_pathways': len(pw['sampling_all'])})
    pw = pd.DataFrame(rows)
    pw.to_csv(os.path.join(OUT, 'pathway_overlap_size_matched.csv'), index=False)
    fe = pd.DataFrame(fea_rows)
    fe.to_csv(os.path.join(OUT, 'fea_size_matched.csv'), index=False)
    print(pw.to_string(index=False))

    sig = fe[fe.p_adj < 0.05]
    print('\n=== pan-cancer signature (>= 5 of 10 contexts), size-matched ===')
    for label in ('FSR', 'sampling_topN', 'sampling_all'):
        cnt = sig[sig.method == label].groupby('group').context.nunique().sort_values(ascending=False)
        keep = cnt[cnt >= 5]
        print('%-14s %2d pathways' % (label, len(keep)))
        for g, c in keep.items():
            print('                 %-52s %d/10' % (g[:52], c))


if __name__ == '__main__':
    main()
