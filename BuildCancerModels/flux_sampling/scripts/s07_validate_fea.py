"""s07_validate_fea.py -- checkpoint on the enrichment step.

Runs the Python re-implementation of COBRA's FEA (fea.py) on the FSR-derived
dysregulated reaction sets reproduced in s03, and compares against the published
BuildCancerModels/flux_results/fea_2025_checked_2.xlsx.

Two separate things are checked:

  (a) COMPOSITION -- the set of subsystems hit, the number of query reactions in
      each (k), and the number of model reactions in each (K). These must match,
      and they do; that confirms the query sets and the background model are the
      same ones the paper used.

  (b) P-VALUES -- these do NOT match. The published values cannot be produced by
      the documented statistic (upper-tail hypergeometric on M = numel(model.rxns),
      N = numel(rxnSet)) under any single parameterisation: solving each group for
      the M or N that would yield the published p gives inconsistent answers, and
      'Exchange/demand reaction' -- hit by 3 of its 1103 reactions, i.e. strongly
      DEPLETED relative to the 25.5 expected -- is reported at p = 2.4e-9, which no
      upper-tail enrichment test can produce.

So this script reports the consequence instead of asserting equality: which
pathways pass FDR < 0.05 under the published numbers versus under the correct
one-sided hypergeometric test.
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xlsxio
from fea import fea

ROOT = '/home/pavan/scratch/SPECTRA_review_2'
REPO = '/home/pavan/scratch/spectra/BuildCancerModels'
CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
CANCER = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
          'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961']


def model_subsystems(ach):
    d = sio.loadmat(os.path.join(ROOT, 'models', 'sprintcore', ach + '.mat'),
                    squeeze_me=True, struct_as_record=False)
    m = d['m']
    rxns = [str(x) for x in np.atleast_1d(m.rxns)]
    ss = []
    for x in np.atleast_1d(m.subSystems):
        ss.append(' | '.join(str(y) for y in np.atleast_1d(x))
                  if isinstance(x, np.ndarray) else str(x))
    return rxns, ss


def main():
    fsr = pd.read_csv(os.path.join(ROOT, 'results', 'fsr_reproduction_per_reaction.csv'))
    pub = xlsxio.read_sheets(os.path.join(REPO, 'flux_results', 'fea_2025_checked_2.xlsx'))

    comp_ok = True
    rows_out, detail = [], []
    for ctx, ach in zip(CONTEXTS, CANCER):
        sub = fsr[fsr.context == ctx]
        query = list(sub.loc[sub.call_down, 'rxn']) + list(sub.loc[sub.call_up, 'rxn'])

        rxns, ss = model_subsystems(ach)
        pos = {r: i for i, r in enumerate(rxns)}
        got = fea(ss, [pos[r] for r in query])
        gmap = {g['group']: g for g in got}

        prows = [r for r in pub[ctx][2:] if r and r[2]]
        pmap = {r[2]: r for r in prows}

        same_groups = set(gmap) == set(pmap)
        same_k = all(gmap[g]['enriched_set_size'] == int(pmap[g][3]) for g in pmap if g in gmap)
        same_K = all(gmap[g]['total_set_size'] == int(pmap[g][4]) for g in pmap if g in gmap)
        comp_ok &= same_groups and same_k and same_K

        pub_sig = {g for g, r in pmap.items() if float(r[1]) < 0.05}
        new_sig = {g['group'] for g in got if g['p_adj'] < 0.05}
        # groups the published run calls enriched while they are actually depleted
        depleted = {g for g in pub_sig
                    if gmap[g]['enriched_set_size']
                    < len(query) * gmap[g]['total_set_size'] / len(rxns)}

        rows_out.append({
            'context': ctx, 'query_size': len(query),
            'groups_hit': len(got), 'same_groups': same_groups,
            'same_k': same_k, 'same_K': same_K,
            'pub_sig': len(pub_sig), 'correct_sig': len(new_sig),
            'lost': len(pub_sig - new_sig), 'gained': len(new_sig - pub_sig),
            'pub_sig_but_depleted': len(depleted),
        })
        for g in sorted(pub_sig - new_sig):
            e = len(query) * gmap[g]['total_set_size'] / len(rxns)
            detail.append({'context': ctx, 'group': g, 'status': 'published-only',
                           'k': gmap[g]['enriched_set_size'],
                           'K': gmap[g]['total_set_size'], 'expected': round(e, 2),
                           'pub_p_adj': float(pmap[g][1]), 'correct_p_adj': gmap[g]['p_adj'],
                           'direction': 'DEPLETED' if gmap[g]['enriched_set_size'] < e else 'enriched'})
        for g in sorted(new_sig - pub_sig):
            e = len(query) * gmap[g]['total_set_size'] / len(rxns)
            detail.append({'context': ctx, 'group': g, 'status': 'correct-only',
                           'k': gmap[g]['enriched_set_size'],
                           'K': gmap[g]['total_set_size'], 'expected': round(e, 2),
                           'pub_p_adj': float(pmap[g][1]), 'correct_p_adj': gmap[g]['p_adj'],
                           'direction': 'DEPLETED' if gmap[g]['enriched_set_size'] < e else 'enriched'})

    sm = pd.DataFrame(rows_out)
    dt = pd.DataFrame(detail)
    sm.to_csv(os.path.join(ROOT, 'results', 'fea_published_vs_correct_summary.csv'), index=False)
    dt.to_csv(os.path.join(ROOT, 'results', 'fea_published_vs_correct_detail.csv'), index=False)

    print('=== (a) composition: does the query set / background reproduce? ===')
    print(sm[['context', 'query_size', 'groups_hit', 'same_groups', 'same_k', 'same_K']]
          .to_string(index=False))
    print('\nPASS: subsystem membership, k and K reproduce exactly for all 10 contexts.'
          if comp_ok else '\nFAIL: composition does not reproduce.')

    print('\n=== (b) p-values: pathways passing FDR < 0.05 ===')
    print(sm[['context', 'pub_sig', 'correct_sig', 'lost', 'gained', 'pub_sig_but_depleted']]
          .to_string(index=False))
    print('\ntotals: published %d significant, correct one-sided test %d significant; '
          '%d published calls are of DEPLETED subsystems'
          % (sm.pub_sig.sum(), sm.correct_sig.sum(), sm.pub_sig_but_depleted.sum()))
    if len(dt):
        print('\nlargest published-only calls that are actually depleted:')
        d2 = dt[(dt.status == 'published-only') & (dt.direction == 'DEPLETED')]
        print(d2.sort_values('K', ascending=False)
                .head(12)[['context', 'group', 'k', 'K', 'expected', 'pub_p_adj']]
                .to_string(index=False))


if __name__ == '__main__':
    main()
