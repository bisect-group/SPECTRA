"""s03_fsr_diagnostics.py  --  Reviewer-2 response, Phase 5 (part 1).

Reproduces the published Flux Span Ratio (FSR) calls from the archived FVA
sheets, then quantifies the specific failure mode the reviewer raised:

    "The FSR can be mathematically large simply due to a much lower minimum
     flux, even if the maximum flux is slightly higher than in the normal
     state."

Needs no metabolic models -- only BuildCancerModels/fva/*.xlsx -- so it runs
independently of the model rebuild.

Implemented FSR (from fsranalysis_ccle_sprint.m):

    if  max_c != max_n  and  min_c <= min_n :
        FSR = (max_c - min_c) / (max_n - min_n)
    elif max_c == max_n and min_c == min_n :
        FSR = 0
    else:
        (never assigned -- see NOTE below)

    up   : FSR != 0, finite, 2    <= FSR < 20000
    down : FSR != 0, finite, 0.01 <= FSR <= 0.5

NOTE: reactions falling in the third branch are left unassigned. In MATLAB an
unassigned element of a growing array is zero-filled, so they are dropped by
the `FSR != 0` filter -- i.e. any reaction whose minimum flux *increased* in
cancer is silently excluded from consideration, and any such reactions at the
tail of the list are dropped by the resulting length mismatch. Both effects
are reproduced here so the comparison against the published output is exact.
"""
import os
import sys
import json

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xlsxio

REPO = '/home/pavan/scratch/spectra/BuildCancerModels'
ROOT = '/home/pavan/scratch/SPECTRA_review_2'
FVA_DIR = os.path.join(REPO, 'fva')
OUT_DIR = os.path.join(ROOT, 'results')

CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
CANCER = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
          'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961']
NORMAL = ['ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466',
          'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854']


def load_fva(ach):
    rows = xlsxio.read_sheets(os.path.join(FVA_DIR, ach + '.xlsx'))['Sheet1']
    recs = [(r[0], r[1], r[2]) for r in rows[1:]
            if r and r[0] not in (None, '')]
    df = pd.DataFrame(recs, columns=['rxn', 'min', 'max'])
    df['min'] = pd.to_numeric(df['min'], errors='coerce').fillna(0.0)
    df['max'] = pd.to_numeric(df['max'], errors='coerce').fillna(0.0)
    return df.set_index('rxn')


def fsr_for_pair(nrm, can):
    """Replicate fsranalysis_ccle_sprint.m for one cancer/normal pair."""
    common = sorted(set(nrm.index) & set(can.index))          # MATLAB intersect sorts
    common = [r for r in common
              if not (r.startswith('EX_') or r.startswith('sink_')
                      or r.startswith('DM_'))]

    max_n = nrm.loc[common, 'max'].to_numpy(float)
    min_n = nrm.loc[common, 'min'].to_numpy(float)
    max_c = can.loc[common, 'max'].to_numpy(float)
    min_c = can.loc[common, 'min'].to_numpy(float)

    n = len(common)
    fsr = np.zeros(n)
    assigned = np.zeros(n, dtype=bool)      # branch 1 or 2 fired
    last = -1
    with np.errstate(divide='ignore', invalid='ignore'):
        for i in range(n):
            if max_c[i] != max_n[i] and min_c[i] <= min_n[i]:
                fsr[i] = (max_c[i] - min_c[i]) / (max_n[i] - min_n[i])
                assigned[i] = True
                last = i
            elif max_c[i] == max_n[i] and min_c[i] == min_n[i]:
                fsr[i] = 0.0
                assigned[i] = True
                last = i

    # MATLAB leaves the array at length `last+1`; logical indexing past that
    # point silently truncates, so reactions beyond it are never called.
    in_range = np.zeros(n, dtype=bool)
    in_range[:last + 1] = True

    ok = in_range & (fsr != 0) & np.isfinite(fsr)
    up = ok & (fsr >= 2) & (fsr < 20000.0)
    down = ok & (fsr <= 0.5) & (fsr >= 0.01)

    return pd.DataFrame({
        'rxn': common,
        'min_normal': min_n, 'max_normal': max_n,
        'min_cancer': min_c, 'max_cancer': max_c,
        'span_normal': max_n - min_n,
        'span_cancer': max_c - min_c,
        'fsr': fsr,
        'assigned': assigned,
        'in_range': in_range,
        'call_up': up,
        'call_down': down,
    })


def published_calls(fname, sheet):
    sheets = xlsxio.read_sheets(os.path.join(REPO, 'flux_results', fname))
    rows = sheets[sheet]
    return {r[0] for r in rows[1:] if r and r[0]}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    fva = {a: load_fva(a) for a in set(CANCER) | set(NORMAL)}

    per_rxn, summary = [], []
    for ctx, c, nn in zip(CONTEXTS, CANCER, NORMAL):
        df = fsr_for_pair(fva[nn], fva[c])
        df.insert(0, 'context', ctx)
        per_rxn.append(df)

        pub_up = published_calls('fsr_up_2025_checked_2.xlsx', ctx)
        pub_dn = published_calls('fsr_down_2025_checked_2.xlsx', ctx)
        got_up = set(df.loc[df.call_up, 'rxn'])
        got_dn = set(df.loc[df.call_down, 'rxn'])

        up = df[df.call_up]
        dmax = up.max_cancer - up.max_normal
        dmin = up.min_cancer - up.min_normal

        # A call is "min-driven" when the widening of the span comes mostly from
        # the lower bound dropping rather than the upper bound rising.
        min_driven = (-dmin) > dmax
        # The paper's stated criterion is max_cancer > max_normal. Reactions
        # failing it are called up on span alone.
        violates_paper = up.max_cancer <= up.max_normal
        # Max essentially unchanged (<=5% of the normal max) yet called up.
        denom = np.maximum(np.abs(up.max_normal), 1e-9)
        max_flat = (np.abs(dmax) / denom) <= 0.05

        summary.append({
            'context': ctx, 'cancer': c, 'normal': nn,
            'n_common': len(df),
            'n_dropped_min_increased': int((~df.assigned).sum()),
            'n_dropped_truncation': int((~df.in_range).sum()),
            'n_up': len(got_up), 'n_down': len(got_dn),
            'pub_n_up': len(pub_up), 'pub_n_down': len(pub_dn),
            'up_match': got_up == pub_up, 'down_match': got_dn == pub_dn,
            'up_jaccard': len(got_up & pub_up) / max(len(got_up | pub_up), 1),
            'down_jaccard': len(got_dn & pub_dn) / max(len(got_dn | pub_dn), 1),
            'up_min_driven': int(min_driven.sum()),
            'up_min_driven_pct': 100 * min_driven.mean() if len(up) else np.nan,
            'up_violates_paper_criterion': int(violates_paper.sum()),
            'up_violates_pct': 100 * violates_paper.mean() if len(up) else np.nan,
            'up_max_flat_le5pct': int(max_flat.sum()),
            'up_max_flat_pct': 100 * max_flat.mean() if len(up) else np.nan,
        })
        print(f'{ctx:15s} common={len(df):5d} up={len(got_up):4d} (published {len(pub_up):4d}) '
              f'down={len(got_dn):4d} (published {len(pub_dn):4d})  '
              f'match up/down = {got_up == pub_up}/{got_dn == pub_dn}')

    allrxn = pd.concat(per_rxn, ignore_index=True)
    allrxn.to_csv(os.path.join(OUT_DIR, 'fsr_reproduction_per_reaction.csv'), index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(OUT_DIR, 'fsr_reproduction_summary.csv'), index=False)

    print('\n=== reproduction of the published calls ===')
    print(sm[['context', 'n_up', 'pub_n_up', 'up_match', 'up_jaccard',
              'n_down', 'pub_n_down', 'down_match', 'down_jaccard']].to_string(index=False))
    print('\n=== reviewer-2 failure mode, among reactions called UPREGULATED ===')
    print(sm[['context', 'n_up', 'up_min_driven', 'up_min_driven_pct',
              'up_violates_paper_criterion', 'up_violates_pct',
              'up_max_flat_le5pct', 'up_max_flat_pct']].to_string(index=False))
    print('\n=== reactions never evaluated at all ===')
    print(sm[['context', 'n_common', 'n_dropped_min_increased',
              'n_dropped_truncation']].to_string(index=False))

    with open(os.path.join(OUT_DIR, 'fsr_reproduction_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2, default=str)


if __name__ == '__main__':
    main()
