"""s08_compare_distributions.py -- Reviewer-2 response, Phase 3.

Compares the sampled flux distributions of every shared reaction between each
cancer CSM and its non-cancerous counterpart, which is what Reviewer 2 asked
for in place of the flux span ratio.

Design decisions worth stating explicitly, since they are what make this more
defensible than the FSR:

1. TEST STATISTIC. Two-sided Mann-Whitney U on the reaction's sampled flux.
   Rank-based, so it needs no distributional assumption -- flux marginals from a
   polytope are typically skewed, bounded and often multi-modal.

2. WHAT IS COMPARED. Primary analysis uses |v|, the magnitude of flux carried,
   because "up-/down-regulated reaction" is a statement about how much flux the
   reaction carries, and a reversible reaction running harder in reverse is more
   active, not less. The signed-flux comparison is computed alongside as a
   secondary column so directionality changes remain visible.

3. SAMPLE SIZE. With 10,000 MCMC draws per model, a Mann-Whitney p-value is
   significant for arbitrarily tiny shifts -- the same "large number that means
   little" trap the FSR falls into. Two guards:
     (a) each reaction is thinned to min(ESS_cancer, ESS_normal) draws before
         testing, so the nominal n never exceeds the information actually
         present in the chains;
     (b) significance ALONE never makes a call. A reaction is called
         dysregulated only if it also clears an effect-size threshold,
         Cliff's delta |d| >= 0.33 (Romano's "medium"). Cliff's delta is the
         rank equivalent of "how often is cancer flux above normal flux", so it
         is on the same scale for every reaction regardless of flux magnitude.

4. THRESHOLD SENSITIVITY. Calls are recomputed over a grid of delta cutoffs and
   written out, so the reaction sets can be shown not to hinge on 0.33.

Reaction universe is identical to the FSR analysis (reactions shared by both
models, excluding EX_/sink_/DM_) so the two methods are compared like for like.
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fea import bh_fdr

ROOT = '/home/pavan/scratch/SPECTRA_review_2'
CHAINS = os.path.join(ROOT, 'sampling', 'chains')
OUT = os.path.join(ROOT, 'results')

CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
CANCER = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
          'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961']
NORMAL = ['ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466',
          'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854']

NCHAIN = 4
DELTA_CUT = 0.33
DELTA_GRID = [0.147, 0.33, 0.474, 0.6]   # Romano small / medium / large, plus stricter
FDR_CUT = 0.05
CHUNK = 400


def load_chains(ach):
    """Return (rxns, samples nRxn x (NCHAIN*N), per-chain sample counts)."""
    mats, rxns, counts = [], None, []
    for c in range(1, NCHAIN + 1):
        f = os.path.join(CHAINS, '%s_chain%d.mat' % (ach, c))
        if not os.path.exists(f):
            raise FileNotFoundError(f)
        d = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
        s = np.asarray(d['samples'], dtype=np.float32)
        r = [str(x) for x in np.atleast_1d(d['rxns'])]
        if rxns is None:
            rxns = r
        elif r != rxns:
            raise ValueError('reaction order differs between chains of %s' % ach)
        mats.append(s)
        counts.append(s.shape[1])
    return rxns, np.hstack(mats), counts


def split_rhat(mats):
    """Gelman-Rubin R-hat per reaction across chains (equal-length truncation)."""
    n = min(m.shape[1] for m in mats)
    x = np.stack([m[:, :n] for m in mats], axis=0).astype(np.float64)  # chains x rxn x n
    m_chain = x.mean(axis=2)
    v_chain = x.var(axis=2, ddof=1)
    W = v_chain.mean(axis=0)
    B = n * m_chain.var(axis=0, ddof=1)
    var_hat = (n - 1) / n * W + B / n
    with np.errstate(divide='ignore', invalid='ignore'):
        rh = np.sqrt(var_hat / W)
    rh[~np.isfinite(rh)] = 1.0          # constant reactions: no between-chain variance
    return rh


def ess_1d(x):
    """Effective sample size of a single chain via initial positive sequence."""
    n = x.size
    x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return float(n)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n] / (n * v)
    s, t = 0.0, 1
    while t + 1 < n:
        pair = ac[t] + ac[t + 1]
        if pair <= 0:
            break
        s += pair
        t += 2
    return float(min(n, n / (1 + 2 * s)))


def chain_ess(mats, idx):
    """Summed across-chain ESS for the reactions in idx."""
    out = np.zeros(len(idx))
    for m in mats:
        for j, i in enumerate(idx):
            out[j] += ess_1d(m[i, :].astype(np.float64))
    return out


def thin_to(x, k):
    """Evenly spaced k draws from x (columns)."""
    n = x.shape[1]
    if k >= n:
        return x
    return x[:, np.linspace(0, n - 1, k).astype(int)]


def compare_pair(ctx, ach_c, ach_n):
    rxn_c, Xc, cnt_c = load_chains(ach_c)
    rxn_n, Xn, cnt_n = load_chains(ach_n)
    # chains stop at different lengths (RHMC runs until ESS >= N), so split on the
    # recorded per-chain counts -- an equal split would straddle chain boundaries
    # and silently corrupt the between-chain R-hat.
    mats_c = np.split(Xc, np.cumsum(cnt_c)[:-1], axis=1)
    mats_n = np.split(Xn, np.cumsum(cnt_n)[:-1], axis=1)
    assert len(mats_c) == NCHAIN and len(mats_n) == NCHAIN
    assert [m.shape[1] for m in mats_c] == list(cnt_c)

    rhat_c = split_rhat(mats_c)
    rhat_n = split_rhat(mats_n)

    common = sorted(set(rxn_c) & set(rxn_n))
    common = [r for r in common
              if not (r.startswith('EX_') or r.startswith('sink_') or r.startswith('DM_'))]
    ic = np.array([rxn_c.index(r) for r in common])
    inn = np.array([rxn_n.index(r) for r in common])

    ess_c = chain_ess(mats_c, ic)
    ess_n = chain_ess(mats_n, inn)
    n_eff = np.minimum(ess_c, ess_n).astype(int)
    n_eff = np.clip(n_eff, 10, min(Xc.shape[1], Xn.shape[1]))

    recs = []
    for start in range(0, len(common), CHUNK):
        sl = slice(start, min(start + CHUNK, len(common)))
        A = np.abs(Xc[ic[sl], :]).astype(np.float64)
        B = np.abs(Xn[inn[sl], :]).astype(np.float64)
        Asig = Xc[ic[sl], :].astype(np.float64)
        Bsig = Xn[inn[sl], :].astype(np.float64)
        for j in range(A.shape[0]):
            gi = start + j
            k = int(n_eff[gi])
            a = thin_to(A[j][None, :], k)[0]
            b = thin_to(B[j][None, :], k)[0]
            asg = thin_to(Asig[j][None, :], k)[0]
            bsg = thin_to(Bsig[j][None, :], k)[0]

            if np.ptp(a) == 0 and np.ptp(b) == 0 and a[0] == b[0]:
                p, delta = 1.0, 0.0
            else:
                u, p = mannwhitneyu(a, b, alternative='two-sided')
                delta = 2.0 * u / (a.size * b.size) - 1.0
            if np.ptp(asg) == 0 and np.ptp(bsg) == 0 and asg[0] == bsg[0]:
                p_s, delta_s = 1.0, 0.0
            else:
                us, p_s = mannwhitneyu(asg, bsg, alternative='two-sided')
                delta_s = 2.0 * us / (asg.size * bsg.size) - 1.0

            recs.append({
                'context': ctx, 'rxn': common[gi],
                'n_eff': k,
                'ess_cancer': ess_c[gi], 'ess_normal': ess_n[gi],
                'rhat_cancer': rhat_c[ic[gi]], 'rhat_normal': rhat_n[inn[gi]],
                'median_abs_cancer': float(np.median(a)), 'median_abs_normal': float(np.median(b)),
                'mean_abs_cancer': float(a.mean()), 'mean_abs_normal': float(b.mean()),
                'median_signed_cancer': float(np.median(asg)),
                'median_signed_normal': float(np.median(bsg)),
                'p_abs': float(p), 'cliffs_delta_abs': float(delta),
                'p_signed': float(p_s), 'cliffs_delta_signed': float(delta_s),
            })
    df = pd.DataFrame(recs)
    df['fdr_abs'] = bh_fdr(df.p_abs.values)
    df['fdr_signed'] = bh_fdr(df.p_signed.values)
    df['call_up'] = (df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs >= DELTA_CUT)
    df['call_down'] = (df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs <= -DELTA_CUT)
    return df


def main():
    os.makedirs(OUT, exist_ok=True)
    frames, summary = [], []
    for ctx, c, n in zip(CONTEXTS, CANCER, NORMAL):
        print('=== %s (%s vs %s) ===' % (ctx, c, n), flush=True)
        df = compare_pair(ctx, c, n)
        frames.append(df)
        row = {'context': ctx, 'cancer': c, 'normal': n, 'n_common': len(df),
               'n_up': int(df.call_up.sum()), 'n_down': int(df.call_down.sum()),
               'n_fdr_sig': int((df.fdr_abs < FDR_CUT).sum()),
               'max_rhat': float(max(df.rhat_cancer.max(), df.rhat_normal.max())),
               'frac_rhat_gt_1p05': float(((df.rhat_cancer > 1.05) | (df.rhat_normal > 1.05)).mean()),
               'min_ess': float(min(df.ess_cancer.min(), df.ess_normal.min())),
               'median_ess': float(np.median(np.r_[df.ess_cancer, df.ess_normal]))}
        for g in DELTA_GRID:
            row['up_d%.3f' % g] = int(((df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs >= g)).sum())
            row['down_d%.3f' % g] = int(((df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs <= -g)).sum())
        summary.append(row)
        print('   common %d | FDR<0.05 %d | up %d | down %d | max Rhat %.4f'
              % (len(df), row['n_fdr_sig'], row['n_up'], row['n_down'], row['max_rhat']), flush=True)

    allf = pd.concat(frames, ignore_index=True)
    allf.to_csv(os.path.join(OUT, 'sampling_per_reaction.csv'), index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(OUT, 'sampling_summary.csv'), index=False)

    print('\n=== convergence ===')
    print(sm[['context', 'max_rhat', 'frac_rhat_gt_1p05', 'min_ess', 'median_ess']].to_string(index=False))
    print('\n=== calls ===')
    print(sm[['context', 'n_common', 'n_fdr_sig', 'n_up', 'n_down']].to_string(index=False))
    print('\n=== sensitivity to the Cliff\'s delta cutoff ===')
    cols = ['context'] + ['%s_d%.3f' % (d, g) for g in DELTA_GRID for d in ('up', 'down')]
    print(sm[cols].to_string(index=False))


if __name__ == '__main__':
    main()
