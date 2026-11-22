"""fea.py -- Python re-implementation of the COBRA Toolbox FEA function.

Mirrors src/analysis/fluxEnrichmentAnalysis/FEA.m so that the sampling-based
reaction sets are enriched with exactly the same statistics the paper used:
an upper-tail hypergeometric test per subsystem, followed by Benjamini-Hochberg
FDR correction applied only to the subsystems actually hit by the query set.

Validated against the published fea_2025_checked_2.xlsx in s07_fea.py.
"""
import numpy as np
from scipy.stats import hypergeom


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values, matching MATLAB mafdr(...,'BHFDR',true)."""
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]   # enforce monotonicity
    adj = np.minimum(adj, 1.0)
    out = np.empty(n)
    out[order] = adj
    return out


def fea(subsystems, rxn_idx):
    """Flux enrichment analysis over `subsystems` for the reactions in `rxn_idx`.

    subsystems : sequence of str, one per reaction in the model (the background)
    rxn_idx    : integer indices into `subsystems` forming the query set

    Returns a list of dicts sorted by ascending p-value, one per subsystem that
    the query set hits.
    """
    subsystems = list(subsystems)
    total = len(subsystems)
    rxn_idx = np.asarray(list(rxn_idx), int)

    uniq = sorted(set(subsystems))                       # MATLAB unique() sorts
    idx_of = {s: i for i, s in enumerate(uniq)}

    n_per_group = np.zeros(len(uniq), int)               # K: successes in population
    for s in subsystems:
        n_per_group[idx_of[s]] += 1

    hits = np.zeros(len(uniq), int)                      # X: successes drawn
    for i in rxn_idx:
        hits[idx_of[subsystems[i]]] += 1

    draws = len(rxn_idx)                                 # N
    # hygecdf(X-1, M, K, N, 'upper') == P(x >= X)
    pvals = hypergeom.sf(hits - 1, total, n_per_group, draws)

    nonzero = np.flatnonzero(hits)
    order = nonzero[np.argsort(pvals[nonzero], kind='stable')]
    adj = bh_fdr(pvals[order])

    return [{'p': float(pvals[g]), 'p_adj': float(adj[k]), 'group': uniq[g],
             'enriched_set_size': int(hits[g]), 'total_set_size': int(n_per_group[g])}
            for k, g in enumerate(order)]
