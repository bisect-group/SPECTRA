"""Supplementary figure: SPECTRA_CC on expanded Recon3D reconstructions.

Expanded models are built by the scFASTCORMICS expansion, whose reaction
count is affine in the cluster count k: n(k) = 9040k + 3120. Five sizes from
47,820 to 454,620 reactions, each with four removal levels (500, 2,000,
5,000, 10,000 reactions deleted at random), three draws and three timings.

WHAT THE FIGURE SHOWS. SPECTRA_CC solves 10 LPs at every size from 93,020
reactions upward - the count does not move across a 4.9-fold change in model
size, nor across a 20-fold change in how many reactions were deleted.
FASTERCC has no bar beyond 93,020 reactions because it returns no result
there.

WHY FASTERCC STOPS. Not the consistency algorithm: an instrumented run
showed it completing its entire main loop on the 183,420-reaction model at
3.23 GB. It fails afterwards, in restoreBounds, on

    SameDirr = find(all(S_orig == S_cons, 1))

where `sparse == sparse` is true wherever BOTH operands are zero, so the
result carries one entry per matrix ELEMENT rather than per nonzero. At
k=50 MATLAB names the array it refuses to allocate - 216,093 x 453,736,
821.8 GB - and those dimensions are exactly (metabolites x reactions) of the
consistent model. The cost therefore scales with the size of the ANSWER, not
of the input, which is why it bites hardest here: these models are ~99%
consistent.

Bars are means over all runs in each cell, as in the other panels, and carry
no error bars for the same reason.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "RuntimeComparisons", "FASTERCC_comparison",
                   "results_csv", "exp2_expanded_all.csv")

ARM_LABEL = {"fastercc_r0": "FASTERCC", "spectracc": "SPECTRA_CC"}
ORDER = ["FASTERCC", "SPECTRA_CC"]
COLORS = {"FASTERCC": "#00a0a0", "SPECTRA_CC": "#80c080"}


def load():
    d = pd.read_csv(CSV)
    d["Algorithm"] = d["arm"].map(ARM_LABEL)
    ok = d[d["runtimeSec"].notna()]
    g = (ok.groupby(["k", "Algorithm"])
           .agg(n_rxns=("nRxns", "max"), runtime=("runtimeSec", "mean"),
                LPs=("nLPs", "mean"), runs=("runtimeSec", "size"))
           .reset_index())
    return g, d


def plot_bar_plot(df, variable, algos, ax):
    ks = sorted(df["k"].unique())
    labels = [int(df[df.k == k].n_rxns.max()) for k in ks]
    x = np.arange(len(ks))
    width = 0.8 / len(algos)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, linestyle="-", alpha=0.5)
    ax.set_axisbelow(True)

    for i, algo in enumerate(algos):
        sub = df[df["Algorithm"] == algo].set_index("k").reindex(ks)
        ax.bar(x + i * width, sub[variable].values, width,
               color=COLORS[algo], label=algo)

    ax.set_xticks(x + width * (len(algos) - 1) / 2)
    ax.set_xticklabels([f"{v:,}" for v in labels], fontsize=14 * 2,
                       rotation=45, ha="right")
    ax.set_yscale("log")
    ax.tick_params(axis="y", labelsize=14 * 2)
    ax.legend(fontsize=14 * 2, loc="upper center", bbox_to_anchor=(0.5, -0.25),
              ncol=len(algos), frameon=False)

    ylabel = {"runtime": "Runtime in seconds", "LPs": "Number of LPs solved"}
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Number of reactions in the input model",
                  fontsize=16 * 2, fontweight="bold")
    ax.set_ylabel(ylabel[variable], fontsize=16 * 2, fontweight="bold")


def main():
    g, d = load()

    fig, ax = plt.subplots(1, 2, figsize=(30, 10))
    plot_bar_plot(g, "runtime", ORDER, ax[0])
    plot_bar_plot(g, "LPs", ORDER, ax[1])
    plt.tight_layout()
    out = os.path.join(HERE, "SPECTRA_CC_clusters_runtime.png")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    print(f"wrote {out}")

    print("\n=== expanded Recon3D, mean over all runs per cell ===")
    print(f"{'k':>4} {'reactions':>11} {'arm':<12} {'runs':>5} {'runtime s':>11} {'LPs':>7}")
    for _, r in g.sort_values(["k", "Algorithm"]).iterrows():
        print(f"{int(r.k):>4} {int(r.n_rxns):>11,} {r.Algorithm:<12} "
              f"{int(r.runs):>5} {r.runtime:>11.1f} {r.LPs:>7.1f}")

    bad = d[d["runtimeSec"].isna()]
    print("\n=== runs with no result ===")
    for (k, arm), grp in bad.groupby(["k", "arm"]):
        print(f"  k={int(k):<3} {arm:<12} {len(grp):>3} of "
              f"{len(d[(d.k == k) & (d.arm == arm)]):>3}")


if __name__ == "__main__":
    main()
