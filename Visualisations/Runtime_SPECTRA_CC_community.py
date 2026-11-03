"""Supplementary figure: SPECTRA_CC vs FASTERCC on AGORA2 community models.

Communities of 50 and 100 AGORA2 reconstructions joined through a shared [u]
compartment, under a random medium (each community uptake opened with
probability 0.30, open uptakes drawn log-uniformly with a 1e-4 floor).
Reactions become inconsistent through the medium rather than through
deletion, which is how blocked content arises in practice and blocks in
correlated clusters rather than at random.

CONFIGURATION, and why. Media floor 1e-4 with epsilon 1e-4 and feasTol 1e-9.
At epsilon 1e-8 the consistent sets move by at most 0.2%, but FASTERCC fails
on five of the six instances with an infeasible LP7, leaving nothing to
compare; epsilon 1e-4 is the setting at which both methods return an answer
on every instance. It is also the conservative choice for SPECTRA_CC, whose
0.99*epsilon detection cutoff then sits level with the media floor and costs
it LPs it does not need at 1e-8.

The x axis is the reaction count of each joined community rather than the
species count, so the six instances order by size: three draws at 50 species
span 79,948-101,965 reactions and three at 100 species span
169,213-180,874.

Bars are the mean over the timing repeats, as in the other panels, and carry
no error bars for the same reason.
"""

import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(
    HERE, "..", "RuntimeComparisons", "FASTERCC_comparison", "results_csv"
)

ARM_LABEL = {"fastercc_r0": "FASTERCC", "spectracc": "SPECTRA_CC"}
ORDER = ["FASTERCC", "SPECTRA_CC"]
COLORS = {"FASTERCC": "#00a0a0", "SPECTRA_CC": "#80c080"}


def load():
    """The floor-1e-4, epsilon-1e-4 community runs, one row per instance/arm."""
    files = sorted(glob.glob(os.path.join(CSV_DIR, "*m4_e4_n*_a2.csv")))
    # The single-arm reruns: one FASTERCC cell lost its first two attempts to
    # a licence refusal and a wall-clock timeout, and was re-run alone.
    files += sorted(glob.glob(os.path.join(CSV_DIR, "*m4_e4_n*fccONLY*.csv")))
    if not files:
        raise FileNotFoundError(f"no community CSVs under {CSV_DIR}")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)

    d = d[d["error"].isna() | (d["error"].astype(str).str.strip() == "")]
    d["Algorithm"] = d["arm"].map(ARM_LABEL)
    d = d[d["Algorithm"].notna()]

    g = (d.groupby(["nRxns", "Algorithm"])
           .agg(runtime=("runtimeSec", "mean"), LPs=("nLPs", "mean"),
                cons=("nConsistent", "mean"))
           .reset_index())
    return g, d


def plot_bar_plot(df, variable, algos, ax):
    labels = sorted(df["n_rxns"].unique())
    x = np.arange(len(labels))
    width = 0.8 / len(algos)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, linestyle="-", alpha=0.5)
    ax.set_axisbelow(True)

    for k, algo in enumerate(algos):
        sub = df[df["Algorithm"] == algo].set_index("n_rxns").reindex(labels)
        ax.bar(x + k * width, sub[variable].values, width,
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
    ax.set_xlabel("Number of reactions in the community model",
                  fontsize=16 * 2, fontweight="bold")
    ax.set_ylabel(ylabel[variable], fontsize=16 * 2, fontweight="bold")


def main():
    g, raw = load()
    g = g.rename(columns={"nRxns": "n_rxns"})

    fig, ax = plt.subplots(1, 2, figsize=(30, 10))
    plot_bar_plot(g, "runtime", ORDER, ax[0])
    plot_bar_plot(g, "LPs", ORDER, ax[1])
    plt.tight_layout()
    out = os.path.join(HERE, "SPECTRA_CC_community_runtime.png")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    print(f"wrote {out}")

    print("\n=== communities, floor 1e-4, epsilon 1e-4 ===")
    piv = g.pivot(index="n_rxns", columns="Algorithm")
    print(piv.to_string(float_format=lambda v: f"{v:10.1f}"))
    sp = g[g.Algorithm == "SPECTRA_CC"].set_index("n_rxns")
    fc = g[g.Algorithm == "FASTERCC"].set_index("n_rxns")
    print("\n=== speedup and LP ratio, per instance ===")
    for n in sorted(sp.index):
        if n not in fc.index:
            continue
        print(f"  {n:>7,} rxns   {fc.runtime[n]/sp.runtime[n]:5.1f}x faster   "
              f"{fc.LPs[n]/sp.LPs[n]:5.1f}x fewer LPs   "
              f"SPECTRA +{sp.cons[n]-fc.cons[n]:.0f} consistent")


if __name__ == "__main__":
    main()
