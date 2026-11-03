"""Figure 2a: SPECTRA_ME vs Fastcore vs Swiftcore, averaged over three draws.

The published panel came from Results_SPECTRA_ME_Runtime.mat, a single random
core set per size. ForSPECTRA_MEResults.m already drew its core sets at
random - what it did not do was draw more than once, so each bar rested on
one sample. This reads three independent draws instead and plots their MEAN,
matching how the CC and CC_ME panels are now built.

Getting three draws needed two fixes to that script, both recorded there: a
plain `clear` wiped the caller's loop variables, and MATLAB starts every
session from the same RNG state, so three separate runs would have drawn the
IDENTICAL core sets and reported a spread of exactly zero.

NO ERROR BARS, matching the other two panels: a symmetric mean +/- std is
misleading on these axes, and the published figure carries none.

The plotting function, the three arms, the colours, the bar width, the font
sizes and the linear y-axis are all as in Runtime_SPECTRA_ME.ipynb. Only the
source of the numbers changed.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
REP_DIR = os.path.join(HERE, "Results", "ME_reps")
REPS = [1, 2, 3]

# MATLAB variable -> display name, as the published notebook labels them.
ARMS = {
    "fc": "Fastcore",
    "swWR": "Swiftcore-WR",
    "sp": "SPECTRA_ME",
}
ORDER = ["Fastcore", "Swiftcore-WR", "SPECTRA_ME"]


def load():
    """Mean over the three draws, per size per arm."""
    frames = []
    for rep in REPS:
        p = os.path.join(REP_DIR, f"Results_SPECTRA_ME_Runtime_rep{rep}.mat")
        s = sio.loadmat(p)
        n = np.asarray(s["n"]).ravel()
        for key, label in ARMS.items():
            frames.append(
                pd.DataFrame({
                    "n_rxns": n,
                    "Algorithm": label,
                    "runtime": np.asarray(s[f"{key}times"]).ravel(),
                    "LPs": np.asarray(s[f"{key}LPS"]).ravel(),
                    "rep": rep,
                })
            )
    raw = pd.concat(frames, ignore_index=True)
    g = (raw.groupby(["n_rxns", "Algorithm"])
            .agg(runtime=("runtime", "mean"), LPs=("LPs", "mean"),
                 draws=("runtime", "size"))
            .reset_index())
    return g, raw


def plot_bar_plot(df, variable, algos, ax):
    """As in Runtime_SPECTRA_ME.ipynb, unchanged."""
    labels = list(df["n_rxns"].unique())
    labels.sort()  # Ensure labels are sorted for consistent bar placement
    values = [df[df["Algorithm"] == algo][variable].values for algo in algos]

    x = np.arange(len(labels))  # label locations
    width = 0.25  # width of the bars
    # setting only the x-axis grid
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, linestyle="-", alpha=0.5)
    ax.set_axisbelow(True)
    colors = {"Fastcore": "#008080", "Swiftcore-WR": "#004080",
              "Swiftcore-R": "#004080", "SPECTRA_ME": "#80c080",
              "SPECTRA_ME2": "#80c080"}
    # Bars with design-matching colors
    for k in range(len(algos)):
        ax.bar(x + k * width, values[k], width, color=colors[algos[k]], label=algos[k])

    # Set x-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=14 * 2, rotation=45, ha="right")

    ax.tick_params(axis="y", labelsize=14 * 2)
    ax.legend(fontsize=14 * 2, loc="upper center", bbox_to_anchor=(0.5, -0.25),
              ncol=3, frameon=False)

    ylabel_dict = {"runtime": "Runtime in seconds", "LPs": "Number of LPs solved"}

    # Remove chart border (top/right)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.set_xlabel("Number of core reactions", fontsize=16 * 2, fontweight="bold")
    ax.set_ylabel(ylabel_dict[variable], fontsize=16 * 2, fontweight="bold")


def main():
    df, raw = load()

    fig, ax = plt.subplots(1, 2, figsize=(30, 10))
    plot_bar_plot(df, "runtime", ORDER, ax[0])
    plot_bar_plot(df, "LPs", ORDER, ax[1])
    plt.tight_layout()
    out = os.path.join(HERE, "SPECTRA_ME_runtime_reps.png")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    print(f"wrote {out}")

    print("\n=== runtime, mean of 3 draws (s) ===")
    print(df.pivot(index="n_rxns", columns="Algorithm", values="runtime")[ORDER]
            .to_string(float_format=lambda v: f"{v:8.2f}"))
    print("\n=== LPs, mean of 3 draws ===")
    print(df.pivot(index="n_rxns", columns="Algorithm", values="LPs")[ORDER]
            .to_string(float_format=lambda v: f"{v:7.1f}"))
    print("\n=== means over all sizes ===")
    for a in ORDER:
        s = df[df.Algorithm == a]
        print(f"  {a:<14} {s.runtime.mean():8.2f} s   {s.LPs.mean():7.1f} LPs")
    print("\n=== draw-to-draw spread (max-min of the 3 draws), for the text ===")
    sp = (raw.groupby(["n_rxns", "Algorithm"]).runtime
             .agg(lambda v: v.max() - v.min()).reset_index())
    print(sp.pivot(index="n_rxns", columns="Algorithm", values="runtime")[ORDER]
            .to_string(float_format=lambda v: f"{v:8.2f}"))


if __name__ == "__main__":
    main()
