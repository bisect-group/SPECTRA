"""Figure: SPECTRA_CC vs FASTERCC vs FastCC vs SwiftCC on Recon3D.

Regenerates the blocked-reaction panel of the SPECTRA figure with FASTERCC
added as a fourth arm. The plotting function is the one from
Runtime_SPECTRA_CC.ipynb, unchanged except that it now takes four bars per
group instead of three; colours, log scale, fonts, axis labels and legend
placement are identical so the panel still matches the rest of the figure.

WHAT CHANGED IN THE DATA, and why the old numbers are not reused.
The published panel came from Results_SPECTRA_CC_Runtime.mat, a single sweep
of n = 20:20:10600 with one run per size, from which the notebook sampled ten
points with linspace(0, 528, 10). This reads instead the benchmark in
RuntimeComparisons/FASTERCC_comparison: the same construction - sample n
reactions to KEEP from the consistent Recon3D core - on a round
1,000-10,000 axis, with 3 independent draws x 3 timing repeats per size per
arm, on one machine with one solver.

THE TOLERANCE IS THE SOLVER DEFAULT (1e-6), matching what the published
panel was produced under. The rest of this benchmark pins feasTol = 1e-9,
which matters where the question is whether an answer is reproducible at
scale, but costs the LP-based arms several fold and would make this panel
incomparable with the published one. The 1e-9 rows are still in the timings
file under a different tag; RUN_TAG selects which family is plotted.

All four arms therefore have to come from this run. Quoting FASTERCC from here
beside FastCC and SwiftCC numbers produced on different hardware, a different
solver and a different tolerance would not be a comparison.

AGGREGATION. A bar is the plain MEAN of every measurement in its cell: 3
independent draws x 3 timing repeats = 9 runs, read from
exp1_papersizes_timings.csv, which carries the individual timings. The
per-run summary CSVs store only runtimeSec = median(times) per draw, so a
figure built from those could only ever be a mean of medians; this averages
the raw numbers instead.

NO ERROR BARS, deliberately. A symmetric mean +/- std renders asymmetrically
on a log axis and the lower whisker runs toward zero, which cannot be drawn
and reads as a far larger relative uncertainty than there is; with runtimes
spanning four orders of magnitude the bars are misleading rather than
informative. The spread across draws is printed in the table below instead.
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(
    HERE, "..", "RuntimeComparisons", "FASTERCC_comparison", "results_csv"
)

# Display names, in the order the bars should appear within each group.
ARM_LABEL = {
    "fastcc": "FastCC",
    "swiftcc": "SwiftCC",
    "fastercc_r0": "FASTERCC",
    "spectracc": "SPECTRA_CC",
}
ORDER = ["FastCC", "SwiftCC", "FASTERCC", "SPECTRA_CC"]

# Which family of runs to plot. "dtol" is the solver-default tolerance pass on
# the 1,000-10,000 axis, which is what the paper reports; drop this to
# "rep3_x3_n" to plot the earlier feasTol = 1e-9 pass on the published sizes.
RUN_TAG = "dtol"

# The three original colours are untouched. FASTERCC takes a brighter teal
# from the same family, so the panel reads as one palette.
#
# It is the nearest candidate to FastCC's #008080 - same hue, lighter - which
# is tolerable here only because the two are never adjacent: the bar order is
# FastCC, SwiftCC, FASTERCC, SPECTRA_CC, so SwiftCC's dark navy always sits
# between them. Keep that order if the arms are ever re-sorted.
COLORS = {
    "FastCC": "#008080",
    "SwiftCC": "#004080",
    "FASTERCC": "#00a0a0",
    "SPECTRA_CC": "#80c080",
}


def load():
    """Ten paper sizes x four arms, mean over all 9 runs per cell."""
    timings = os.path.join(CSV_DIR, "exp1_papersizes_timings.csv")
    raw = pd.read_csv(timings)
    raw = raw[raw["tag"].str.contains(RUN_TAG, regex=False)]
    if raw.empty:
        raise ValueError(f"no timing rows tagged {RUN_TAG!r} in {timings}")

    # An arm that failed on an instance contributes no timing rows at all, so
    # it simply has no bar at that size - the truthful rendering of "it did
    # not run". The per-run summaries carry the reason.
    summ = pd.concat(
        [pd.read_csv(f) for f in
         sorted(glob.glob(os.path.join(CSV_DIR, f"*{RUN_TAG}*_a4.csv")))],
        ignore_index=True)
    failed = summ[summ["error"].notna()]

    raw["Algorithm"] = raw["arm"].map(ARM_LABEL)
    g = (
        raw.groupby(["nRxns", "Algorithm"])
        .agg(
            runtime_mean=("runtimeSec", "mean"),
            runtime_std=("runtimeSec", "std"),
            LPs_mean=("nLPs", "mean"),
            LPs_std=("nLPs", "std"),
            nCons=("nConsistent", "mean"),
            runs=("runtimeSec", "size"),
        )
        .reset_index()
        .rename(columns={"nRxns": "n_rxns"})
    )
    return g, failed


def plot_bar_plot(df, variable, algos, ax):
    """Unchanged from Runtime_SPECTRA_CC.ipynb apart from the bar width."""
    labels = list(df["n_rxns"].unique())
    labels.sort()  # Ensure labels are sorted for consistent bar placement

    x = np.arange(len(labels))  # label locations
    width = 0.8 / len(algos)  # was 0.25 for three bars; keeps groups the same width
    # setting only the x-axis grid
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, linestyle="-", alpha=0.5)
    ax.set_axisbelow(True)

    # Bars with design-matching colors. Reindexing on the full label list puts
    # a NaN - and so no bar - wherever an arm has no result for that size,
    # instead of silently shifting the remaining bars left.
    for k, algo in enumerate(algos):
        sub = df[df["Algorithm"] == algo].set_index("n_rxns").reindex(labels)
        mean_values = sub[variable + "_mean"].values
        ax.bar(x + k * width, mean_values, width, color=COLORS[algo], label=algo)

    # Set x-axis labels
    ax.set_xticks(x + width * (len(algos) - 1) / 2)
    ax.set_xticklabels(labels, fontsize=14 * 2, rotation=45, ha="right")

    # setting the y-axis to log scale
    ax.set_yscale("log")
    ax.tick_params(axis="y", labelsize=14 * 2)
    ax.legend(
        fontsize=14 * 2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=len(algos),
        frameon=False,
    )

    ylabel_dict = {"runtime": "Runtime in seconds", "LPs": "Number of LPs solved"}

    # Remove chart border (top/right)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.set_xlabel("Number of reactions in the input model", fontsize=16 * 2, fontweight="bold")
    ax.set_ylabel(ylabel_dict[variable], fontsize=16 * 2, fontweight="bold")


def main():
    df, failed = load()

    fig, ax = plt.subplots(1, 2, figsize=(30, 10))
    plot_bar_plot(df, "runtime", ORDER, ax[0])
    # SwiftCC is excluded from the LP panel, as in the published figure: it
    # settles reversible reactions by QR factorisation, so its "1 LP" is not
    # comparable. FASTERCC is included - it is an LP method throughout.
    plot_bar_plot(df, "LPs", ["FastCC", "FASTERCC", "SPECTRA_CC"], ax[1])
    plt.tight_layout()
    out = os.path.join(HERE, "SPECTRA_CC_FASTERCC_runtime.png")
    plt.savefig(out, bbox_inches="tight", dpi=150)
    print(f"wrote {out}")

    # The numbers behind the bars, for the manuscript text.
    print("\n=== mean over all runs per cell (3 draws x 3 timings) ===")
    piv = df.pivot(index="n_rxns", columns="Algorithm", values="runtime_mean")[ORDER]
    print(piv.to_string(float_format=lambda v: f"{v:9.3f}"))
    print("\n=== LPs ===")
    print(
        df.pivot(index="n_rxns", columns="Algorithm", values="LPs_mean")[
            ["FastCC", "FASTERCC", "SPECTRA_CC"]
        ].to_string(float_format=lambda v: f"{v:7.0f}")
    )
    print("\n=== totals over the ten sizes (s) ===")
    tot = piv.sum()
    for a in ORDER:
        print(f"  {a:<12} {tot[a]:10.1f}   {tot[a]/tot['SPECTRA_CC']:5.1f}x SPECTRA_CC")
    print("\n=== fastest arm per size ===")
    for n, row in piv.iterrows():
        print(f"  n={n:<6} {row.idxmin():<12} {row.min():8.3f} s")
    if len(failed):
        print("\n=== arms that did not complete (no bar drawn) ===")
        for _, r in failed.iterrows():
            print(f"  n={r.nRxns:<6} {r.arm:<12} {str(r.error).strip()[:70]}")


if __name__ == "__main__":
    main()
