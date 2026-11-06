"""Step 19: draw the main figure: sensitivity and specificity at each step (G1 to G4), with 95% confidence intervals.
Uses the 622 isolates where all four gene groups gave a definite answer, so every point is on the same isolates.
Reads:  results/results_summary.csv
Writes: figures/sensitivity_specificity_by_step.png
Usage:  python scripts/19_make_figure.py"""
import os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "figures"); os.makedirs(OUT, exist_ok=True)
SET = "common set, evaluable for all of R1-R4"

# colours: two series, validated for colour-blind separation and contrast
SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
BLUE, ORANGE = "#2a78d6", "#eb6834"

r = pd.read_csv(os.path.join(ROOT, "results", "results_summary.csv"))
r = r[r.analysis_set == SET].set_index("rule").loc[["G1", "G2", "G3", "G4"]]
def parse(s):
    """'81.1% (74.3-86.8)' -> (81.1, 74.3, 86.8)"""
    v, lo, hi = map(float, re.match(r"([\d.]+)% \(([\d.]+)-([\d.]+)\)", s).groups())
    return v, lo, hi
sens = [parse(s) for s in r["sensitivity_pct_(95%_cluster_bootstrap_CI)"]]
spec = [parse(s) for s in r["specificity_pct_(95%_cluster_bootstrap_CI)"]]
n, n_r, n_s = int(r.evaluable_n.iloc[0]), int(r.resistant.iloc[0]), int(r.susceptible.iloc[0])

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
                     "axes.labelcolor": INK_2, "xtick.color": INK_2, "ytick.color": INK_2})
fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
x = [0, 1, 2, 3]
for vals, colour, name, off in ((sens, BLUE, "Sensitivity", -0.06), (spec, ORANGE, "Specificity", 0.06)):
    xs = [i + off for i in x]   # small sideways offset so the two sets of confidence bars do not overlap
    v = [a for a, _, _ in vals]; lo = [a - b for a, b, _ in vals]; hi = [c - a for a, _, c in vals]
    ax.errorbar(xs, v, yerr=[lo, hi], color=colour, lw=2, capsize=0, elinewidth=1.2, alpha=1,
                marker="o", ms=8, mfc=colour, mec=SURFACE, mew=2, solid_capstyle="round", solid_joinstyle="round",
                label=name, zorder=3)
    # direct labels at the two ends only
    for i, ha, dx in ((0, "right", -0.2), (3, "left", 0.2)):
        ax.text(x[i] + dx, v[i], f"{v[i]:.1f}%", ha=ha, va="center", fontsize=10, color=INK, zorder=4)

# where the specificity drop comes from
ax.annotate("$\\it{ampC}$ promoter mutations\n(mostly T-32A)", xy=(3.06, spec[3][0] - 1.0), xytext=(2.35, 78.5),
            fontsize=9, color=INK_2, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8, shrinkA=2, shrinkB=6))

ax.set_xticks(x, ["G1\nCTX-M", "G2\n+ other ESBL", "G3\n+ acquired AmpC", "G4\n+ $\\it{ampC}$ promoter"])
ax.set_xlim(-0.6, 3.6); ax.set_ylim(70, 100)
ax.set_yticks(range(70, 101, 5), [f"{t}%" for t in range(70, 101, 5)])
ax.yaxis.grid(True, color=GRID, lw=1); ax.set_axisbelow(True)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color(GRID); ax.tick_params(length=0, pad=6)
fig.legend(loc="upper left", bbox_to_anchor=(0.055, 0.885), ncol=2, frameon=False, fontsize=10, handlelength=1.6, columnspacing=2.0)
fig.text(0.06, 0.95, "Adding genes: sensitivity rises, specificity falls at the last step", fontsize=13, color=INK, weight="bold")
fig.text(0.06, 0.905, f"Ceftazidime, {n} $\\it{{E. coli}}$ isolates ({n_r} resistant, {n_s} susceptible). Bars show 95% confidence intervals.",
         fontsize=9.5, color=INK_2)
fig.subplots_adjust(left=0.08, right=0.95, top=0.80, bottom=0.14)
path = os.path.join(OUT, "sensitivity_specificity_by_step.png")
fig.savefig(path, facecolor=SURFACE)
print("saved", path)
