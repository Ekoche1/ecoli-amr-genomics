"""Step 18: the headline results table, plus best- and worst-case bounds for the INDETERMINATE calls.
Also recounts every 2x2 table directly from the linked table, as an independent check on step 16.
Reads:  results/linked_isolates_749.csv, results/results_all_sets_long.csv
Writes: results/results_summary.csv, results/indeterminate_bounds.csv
Usage:  python scripts/18_summary_tables.py"""
import os
import pandas as pd
T = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
L = pd.read_csv(os.path.join(T, "linked_isolates_749.csv"), dtype={"genome_id": str}).set_index("genome_id")
R = pd.read_csv(os.path.join(T, "results_all_sets_long.csv"))
# independent recount from the saved linked table (does not use the analysis script's code path)
ana = L[L.caz_category.isin(["R", "S"])]
def counts(rule, mask):
    s = ana[mask(ana)]; pos = s[rule] == "POSITIVE"; r = s.caz_category == "R"
    return int((pos & r).sum()), int((~pos & r).sum()), int((pos & ~r).sum()), int((~pos & ~r).sum())
det = lambda k: (lambda d: d[k] != "INDETERMINATE")
sets = {"own evaluable set": lambda k: det(k), "common set C(R1,R4)": lambda k: (lambda d: (d.R1 != "INDETERMINATE") & (d.R4 != "INDETERMINATE")),
        "common set, evaluable for all of R1-R4": lambda k: (lambda d: (d.R1 != "INDETERMINATE") & (d.R2 != "INDETERMINATE") & (d.R3 != "INDETERMINATE") & (d.R4 != "INDETERMINATE")),
        "common set C(R4,R5)": lambda k: (lambda d: (d.R4 != "INDETERMINATE") & (d.R5 != "INDETERMINATE"))}
bad = 0
for _, r in R.iterrows():
    if counts(r.rule, sets[r["set"]](r.rule)) != (r.TP, r.FN, r.FP, r.TN): bad += 1; print("MISMATCH", r["set"], r.rule)
print("independent recount of", len(R), "2x2 tables from linked CSV: mismatches =", bad)
# concise table
f = lambda v, lo, hi: f"{100*v:.1f}% ({100*lo:.1f}-{100*hi:.1f})"
rows = []
for _, r in R.iterrows():
    rows.append({"analysis_set": r["set"], "rule": {"R1": "G1", "R2": "G2", "R3": "G3", "R4": "G4", "R5": "G5 (supplementary)"}[r.rule], "evaluable_n": r.evaluable_n, "resistant": r.evaluable_R, "susceptible": r.evaluable_S, "TP": r.TP, "FN": r.FN, "FP": r.FP, "TN": r.TN,
                 "sensitivity_pct_(95%_cluster_bootstrap_CI)": f(r.sens, r.sens_boot_lo, r.sens_boot_hi), "specificity_pct_(95%_cluster_bootstrap_CI)": f(r.spec, r.spec_boot_lo, r.spec_boot_hi),
                 "sensitivity_Wilson_CI": f"{100*r.sens_wilson_lo:.1f}-{100*r.sens_wilson_hi:.1f}", "specificity_Wilson_CI": f"{100*r.spec_wilson_lo:.1f}-{100*r.spec_wilson_hi:.1f}",
                 "PPV_pct_(sample-specific)": f"{100*r.ppv:.1f}", "NPV_pct_(sample-specific)": f"{100*r.npv:.1f}", "indeterminate_set_aside_(with_R/S_phenotype)": int(r.indeterminate_with_phenotype), "of_which_R": int(r.indet_R), "of_which_S": int(r.indet_S)})
pd.DataFrame(rows).to_csv(os.path.join(T, "results_summary.csv"), index=False)
# arithmetic bounds for indeterminate calls on the full analysable set (644: 207 R, 437 S); S5a/S5b-style bounds, no model
b = []
for k in ["R1", "R2", "R3", "R4", "R5"]:
    s = ana[k]; r = ana.caz_category == "R"; TP = int(((s == "POSITIVE") & r).sum()); FP = int(((s == "POSITIVE") & ~r).sum()); TN = int(((s == "NEGATIVE") & ~r).sum()); iR = int(((s == "INDETERMINATE") & r).sum()); iS = int(((s == "INDETERMINATE") & ~r).sum()); nR, nS = int(r.sum()), int((~r).sum())
    b.append(dict(rule=k, analysable_R=nR, analysable_S=nS, indet_R=iR, indet_S=iS, sens_if_indet_missed=TP / nR, sens_determinate_only=TP / (nR - iR), sens_if_indet_detected=(TP + iR) / nR, spec_if_indet_flagged=TN / nS, spec_determinate_only=TN / (nS - iS), spec_if_indet_negative=(TN + iS) / nS))
B = pd.DataFrame(b).round(4); B.to_csv(os.path.join(T, "indeterminate_bounds.csv"), index=False); print(B.to_string(index=False))
