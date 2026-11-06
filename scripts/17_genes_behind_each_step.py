"""Step 17: which genes make the difference at each step (G1 -> G2 -> G3 -> G4), split by ceftazidime category.
Also prints the MICs of the false negatives and false positives.
Reads:  results/linked_isolates_749.csv, data/gene_hits_749.csv
Writes: results/step_additions_by_determinant.csv
Usage:  python scripts/17_genes_behind_each_step.py"""
import os
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); T = os.path.join(ROOT, "results")
L = pd.read_csv(os.path.join(T, "linked_isolates_749.csv"), dtype={"genome_id": str}).set_index("genome_id")
Hs = pd.read_csv(os.path.join(ROOT, "data", "gene_hits_749.csv"), dtype={"genome": str})
pos = Hs[Hs.state == "POSITIVE"]
alle = pos.groupby(["genome", "component"]).sym.apply(lambda s: "+".join(sorted(set(s)))).unstack()
rows = []
def tab(name, genomes, label):
    sub = L.loc[genomes]; vc = sub.caz_category.value_counts()
    rows.append(dict(group=name, determinant=label, genomes=len(sub), R=int(vc.get("R", 0)), I=int(vc.get("I", 0)), S=int(vc.get("S", 0)), missing_or_no_record=int(vc.get("missing", 0) + vc.get("no_record", 0))))
# G1 positives by allele
g1 = alle["G1"].dropna()
for a in sorted(set(g1)): tab("G1 positives by CTX-M allele", list(g1[g1 == a].index), a)
# newly positive at each step (not positive in the previous rule)
for prev, cur, comp in (("R1", "R2", "G2"), ("R2", "R3", "G3"), ("R3", "R4", "G4")):
    new = L.index[(L[prev] != "POSITIVE") & (L[cur] == "POSITIVE")]
    a = alle.reindex(new)[comp]
    for v in sorted(a.dropna().unique()): tab(f"newly positive at {cur} (component {comp}), not positive at {prev}", list(a[a == v].index), v)
    other = [g for g in new if pd.isna(alle.loc[g, comp]) if g in alle.index]
    if len(other): tab(f"newly positive at {cur}", other, "(not explained by component %s)" % comp)
    tab(f"newly positive at {cur} (all)", list(new), "ALL")
R = pd.DataFrame(rows); R.to_csv(os.path.join(T, "step_additions_by_determinant.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80); print(R.to_string(index=False))
# false negatives of the broadest primary rule (R4): resistant, determinate NEGATIVE
fn = L[(L.caz_category == "R") & (L.R4 == "NEGATIVE")]; print("\nR4 false negatives (resistant, genotype NEGATIVE):", len(fn), "| MIC values:", sorted(fn.caz_mic_numeric.tolist()))
fn1 = L[(L.caz_category == "R") & (L.R1 == "NEGATIVE")]; print("R1 false negatives:", len(fn1), "| MIC:", sorted(fn1.caz_mic_numeric.tolist()))
tp = L[(L.caz_category == "R") & (L.R1 == "POSITIVE")]; print("R1 true positives MIC distribution:", tp.caz_mic_numeric.value_counts().sort_index().to_dict())
fp = L[(L.caz_category == "S") & (L.R1 == "POSITIVE")]; print("R1 false positives (S but CTX-M positive) MIC distribution:", fp.caz_mic_numeric.value_counts().sort_index().to_dict())
