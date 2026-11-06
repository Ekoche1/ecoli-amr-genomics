"""Step 14: apply the gene-calling rules (gene_rules.py) to the AMRFinderPlus hits of the 749 retained genomes.
Prints checks on the hits (methods, coverage, any beta-lactam hit the rules ignore) and summary counts.
Reads:  data/amrfinderplus_hits_749.csv, data/genome_qc_762.csv, data/allele_classes_2026-08-07.1.csv
Writes: data/gene_hits_749.csv (one row per relevant hit), data/gene_calls_749.csv (one row per genome),
        and summary tables in work/
Usage:  python scripts/14_call_genes.py"""
import os, sys, collections
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(os.path.dirname(HERE), "data"); REP = os.path.join(os.path.dirname(HERE), "work"); os.makedirs(REP, exist_ok=True)
sys.path.insert(0, HERE); import gene_rules as H
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 200); pd.set_option("display.max_colwidth", 60)
t26 = H.load_table("2026-08-07.1")
allh = pd.read_csv(os.path.join(DATA, "amrfinderplus_hits_749.csv"), dtype=str)
d = pd.read_csv(os.path.join(DATA, "genome_qc_762.csv"), dtype={"id": str}); retained = sorted(set(d[d.tier != "CONFIRMED"].id))
allh = allh[allh.genome.isin(retained)].copy(); allh["cov"] = pd.to_numeric(allh["% Coverage of reference"], errors="coerce"); allh["ident"] = pd.to_numeric(allh["% Identity to reference"], errors="coerce")
amr = allh[allh.Type == "AMR"]
print("=== A. STRUCTURE OF THE ACTUAL OUTPUT (749 genomes) ===")
print("AMR hit rows:", len(amr), "| Method values (all AMR):", dict(collections.Counter(amr.Method)))
print("Coverage missing (NaN) by method:", dict(amr[amr["cov"].isna()].Method.value_counts()))
key = amr[amr["Element symbol"].str.match(H.FAM) | amr["Element symbol"].str.startswith("ampC_")]
print("\nKey-family (+ampC_) hit rows:", len(key), "in", key.genome.nunique(), "genomes")
print(key.groupby("Method")["cov"].agg(["count", "min", "max"]).round(2))
comp = key[key.Method.isin(H.COMPLETE_METHODS)]; part = key[~key.Method.isin(H.COMPLETE_METHODS) & ~key["Element symbol"].str.startswith("ampC_")]
print("\nCOVERAGE GAP CHECK (non-point hits). complete-method coverage min/max:", comp["cov"].min(), comp["cov"].max(), "| incomplete-method coverage min/max:", part["cov"].min(), part["cov"].max())
print("complete-method hits with coverage <90:", int((comp["cov"] < 90).sum()), "| incomplete-method hits with coverage >=90:", int((part["cov"] >= 90).sum()))
print("hits with coverage strictly between 89.63 and 90.60 (any method):", int(((key["cov"] > 89.63) & (key["cov"] < 90.6)).sum()))
print("incomplete-method hits by method:", dict(part.Method.value_counts()))
print("ampC_ hit methods:", dict(key[key["Element symbol"].str.startswith("ampC_")].Method.value_counts()), "| their coverage:", key[key["Element symbol"].str.startswith("ampC_")]["cov"].describe()[["min", "max"]].to_dict())
# any beta-lactam-relevant AMR hit that the rules would silently ignore?
bl = amr[(amr.Class == "BETA-LACTAM") | amr.Subclass.fillna("").str.contains("CEPHALOSPORIN|CARBAPENEM|BETA-LACTAM")]
ign = bl[~(bl["Element symbol"].str.match(H.FAM) | bl["Element symbol"].str.startswith("ampC_"))]
print("\nBeta-lactam-class AMR hits NOT picked up by the family regex (ignored by rules):", len(ign)); print(ign.groupby(["Element symbol", "Subclass", "Method"]).size().sort_values(ascending=False).head(25))
# alleles absent from the 2026 table
syms = key["Element symbol"][~key["Element symbol"].str.startswith("ampC_")]
miss = syms[~syms.isin(t26.keys())]
print("\nKey-family symbols NOT in the 2026 allele table:", dict(collections.Counter(miss)) or "none")
print("Key symbols in table with class INDETERMINATE:", dict(collections.Counter(s for s in syms if t26.get(s) == "INDETERMINATE")))
# ---- per-genome states
hit_rows = []; comp_rows = []
for g in retained:
    df = amr[amr.genome == g].copy()
    for c in ("Element symbol", "Method", "% Coverage of reference", "Type"): pass
    hits, cs = H.hits_to_states(df, "new", t26)
    rs = H.rule_states(cs)
    comp_rows.append(dict(genome=g, **{f"C_{k}": v for k, v in cs.items()}, **rs))
    if len(hits): hits.insert(0, "genome", g); hit_rows.append(hits)
hits = pd.concat(hit_rows, ignore_index=True); G = pd.DataFrame(comp_rows)
print("\nretained genomes classified:", len(G), "| unique:", G.genome.nunique())
hits.to_csv(os.path.join(DATA, "gene_hits_749.csv"), index=False)
G.to_csv(os.path.join(DATA, "gene_calls_749.csv"), index=False)
# ---- state of every hit & reason for indeterminate
def reason(r):
    if r.state != "INDETERMINATE": return ""
    if r.sym.startswith("ampC_"): return "ampC promoter variant not in pre-specified set"
    if r.method in ("PARTIALX", "PARTIAL_CONTIG_ENDX", "PARTIAL"): return "partial hit (" + r.method + ")"
    if r.method == "INTERNAL_STOP": return "internal stop"
    if r.method in H.COMPLETE_METHODS and r.coverage < 90: return "complete-method hit but coverage <90%"
    if r.sym in t26 and t26[r.sym] == "INDETERMINATE": return "allele class indeterminate in catalogue"
    if r.sym not in t26: return "family-level/unlisted symbol (" + r.method + ")"
    return "other: " + r.method
hits["reason"] = hits.apply(reason, axis=1)
print("\n=== B. PER-HIT STATES (all key hits) ===")
print(hits.groupby(["component", "state"]).size().unstack(fill_value=0))
print("\nINDETERMINATE hits by component and reason:"); print(hits[hits.state == "INDETERMINATE"].groupby(["component", "reason"]).size())
print("\nINDETERMINATE hit symbols:"); print(hits[hits.state == "INDETERMINATE"].groupby(["component", "sym", "method"]).size().sort_values(ascending=False).head(40))
# does any INDETERMINATE hit get lost? check that genomes with an indeterminate hit never come out NEGATIVE for that component
bad = 0
for c in ("G1", "G2", "G3", "G4", "G5"):
    gi = set(hits[(hits.component == c) & (hits.state == "INDETERMINATE")].genome); gp = set(hits[(hits.component == c) & (hits.state == "POSITIVE")].genome)
    neg = set(G[G[f"C_{c}"] == "NEGATIVE"].genome); bad += len((gi - gp) & neg)
print("\ngenomes with an indeterminate (and no positive) hit that were classified NEGATIVE:", bad)
print("genomes with internal-stop key hits:", int((hits.method == "INTERNAL_STOP").sum()))
# ---- tables
N = len(G); rows = []
for k in range(1, 6):
    s = G[f"R{k}"]; p, i, n = int((s == "POSITIVE").sum()), int((s == "INDETERMINATE").sum()), int((s == "NEGATIVE").sum())
    rows.append(dict(rule=f"R{k}", components="+".join(f"G{j}" for j in range(1, k + 1)), positive=p, indeterminate=i, negative=n, evaluable=p + n, evaluable_pct=round(100 * (p + n) / N, 1), total=N))
RT = pd.DataFrame(rows); RT.to_csv(os.path.join(REP, "gene_call_rule_counts.csv"), index=False)
print("\n=== C. RULE TABLE (N =", N, ") ==="); print(RT.to_string(index=False))
crow = []
for c in ("G1", "G2", "G3", "G4", "G5"):
    s = G[f"C_{c}"]; crow.append(dict(component=c, positive=int((s == "POSITIVE").sum()), indeterminate=int((s == "INDETERMINATE").sum()), negative=int((s == "NEGATIVE").sum())))
CT = pd.DataFrame(crow); CT.to_csv(os.path.join(REP, "gene_call_component_counts.csv"), index=False)
print("\n=== D. COMPONENT TABLE ==="); print(CT.to_string(index=False))
# reasons per rule (genome level): genomes indeterminate in the rule, by reasons of their indeterminate hits
rr = []
for k in range(1, 6):
    comps = [f"G{j}" for j in range(1, k + 1)]
    ig = set(G[G[f"R{k}"] == "INDETERMINATE"].genome)
    sub = hits[hits.genome.isin(ig) & hits.component.isin(comps) & (hits.state == "INDETERMINATE")]
    cnt = sub.groupby("reason").genome.nunique()
    for r_, n_ in cnt.items(): rr.append(dict(rule=f"R{k}", reason=r_, genomes=int(n_)))
    rr.append(dict(rule=f"R{k}", reason="(total indeterminate genomes)", genomes=len(ig)))
RR = pd.DataFrame(rr); RR.to_csv(os.path.join(REP, "gene_call_indeterminate_reasons.csv"), index=False)
print("\n=== E. INDETERMINATE REASONS PER RULE (genomes can have >1 reason) ==="); print(RR.to_string(index=False))
# allele counts among positives
pos = hits[hits.state == "POSITIVE"]
AT = pos.groupby(["component", "sym"]).genome.nunique().rename("genomes").reset_index().sort_values(["component", "genomes"], ascending=[True, False]); AT.to_csv(os.path.join(REP, "gene_call_positive_alleles.csv"), index=False)
print("\n=== F. POSITIVE ALLELES (genomes) ==="); print(AT.to_string(index=False))
# overlap structure (genotype only)
print("\nG1 positive & G2 positive:", int(((G.C_G1 == "POSITIVE") & (G.C_G2 == "POSITIVE")).sum()), "| G1 pos & G3 pos:", int(((G.C_G1 == "POSITIVE") & (G.C_G3 == "POSITIVE")).sum()), "| G1 neg & G2 pos:", int(((G.C_G1 == "NEGATIVE") & (G.C_G2 == "POSITIVE")).sum()), "| G1 neg & G3 pos:", int(((G.C_G1 == "NEGATIVE") & (G.C_G3 == "POSITIVE")).sum()), "| G1 neg,G2 neg,G3 neg & G4 pos:", int(((G.C_G1 == "NEGATIVE") & (G.C_G2 == "NEGATIVE") & (G.C_G3 == "NEGATIVE") & (G.C_G4 == "POSITIVE")).sum()))
print("isolates indeterminate on >=1 of G1-G3:", int(((G.C_G1 == "INDETERMINATE") | (G.C_G2 == "INDETERMINATE") | (G.C_G3 == "INDETERMINATE")).sum()))
# CTX-M family-level hits
fam = hits[(hits.component == "G1") & (~hits.sym.isin(t26.keys()) | (hits.sym.map(t26) == "INDETERMINATE"))]
print("\nG1 hits with family-level/INDETERMINATE-class CTX-M symbol:", fam.groupby(["sym", "method", "state"]).size().to_dict())
