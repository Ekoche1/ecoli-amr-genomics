"""Step 16 (main analysis): link each genome's gene calls to its ceftazidime MIC and measure how well each rule
identifies resistance: sensitivity, specificity, PPV and NPV, with 95% confidence intervals.
Ceftazidime MIC (mg/L, EUCAST): resistant > 4, susceptible <= 1; 1 < MIC <= 4 (intermediate) is kept separate.
INDETERMINATE gene calls are set aside, never counted as negative. Near-identical genomes (Jaccard >= 0.99) are
grouped so that the bootstrap resamples them together.
Reads:  data/gene_calls_749.csv, data/gene_calls_749_S10a.csv, data/gene_calls_749_S10b.csv, data/meta_data_tidy.txt,
        data/genome_qc_762.csv, data/near_identical_pairs.csv
Writes: results/linked_isolates_749.csv, results/reconciliation_state_by_caz_category.csv, results/results_all_sets_long.csv,
        results/paired_contrasts.csv, results/analysis_log.txt, results/analysis_run_record.json
Usage:  python scripts/16_link_and_analyse.py"""
import os, re, sys, json, hashlib, collections
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data"); T = os.path.join(ROOT, "results"); os.makedirs(T, exist_ok=True)
SEED = 20261018; B = 10000   # fixed random seed and number of bootstrap replicates, so every run gives the same answer
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
log = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)
# ---------------- inputs
IN = {"gene_calls": os.path.join(DATA, "gene_calls_749.csv"), "gene_calls_S10a": os.path.join(DATA, "gene_calls_749_S10a.csv"),
      "gene_calls_S10b": os.path.join(DATA, "gene_calls_749_S10b.csv"), "lab_results": os.path.join(DATA, "meta_data_tidy.txt"),
      "genome_qc": os.path.join(DATA, "genome_qc_762.csv"), "near_identical_pairs": os.path.join(DATA, "near_identical_pairs.csv")}
G = pd.read_csv(IN["gene_calls"], dtype=str).set_index("genome")
A = pd.read_csv(IN["gene_calls_S10a"], dtype=str).set_index("genome")
Bx = pd.read_csv(IN["gene_calls_S10b"], dtype=str).set_index("genome")
tidy = pd.read_csv(IN["lab_results"], sep="\t", dtype=str, keep_default_na=False)[["genome_id", "mic_CAZ"]].set_index("genome_id")
d = pd.read_csv(IN["genome_qc"], dtype={"id": str}).set_index("id")
pairs = pd.read_csv(IN["near_identical_pairs"], dtype={"a": str, "b": str})
retained = sorted(d[d.tier != "CONFIRMED"].index)
V = []   # validations
def check(name, ok, detail=""):
    V.append((name, bool(ok), detail)); say(("PASS " if ok else "FAIL ") + name + (" | " + str(detail) if detail else ""))
check("V1 genotype table covers exactly the 749 retained genomes", set(G.index) == set(retained) and len(G) == 749)
exp = {"R1": (225, 3, 521), "R2": (232, 24, 493), "R3": (248, 23, 478), "R4": (280, 27, 442), "R5": (291, 27, 431)}
check("V2 rule-state counts as expected", all((int((G[k] == "POSITIVE").sum()), int((G[k] == "INDETERMINATE").sum()), int((G[k] == "NEGATIVE").sum())) == v for k, v in exp.items()))
check("V3 S10 tables share the genome set and differ from the main calls only as expected", set(A.index) == set(G.index) == set(Bx.index) and int((G.R1 != A.R1).sum()) == 2 and int((G.R5 != Bx.R5).sum()) == 3)
# ---------------- ceftazidime category (EUCAST: R > 4, S <= 1, 1 < MIC <= 4 intermediate kept separate; ">x" with x >= 4 is R, "<=x" or "<x" with x <= 1 is S)
def cat(v):
    v = v.strip()
    if v in ("NA", ""): return "missing", np.nan
    m = re.match(r"^(>=|<=|>|<|≤|≥)?\s*([0-9.]+)$", v)
    if not m: raise ValueError("unparsed MIC: " + v)
    op, x = m.group(1) or "", float(m.group(2))
    if op in (">", ">="):
        if x >= 4: return "R", x
        raise ValueError("ambiguous censored MIC " + v)
    if op in ("<", "<=", "≤"):
        if x <= 1: return "S", x
        raise ValueError("ambiguous censored MIC " + v)
    return ("R" if x > 4 else ("S" if x <= 1 else "I")), x
L = pd.DataFrame(index=retained); L.index.name = "genome_id"
L["in_pheno_file"] = [i in tidy.index for i in L.index]
L["caz_mic_raw"] = [tidy.mic_CAZ.get(i, "NO_RECORD") for i in L.index]
cc = [cat(v) if v != "NO_RECORD" else ("no_record", np.nan) for v in L.caz_mic_raw]
L["caz_category"] = [c[0] for c in cc]; L["caz_mic_numeric"] = [c[1] for c in cc]
L["primary_analysable"] = L.caz_category.isin(["R", "S"])
for c in ["C_G1", "C_G2", "C_G3", "C_G4", "C_G5", "R1", "R2", "R3", "R4", "R5"]: L[c] = G.loc[L.index, c]
L["S10a_R1"] = A.loc[L.index, "R1"]; L["S10a_R4"] = A.loc[L.index, "R4"]; L["S10a_R5"] = A.loc[L.index, "R5"]; L["S10b_R5"] = Bx.loc[L.index, "R5"]
L["qc_tier"] = d.loc[L.index, "tier"]
# five isolates whose lab results disagree with their near-identical genomes on 3 or more antibiotics: flagged for information only, kept in the analysis
cd = {"I81", "E4", "I12", "F62", "H79"}; L["clone_discordant_flag"] = L.index.isin(cd)
# ---------------- resampling units: single-linkage clusters of retained genomes at sketch Jaccard >= 0.99
rs = set(retained); e = pairs[(pairs.jaccard >= 0.99) & pairs.a.isin(rs) & pairs.b.isin(rs)]
par = {g: g for g in retained}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a_, b_ in zip(e.a, e.b): par[f(a_)] = f(b_)
root = {g: f(g) for g in retained}; size = collections.Counter(root.values())
L["unit_id"] = [root[g] for g in L.index]; L["unit_size"] = [size[root[g]] for g in L.index]
multi = [r for r, s in size.items() if s > 1]
say(f"resampling units among 749 retained genomes: {len(size)} (multi-genome clusters {len(multi)} holding {sum(size[r] for r in multi)} genomes; singletons {sum(1 for s in size.values() if s == 1)}); max cluster size {max(size.values())}; lowest within-cluster pair Jaccard {e.jaccard.min():.4f}")
say("cluster size distribution:", dict(sorted(collections.Counter(size[r] for r in multi).items())))
L.to_csv(os.path.join(T, "linked_isolates_749.csv"))
# ---------------- reconciliation
say("\n=== RECONCILIATION (primary set) ===")
say("genomes in primary set:", len(L), "| with a phenotype record:", int(L.in_pheno_file.sum()), "| without record:", int((~L.in_pheno_file).sum()))
rc = L.caz_category.value_counts().to_dict(); say("ceftazidime MIC category among the 749:", rc)
check("V4 lab records matched = 724", int(L.in_pheno_file.sum()) == 724)
check("V5 CAZ missing 14, intermediate 66, R 207, S 437", (rc.get("missing", 0), rc.get("I", 0), rc.get("R", 0), rc.get("S", 0)) == (14, 66, 207, 437), rc)
rec = []
for k in ["R1", "R2", "R3", "R4", "R5"]:
    ct = pd.crosstab(L[k], L.caz_category).reindex(index=["POSITIVE", "INDETERMINATE", "NEGATIVE"], columns=["R", "I", "S", "missing", "no_record"], fill_value=0); ct.insert(0, "rule", k); rec.append(ct.reset_index().rename(columns={k: "genotype_state"}))
REC = pd.concat(rec); REC.to_csv(os.path.join(T, "reconciliation_state_by_caz_category.csv"), index=False); say(REC.to_string(index=False))
# ---------------- measures
def wilson(x, n, z=1.959964):
    if n == 0: return (np.nan, np.nan)
    p = x / n; den = 1 + z * z / n; c = (p + z * z / (2 * n)) / den; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den; return (max(0, c - h), min(1, c + h))
units = sorted(size); uidx = {u: i for i, u in enumerate(units)}; ui = np.array([uidx[root[g]] for g in L.index])
rng = np.random.default_rng(SEED); DRAW = [rng.integers(0, len(units), size=len(units)) for _ in range(B)]
def counts_per_unit(mask, pos, res):
    """per-unit TP,FN,FP,TN for isolates in mask. pos: bool array rule positive; res: bool array phenotype R."""
    out = np.zeros((len(units), 4))
    for j, (sel) in enumerate([mask & pos & res, mask & ~pos & res, mask & pos & ~res, mask & ~pos & ~res]): np.add.at(out[:, j], ui[sel], 1)
    return out
def boot_stats(cu):
    tot = np.zeros((B, 4))
    for b in range(B): tot[b] = cu[DRAW[b]].sum(axis=0)
    return tot
def metr(t):
    tp, fn, fp, tn = t.T if t.ndim == 2 else t
    with np.errstate(all="ignore"): return dict(sens=tp / (tp + fn), spec=tn / (tn + fp), ppv=tp / (tp + fp), npv=tn / (tn + fn))
def pct(x): x = x[~np.isnan(x)]; return (float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))) if len(x) else (np.nan, np.nan)
res_ = (L.caz_category == "R").values; ana = L.primary_analysable.values
rows = []
def one(label, ks, evalmask):
    """measures for rule list ks on the isolate mask evalmask (same isolates for all rules in ks)"""
    for k in ks:
        pos = (L[k] == "POSITIVE").values; m = evalmask; cu = counts_per_unit(m, pos, res_); t = cu.sum(axis=0); tp, fn, fp, tn = t; bs = boot_stats(cu); bm = metr(bs); pm = metr(t)
        n_ind = int(((L[k] == "INDETERMINATE") & ana).sum()); n_ind_r = int(((L[k] == "INDETERMINATE") & ana & res_).sum())
        r = dict(set=label, rule=k, evaluable_n=int(t.sum()), evaluable_R=int(tp + fn), evaluable_S=int(fp + tn), TP=int(tp), FN=int(fn), FP=int(fp), TN=int(tn), indeterminate_with_phenotype=n_ind, indet_R=n_ind_r, indet_S=n_ind - n_ind_r)
        for nm, (x, nn) in dict(sens=(tp, tp + fn), spec=(tn, tn + fp), ppv=(tp, tp + fp), npv=(tn, tn + fn)).items():
            w = wilson(x, nn); b_ = pct(bm[nm]); r.update({nm: float(pm[nm]), nm + "_wilson_lo": w[0], nm + "_wilson_hi": w[1], nm + "_boot_lo": b_[0], nm + "_boot_hi": b_[1]})
        rows.append(r)
ev = {k: ana & (L[k] != "INDETERMINATE").values for k in ["R1", "R2", "R3", "R4", "R5"]}
for k in ["R1", "R2", "R3", "R4", "R5"]: one("own evaluable set", [k], ev[k])
C14 = ev["R1"] & ev["R4"]; C1234 = ev["R1"] & ev["R2"] & ev["R3"] & ev["R4"]
one("common set C(R1,R4)", ["R1", "R4"], C14); one("common set, evaluable for all of R1-R4", ["R1", "R2", "R3", "R4"], C1234)
odd = list(L.index[C14 & ~C1234])
check("V6 four-way common set is a subset of C(R1,R4); difference explained by G2-indeterminate isolates made positive by G3/G4", bool(C1234.sum() <= C14.sum()) and all(L.loc[i, "C_G2"] == "INDETERMINATE" and L.loc[i, "R4"] == "POSITIVE" for i in odd), f"C(R1,R4)={int(C14.sum())}, four-way={int(C1234.sum())}, differing isolates={odd}")
one("common set C(R4,R5)", ["R4", "R5"], ev["R4"] & ev["R5"])
R = pd.DataFrame(rows); R.to_csv(os.path.join(T, "results_all_sets_long.csv"), index=False)
# paired contrasts on common evaluable sets
def contrast(j, k):
    m = ev[j] & ev[k]; pj = (L[j] == "POSITIVE").values; pk = (L[k] == "POSITIVE").values
    cj = counts_per_unit(m, pj, res_); ck = counts_per_unit(m, pk, res_); bj, bk = boot_stats(cj), boot_stats(ck); mj, mk = metr(bj), metr(bk); tj, tk = metr(cj.sum(axis=0)), metr(ck.sum(axis=0))
    ds, dp = mk["sens"] - mj["sens"], mk["spec"] - mj["spec"]; sj, sk = metr(cj.sum(axis=0)), metr(ck.sum(axis=0))
    new_r = int((m & ~pj & pk & res_).sum()); new_s_flag = int((m & ~pj & pk & ~res_).sum()); lost = int((m & pj & ~pk).sum())
    p_sens = min(1.0, 2 * 0.5 ** new_r) if new_r > 0 else 1.0; p_spec = min(1.0, 2 * 0.5 ** new_s_flag) if new_s_flag > 0 else 1.0
    return dict(contrast=f"{j} vs {k}", common_n=int(m.sum()), common_R=int((m & res_).sum()), common_S=int((m & ~res_).sum()), dropped_from_j_own_set=int(ev[j].sum() - m.sum()), dropped_from_k_own_set=int(ev[k].sum() - m.sum()),
                sens_j=float(sj["sens"]), sens_k=float(sk["sens"]), d_sens_pp=100 * float(sk["sens"] - sj["sens"]), d_sens_boot_lo_pp=100 * pct(ds)[0], d_sens_boot_hi_pp=100 * pct(ds)[1],
                spec_j=float(sj["spec"]), spec_k=float(sk["spec"]), d_spec_pp=100 * float(sk["spec"] - sj["spec"]), d_spec_boot_lo_pp=100 * pct(dp)[0], d_spec_boot_hi_pp=100 * pct(dp)[1],
                resistant_newly_detected=new_r, susceptible_newly_flagged=new_s_flag, positives_lost=lost, exact_mcnemar_p_sens_unadj=p_sens, exact_mcnemar_p_spec_unadj=p_spec)
CT = pd.DataFrame([contrast("R1", "R4"), contrast("R1", "R2"), contrast("R2", "R3"), contrast("R3", "R4"), contrast("R4", "R5")])
check("V7 nesting holds (no positive lost between successive/primary rules)", int(CT.positives_lost.sum()) == 0)
CT.to_csv(os.path.join(T, "paired_contrasts.csv"), index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
say("\n=== MEASURES ===")
show = R[["set", "rule", "evaluable_n", "evaluable_R", "evaluable_S", "TP", "FN", "FP", "TN", "indeterminate_with_phenotype", "indet_R", "indet_S"]]; say(show.to_string(index=False))
for _, r in R.iterrows():
    say(f"{r['set']:22s} {r.rule}: sens {r.sens:.3f} (Wilson {r.sens_wilson_lo:.3f}-{r.sens_wilson_hi:.3f}; boot {r.sens_boot_lo:.3f}-{r.sens_boot_hi:.3f}) | spec {r.spec:.3f} (Wilson {r.spec_wilson_lo:.3f}-{r.spec_wilson_hi:.3f}; boot {r.spec_boot_lo:.3f}-{r.spec_boot_hi:.3f}) | PPV {r.ppv:.3f} NPV {r.npv:.3f}")
say("\n=== PAIRED CONTRASTS ==="); say(CT.round(4).T.to_string())
# Holm across the four plan contrasts (R1vR4, R1vR2, R2vR3, R3vR4)
pv = CT.iloc[:4]; ps = [("sens", pv.exact_mcnemar_p_sens_unadj.values), ("spec", pv.exact_mcnemar_p_spec_unadj.values)]
for nm, p in ps:
    o = np.argsort(p); adj = np.empty(4); run = 0
    for rnk, i in enumerate(o): run = max(run, min(1, (4 - rnk) * p[i])); adj[i] = run
    say(f"Holm-adjusted exact McNemar p ({nm}) for [R1vR4, R1vR2, R2vR3, R3vR4]:", np.round(adj, 5).tolist())
# bootstrap stability check (Monte-Carlo error only): second seed, headline rows
rng2 = np.random.default_rng(1); DRAW2 = [rng2.integers(0, len(units), size=len(units)) for _ in range(B)]
def chk(k):
    pos = (L[k] == "POSITIVE").values; cu = counts_per_unit(ev["R1"] & ev["R4"], pos, res_); tot = np.array([cu[x].sum(axis=0) for x in DRAW2]); m = metr(tot); return pct(m["sens"]), pct(m["spec"])
say("\nSecond-seed check on C(R1,R4):", {k: chk(k) for k in ["R1", "R4"]})
# the three genomes whose calls change in S10a/S10b
mv = L.loc[["E71", "H12", "I76"], ["caz_mic_raw", "caz_category", "R1", "S10a_R1", "R4", "S10a_R4", "R5", "S10b_R5"]]; say("\nS10 moved genomes and their phenotype:\n" + mv.to_string())
say("\nclone-discordant isolates (retained):", L.loc[sorted(cd), ["caz_mic_raw", "caz_category", "R1", "R4"]].to_string())
say("\nValidation summary:", "ALL PASS" if all(v[1] for v in V) else "FAILURES PRESENT", [v[0] for v in V if not v[1]])
open(os.path.join(T, "analysis_log.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
json.dump(dict(seed=SEED, bootstrap_replicates=B, validations=[dict(name=n, ok=o, detail=str(dd)) for n, o, dd in V], inputs={os.path.basename(p): sha(p) for p in IN.values()}), open(os.path.join(T, "analysis_run_record.json"), "w"), indent=1)
