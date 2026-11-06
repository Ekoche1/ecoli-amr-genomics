"""Step 10: finish the quality-check table and save it.
Adds a plain evidence label for each flagged genome and whether each genome has a lab record. Also prints, for
information, how the exclusions change the number of usable lab results, and summarises near-identical genome pairs.
These printouts do not change which genomes are excluded (that is decided in step 9).
Reads:  work/final_table.pkl, data/meta_data_tidy.txt, data/near_identical_pairs.csv
Writes: data/genome_qc_762.csv, plus work/genome_qc_all762.csv, work/flagged_genomes.csv, work/exclusion_impact.csv"""
import os, json, itertools
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(ROOT, "work"); DATA = os.path.join(ROOT, "data")
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 60); pd.set_option("display.max_colwidth", 110)
d = pd.read_pickle(os.path.join(RES, "final_table.pkl"))
# ---- refine evidence labels
def label(r):
    if r.tier == "CONFIRMED" and r.problem == "wrong species":
        return "CONFIRMED (k-mer + GC + 16S)" if r.s16_call == "non-Escherichia" else "CONFIRMED (k-mer + GC; 16S not recovered)"
    if r.tier == "CONFIRMED" and r.problem == "mixed culture":
        return "CONFIRMED (second genome-sized sequence set + doubled assembly size)"
    if r.tier == "CONFIRMED": return "CONFIRMED"
    if r.tier == "SUSPECTED":
        sig = int(r.lowcov_bp_ge1kb >= 4e5) + int(r.kmer_redundancy > 0.022) + int(r.bp_lt5kb >= 9e5)
        return f"SUSPECTED ({sig} signal{'s' if sig != 1 else ''}: low-coverage excess{', duplicated content' if r.kmer_redundancy > 0.022 else ''}{', short-contig excess' if r.bp_lt5kb >= 9e5 else ''})"
    return "RETAIN"
d["evidence_level"] = d.apply(label, axis=1)
flag = d[d.tier != "RETAIN"].copy()
# ---- phenotype linkage
tidy = pd.read_csv(os.path.join(DATA, "meta_data_tidy.txt"), sep="\t", dtype=str).replace("NA", np.nan)
ph_ids = set(tidy.genome_id); gen_ids = set(d.id)
print("=== ID LINKAGE ===")
print("genomes:", len(gen_ids), "| phenotype rows:", len(ph_ids), "| linked (both):", len(gen_ids & ph_ids), "| genome without phenotype:", len(gen_ids - ph_ids), "| phenotype without genome:", len(ph_ids - gen_ids))
nophen = sorted(gen_ids - ph_ids)
print("genomes without phenotype (26):", nophen)
print("  of which flagged:", sorted(set(nophen) & set(flag.id)))
d["has_phenotype"] = d.id.isin(ph_ids); flag["has_phenotype"] = flag.id.isin(ph_ids)
conf_ids = set(flag[flag.tier == "CONFIRMED"].id); susp_ids = set(flag[flag.tier == "SUSPECTED"].id)
print("\nflagged with phenotype: confirmed", len(conf_ids & ph_ids), "/", len(conf_ids), "| suspected", len(susp_ids & ph_ids), "/", len(susp_ids))
tiers = {"T0 no exclusions": set(), "T1 exclude CONFIRMED": conf_ids, "T2 exclude CONFIRMED + SUSPECTED": conf_ids | susp_ids}
disk = [c for c in tidy.columns if c.startswith("disk_")]; mic = [c for c in tidy.columns if c.startswith("mic_")]
def mic_num(s):
    try: return float(str(s).replace(">", "").replace("<", "").replace("≤", "").replace("≥", "").replace("=", ""))
    except: return np.nan
tidy["caz"] = tidy.mic_CAZ.map(mic_num)
rows = []
for tname, ex in tiers.items():
    t = tidy[~tidy.genome_id.isin(ex)]
    r = {"tier": tname, "genomes_in_set": 762 - len(ex), "isolates_with_phenotype": len(t)}
    for key, cond in [("CPD_disk_R+S", t.disk_CPD.isin(["R", "S"])), ("CPD_R", t.disk_CPD == "R"), ("CPD_S", t.disk_CPD == "S"),
                      ("CAZ_MIC_R(>4)", t.caz > 4), ("CAZ_MIC_S(<=1)", t.caz <= 1),
                      ("CIP_disk_R+S", t.disk_CIP.isin(["R", "S"])), ("CIP_R", t.disk_CIP == "R"), ("GEN_disk_R+S", t.disk_GEN.isin(["R", "S"])), ("GEN_R", t.disk_GEN == "R"),
                      ("MEM_disk_R+S", t.disk_MEM.isin(["R", "S"])), ("MEM_R", t.disk_MEM == "R")]:
        r[key] = int(cond.sum())
    rows.append(r)
imp = pd.DataFrame(rows).set_index("tier").T
print("\n=== USABLE ISOLATES AND OUTCOMES BY EXCLUSION TIER ===")
print(imp.to_string())
# per-drug loss
print("\n=== PER-DRUG AST RESULTS LOST (disk S/I/R calls and MICs) ===")
lost = []
for c in disk + mic:
    t0 = tidy[c].notna().sum(); t1 = tidy[~tidy.genome_id.isin(conf_ids)][c].notna().sum(); t2 = tidy[~tidy.genome_id.isin(conf_ids | susp_ids)][c].notna().sum()
    extra = ""
    if c.startswith("disk_"):
        r0 = (tidy[c] == "R").sum(); r1 = (tidy[~tidy.genome_id.isin(conf_ids)][c] == "R").sum(); extra = f"R: {r0}->{r1}"
    lost.append((c, t0, t1, t2, extra))
print(pd.DataFrame(lost, columns=["outcome", "available_T0", "after_T1", "after_T2", "resistant_calls_T0->T1"]).to_string(index=False))
# phenotypes of the flagged isolates
print("\n=== FLAGGED ISOLATES WITH THEIR CEPHALOSPORIN/QUINOLONE PHENOTYPES ===")
fp = flag.merge(tidy[["genome_id", "genome_name", "disk_CPD", "mic_CAZ", "disk_CIP", "disk_GEN", "disk_MEM"]], left_on="id", right_on="genome_id", how="left")
print(fp[["id", "evidence_level", "problem", "genome_name", "disk_CPD", "mic_CAZ", "disk_CIP", "disk_GEN", "disk_MEM"]].sort_values(["evidence_level", "id"]).to_string(index=False))
# ---- duplicates / clonal clusters (information, not exclusion)
pairs = pd.read_csv(os.path.join(DATA, "near_identical_pairs.csv"))
dd = d.set_index("id")
hi = pairs[pairs.jaccard >= 0.999].copy()
hi["same_size"] = [dd.loc[a, "total_bp"] == dd.loc[b, "total_bp"] and dd.loc[a, "n_contigs"] == dd.loc[b, "n_contigs"] for a, b in zip(hi.a, hi.b)]
print("\n=== NEAR-IDENTICAL GENOME PAIRS (sketch Jaccard >= 0.999, ANI >= 99.997%) ===")
print("pairs:", len(hi), "| of which identical total length AND contig count (possible duplicate file):", int(hi.same_size.sum()))
print(hi[hi.same_size][["a", "b", "jaccard"]].to_string(index=False))
# clusters at J>=0.99
import collections
adj = collections.defaultdict(set)
for a, b in zip(pairs[pairs.jaccard >= 0.99].a, pairs[pairs.jaccard >= 0.99].b): adj[a].add(b); adj[b].add(a)
seen = set(); comps = []
for g in adj:
    if g in seen: continue
    stack = [g]; comp = set()
    while stack:
        x = stack.pop()
        if x in comp: continue
        comp.add(x); stack.extend(adj[x] - comp)
    seen |= comp; comps.append(comp)
sizes = sorted((len(c) for c in comps), reverse=True)
print("clusters of genomes linked at Jaccard>=0.99:", len(comps), "| sizes (top 10):", sizes[:10], "| genomes in any cluster:", sum(sizes))
# ---- save outputs
d.drop(columns=[c for c in d.columns if c.startswith("cl_")]).to_csv(os.path.join(RES, "genome_qc_all762.csv"), index=False)
flag[["id", "tier", "evidence_level", "problem", "has_phenotype", "total_bp", "gc", "n_contigs", "n50", "evidence"]].sort_values(["tier", "id"]).to_csv(os.path.join(RES, "flagged_genomes.csv"), index=False)
imp.to_csv(os.path.join(RES, "exclusion_impact.csv"))
d.to_csv(os.path.join(DATA, "genome_qc_762.csv"), index=False)
print("\nsaved: data/genome_qc_762.csv; work/genome_qc_all762.csv, work/flagged_genomes.csv, work/exclusion_impact.csv")
