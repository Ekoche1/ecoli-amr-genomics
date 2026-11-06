"""Step 4: combine the per-genome measurements into one table and print their distributions.
Reads:  work/genome_qc.jsonl
Writes: work/qc_table.pkl"""
import json, os, pandas as pd, numpy as np
WORK = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "work")
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
rows = [json.loads(l) for l in open(os.path.join(WORK, "genome_qc.jsonl"))]
d = pd.DataFrame(rows)
cb = pd.json_normalize(d.class_bp).fillna(0); cb.columns = ["cl_" + c for c in cb.columns]
d = pd.concat([d.drop(columns="class_bp"), cb], axis=1)
ecl = [c for c in d.columns if c in ("cl_Ecoli_core", "cl_Ecoli_accessory", "cl_EschRelated")]
oth = [c for c in d.columns if c.startswith("cl_other:")]
d["bp_ecoli_like"] = d[ecl].sum(axis=1); d["bp_other_genus"] = d[oth].sum(axis=1); d["bp_unassigned"] = d.get("cl_unassigned", 0)
d["bp_ge5kb"] = d.bp_ecoli_like + d.bp_other_genus + d.bp_unassigned
d["frac_ecoli_like"] = d.bp_ecoli_like / d.bp_ge5kb
print("genomes:", len(d), "| unique ids:", d.id.nunique())
def mad_out(col, k=5):
    x = d[col].astype(float); med = x.median(); mad = (x - med).abs().median() * 1.4826
    return med, mad
for col in ["total_bp", "n_contigs", "n50", "gc", "bp_lt5kb", "kmer_redundancy", "k12_recovery", "frac_bp_gc_dev_ge20kb", "frac_bp_cov_dev_ge10kb", "median_cov_ge10kb", "nonecoli_bp_ge20kb", "bp_unassigned", "bp_other_genus", "n_bases"]:
    med, mad = mad_out(col); x = d[col].astype(float)
    print(f"{col:26s} median {med:12.4f} robustSD {mad:12.4f} | min {x.min():12.4f} p1 {x.quantile(.01):12.4f} p99 {x.quantile(.99):12.4f} max {x.max():12.4f}")
print("\nfrac_ecoli_like distribution:", d.frac_ecoli_like.describe(percentiles=[.001,.01,.05,.5]).round(3).to_dict())
print("\nGenomes with frac_ecoli_like < 0.8:")
print(d[d.frac_ecoli_like < 0.8][["id", "total_bp", "gc", "frac_ecoli_like", "bp_ecoli_like", "bp_other_genus", "bp_unassigned", "n_contigs", "n50"]].sort_values("frac_ecoli_like").to_string(index=False))
print("\nnon-E. coli bp in contigs >=20 kb, top 40:")
cols = ["id", "total_bp", "gc", "nonecoli_bp_ge20kb", "nonecoli_largest_contig", "bp_ecoli_like", "bp_other_genus", "bp_unassigned", "frac_bp_gc_dev_ge20kb", "kmer_redundancy", "n_contigs", "n50"] + oth
print(d.sort_values("nonecoli_bp_ge20kb", ascending=False)[cols].head(40).to_string(index=False))
d.to_pickle(os.path.join(WORK, "qc_table.pkl"))
