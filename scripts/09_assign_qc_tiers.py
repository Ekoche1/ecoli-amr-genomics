"""Step 9: decide which genomes to exclude, using genome measurements only.
CONFIRMED problem (excluded): wrong species, mixed culture, failed assembly, or foreign contamination.
SUSPECTED (kept): a large amount of low-coverage, mostly E. coli-like sequence. Everything else: RETAIN.
Reads:  work/qc_table.pkl, work/lowcov.pkl, work/lowcov_classes.pkl, work/genome_16s.jsonl
Writes: work/final_table.pkl, work/genome_qc_all762.csv"""
import os, json
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(ROOT, "work")
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 60); pd.set_option("display.max_colwidth", 120)
d = pd.read_pickle(os.path.join(RES, "qc_table.pkl"))
lc = pd.read_pickle(os.path.join(RES, "lowcov.pkl")); d = d.merge(lc[["id", "lowcov_bp_ge1kb", "highcov_bp_ge1kb"]], on="id")
lcc = pd.read_pickle(os.path.join(RES, "lowcov_classes.pkl")).set_index("id")
s16 = {}
for l in open(os.path.join(RES, "genome_16s.jsonl")):
    r = json.loads(l); s16[r["id"]] = r["regions"]
MEDGC = d.gc.median()
# ---------- 16S call
def call16(regs):
    good = [x for x in regs if x["aln_len"] >= 800 and x["identity"] >= 0.98]
    if not good: return ("none", "")
    groups = sorted({x["best_group"] for x in good})
    esch = [g for g in groups if g in ("Ecoli", "Shigella", "E_fergusonii")]
    other = [g for g in groups if g not in ("Ecoli", "Shigella", "E_fergusonii")]
    desc = "; ".join(f"{x['best_group']} {x['identity']*100:.1f}% ({x['aln_len']} bp)" for x in good)
    if esch and other: return ("both", desc)
    if esch: return ("Escherichia/Shigella", desc)
    return ("non-Escherichia", desc)
c16 = {g: call16(r) for g, r in s16.items()}
d["s16_call"] = d.id.map(lambda g: c16[g][0]); d["s16_detail"] = d.id.map(lambda g: c16[g][1])
# ---------- k-mer species call
def kmer_call(r):
    if r.frac_ecoli_like >= 0.9: return "Escherichia-like"
    oth = {c[6:]: r[c] for c in r.index if c.startswith("cl_other:") and r[c] > 0}
    top = max(oth, key=oth.get) if oth else "unassigned"
    if r.frac_ecoli_like < 0.1: return f"NOT Escherichia (k-mer: {top})"
    return f"MIXED (Escherichia-like {r.frac_ecoli_like:.2f}; other: {top})"
d["kmer_call"] = d.apply(kmer_call, axis=1)
# ---------- concordance of the two independent species lines
cov16 = d.s16_call != "none"
print("=== 16S COVERAGE ===")
print("genomes with a usable 16S region (>=800 bp aligned, >=98% identity):", int(cov16.sum()), "of", len(d), "|", d.s16_call.value_counts().to_dict())
kmer_esch = d.kmer_call == "Escherichia-like"
print("concordance among genomes with 16S: k-mer Escherichia-like & 16S Escherichia:", int((kmer_esch & (d.s16_call == "Escherichia/Shigella")).sum()),
      "| k-mer NOT Escherichia & 16S non-Escherichia:", int((d.kmer_call.str.startswith("NOT") & (d.s16_call == "non-Escherichia")).sum()))
disc = d[cov16 & ((kmer_esch & (d.s16_call == "non-Escherichia")) | (d.kmer_call.str.startswith("NOT") & (d.s16_call == "Escherichia/Shigella")))]
print("DISCORDANT species calls (k-mer vs 16S):", disc[["id", "kmer_call", "s16_call", "s16_detail"]].to_dict("records"))
print("\n16S of every genome that is not 'Escherichia-like' by k-mer:")
print(d[~kmer_esch][["id", "total_bp", "gc", "kmer_call", "s16_call", "s16_detail"]].sort_values("kmer_call").to_string(index=False))
# ---------- tiers
conf = {}
for r in d.itertuples():
    ev = []
    gc_dev = abs(r.gc - MEDGC)
    if r.frac_ecoli_like < 0.1:
        lines = ["k-mer: <10% of assembly resembles E. coli"]
        if gc_dev > 0.03: lines.append(f"GC {r.gc:.3f} vs {MEDGC:.3f}")
        if r.s16_call == "non-Escherichia": lines.append("16S non-Escherichia: " + r.s16_detail)
        lvl = "CONFIRMED" if len(lines) >= 2 else "STRONG (single line)"
        conf[r.id] = (lvl, "wrong species", "; ".join(lines))
    elif r.frac_ecoli_like < 0.9 and r.nonecoli_bp_ge20kb >= 1_000_000 and r.bp_ecoli_like >= 3_000_000:
        lines = [f"k-mer: {r.nonecoli_bp_ge20kb/1e6:.2f} Mb of >=20 kb contigs not E. coli-like plus {r.bp_ecoli_like/1e6:.2f} Mb E. coli-like"]
        if r.total_bp > 1.5 * d.total_bp.median(): lines.append(f"assembly {r.total_bp/1e6:.2f} Mb vs median {d.total_bp.median()/1e6:.2f}")
        if gc_dev > 0.03: lines.append(f"GC {r.gc:.3f}")
        lvl = "CONFIRMED" if len(lines) >= 2 else "STRONG (single line)"
        conf[r.id] = (lvl, "mixed culture", "; ".join(lines))
    elif r.n_contigs > 10000 or r.n50 < 5000:
        conf[r.id] = ("CONFIRMED", "failed assembly", f"{r.n_contigs} contigs, N50 {r.n50} bp, {r.bp_lt5kb/1e6:.2f} Mb in contigs <5 kb")
    elif r.id in lcc.index and lcc.loc[r.id, "other_genus"] >= 1_000_000:
        x = lcc.loc[r.id]
        conf[r.id] = ("CONFIRMED", "foreign contamination (low coverage)", f"{x.other_genus/1e6:.2f} Mb of low-coverage contigs assign to {x.top_other}-like k-mers; their GC {x.gc_lowcov:.3f} vs {MEDGC:.3f}")
    elif r.lowcov_bp_ge1kb >= 400_000:
        x = lcc.loc[r.id] if r.id in lcc.index else None
        det = f"{r.lowcov_bp_ge1kb/1e6:.2f} Mb in low-coverage (<0.3x median) contigs >=1 kb (99th percentile of all genomes: 0.49 Mb)"
        if x is not None: det += f"; {x.frac_ecoli_like*100:.0f}% of it E. coli-like"
        if r.kmer_redundancy > 0.022: det += f"; duplicated-content index {r.kmer_redundancy:.3f} (clean 99th pct 0.022)"
        conf[r.id] = ("SUSPECTED", "possible minor co-contaminant (E. coli-like, low coverage)", det)
d["tier"] = d.id.map(lambda g: conf.get(g, ("RETAIN", "", ""))[0]); d["problem"] = d.id.map(lambda g: conf.get(g, ("", "", ""))[1]); d["evidence"] = d.id.map(lambda g: conf.get(g, ("", "", ""))[2])
# secondary, retain-with-note flags
note = []
for r in d.itertuples():
    n = []
    if r.tier == "RETAIN":
        if r.median_cov_ge10kb > 500: n.append(f"very high coverage ({r.median_cov_ge10kb:.0f}x) with fragmentation (n_contigs {r.n_contigs})")
        if r.n_contigs > 2000: n.append(f"fragmented ({r.n_contigs} contigs)")
        if r.nonecoli_bp_ge20kb >= 250_000: n.append(f"{r.nonecoli_bp_ge20kb/1e3:.0f} kb non-E. coli-like in large contigs (likely plasmids/mobile elements)")
    note.append("; ".join(n))
d["note_retained"] = note
print("\n=== TIER COUNTS (all 762) ===")
print(d.groupby(["tier", "problem"]).size().to_string())
show = d[d.tier != "RETAIN"].sort_values(["tier", "problem", "id"])
print(show[["id", "tier", "problem", "total_bp", "gc", "evidence"]].to_string(index=False))
print("\nRETAINED but noted:", int((d.note_retained != "").sum()), d[d.note_retained != ""][["id", "note_retained"]].head(30).to_string(index=False))
d.to_pickle(os.path.join(RES, "final_table.pkl"))
d.drop(columns=[c for c in d.columns if c.startswith("cl_")]).to_csv(os.path.join(RES, "genome_qc_all762.csv"), index=False)
print("\nfinal table saved:", len(d), "rows")
