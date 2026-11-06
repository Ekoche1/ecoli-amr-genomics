"""Step 6: for the 11 genomes with at least 200 kb of low-coverage sequence (from step 5), work out what that
low-coverage sequence resembles: E. coli, another genus, or nothing in the reference set.
Reads:  genomes/, reference_genomes/
Writes: work/lowcov_classes.pkl"""
import os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclib as Q
SETS = Q.load_ref_sets()
ids = ["H63", "I42", "H69", "H75", "I53", "B53", "L15", "F60", "G39", "J80", "D50"]
rows = []
for gid in ids:
    contigs = []
    for rec in Q.read_fasta_text(os.path.join(Q.GENOMES, gid + ".fna")).split(">")[1:]:
        head, _, body = rec.partition("\n"); m = re.search(r"length_(\d+)_cov_([\d.]+)", head)
        contigs.append((len(body.replace("\n", "")), float(m.group(2)), body.replace("\n", "").upper()))
    med = np.median([c[1] for c in contigs if c[0] >= 10000])
    cls = {"ecoli_like": 0, "other_genus": 0, "unassigned": 0}; best_g = {}; gcs = []
    for L, cv, s in contigs:
        if L < 1000 or cv >= 0.3 * med: continue
        km = Q.kmers_sub(s.encode()); se = Q.share(km, SETS["Ecoli"])
        if se >= 0.2: cls["ecoli_like"] += L
        else:
            sh = {g: Q.share(km, r) for g, r in SETS.items() if g not in ("Ecoli", "K12")}
            g, v = max(sh.items(), key=lambda kv: kv[1])
            if v >= 0.05: cls["other_genus"] += L; best_g[g] = best_g.get(g, 0) + L
            else: cls["unassigned"] += L
        gcs.append(((s.count("G") + s.count("C")) / L, L))
    tot = sum(cls.values()); gcw = sum(g * L for g, L in gcs) / max(1, sum(L for _, L in gcs))
    rows.append(dict(id=gid, lowcov_bp=tot, **cls, frac_ecoli_like=round(cls["ecoli_like"] / max(1, tot), 3), gc_lowcov=round(gcw, 3), top_other=max(best_g, key=best_g.get) if best_g else "-"))
    print(gid, tot, flush=True)
t = pd.DataFrame(rows); t.to_pickle(os.path.join(Q.WORK, "lowcov_classes.pkl"))
pd.set_option("display.width", 250); print(t.to_string(index=False))
