"""Step 5: for each genome, measure how much sequence sits in low-coverage contigs (below 0.3 x the genome's
median coverage), a possible sign of a second organism in the sample.
Reads:  genomes/   (contig coverage is part of each SPAdes contig name)
Writes: work/lowcov.pkl"""
import os, re, json
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); GEN = os.path.join(ROOT, "genomes"); WORK = os.path.join(ROOT, "work")
rows = []
for f in sorted(os.listdir(GEN)):
    if not f.endswith(".fna"): continue
    gid = f[:-4]; lens, covs = [], []
    with open(os.path.join(GEN, f), encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.search(r"length_(\d+)_cov_([\d.]+)", line)
                lens.append(int(m.group(1))); covs.append(float(m.group(2)))
    L = np.array(lens); C = np.array(covs); big = L >= 10000
    med = float(np.median(C[big])) if big.any() else float("nan")
    low = C < 0.3 * med
    rows.append(dict(id=gid, median_cov=round(med, 1), lowcov_bp_ge1kb=int(L[low & (L >= 1000)].sum()), lowcov_bp_lt1kb=int(L[low & (L < 1000)].sum()),
                     lowcov_contigs_ge1kb=int((low & (L >= 1000)).sum()), bp_ge1kb=int(L[L >= 1000].sum()), highcov_bp_ge1kb=int(L[(C > 3 * med) & (L >= 1000)].sum())))
d = pd.DataFrame(rows); d.to_pickle(os.path.join(WORK, "lowcov.pkl"))
print("genomes:", len(d))
print(d[["lowcov_bp_ge1kb", "lowcov_bp_lt1kb", "highcov_bp_ge1kb"]].describe(percentiles=[.5, .9, .95, .99]).round(0).T[["50%", "90%", "95%", "99%", "max"]].to_string())
print("genomes with >= 400 kb low-coverage sequence in contigs >=1 kb:", int((d.lowcov_bp_ge1kb >= 4e5).sum()))
print(d[d.lowcov_bp_ge1kb >= 2e5].sort_values("lowcov_bp_ge1kb", ascending=False).to_string(index=False))
