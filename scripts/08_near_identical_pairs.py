"""Step 8: find pairs of genomes that are near-identical (for example, repeat isolates of the same strain).
Compares the k-mer sketches from step 3 (1,000 k-mers per genome): Jaccard similarity and a Mash-style identity (ANI).
Reads:  work/sketches/
Writes: data/near_identical_pairs.csv   (every pair with Jaccard >= 0.5)"""
import os, sys, math, itertools
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SK = os.path.join(ROOT, "work", "sketches"); C3 = np.uint64(0xD6E8FEB86659FD93)
ids = sorted(f[:-4] for f in os.listdir(SK) if f.endswith(".npy"))
S = 1000
sk = {}
for g in ids:
    a = np.load(os.path.join(SK, g + ".npy"))
    h = np.sort(a * C3)               # hashed values
    sk[g] = h[:S]
print("genomes sketched:", len(sk), "| sketch size:", S)
def jaccard(a, b):
    u = np.union1d(a, b)[:S]                       # bottom-S of the union
    inter = np.intersect1d(a, b, assume_unique=True)
    return float(np.isin(u, inter).mean())
rows = []
arr = list(sk.items())
for i in range(len(arr)):
    gi, ai = arr[i]
    for j in range(i + 1, len(arr)):
        gj, aj = arr[j]
        J = jaccard(ai, aj)
        if J >= 0.5: rows.append((gi, gj, J))
df = pd.DataFrame(rows, columns=["a", "b", "jaccard"])
K = 17
df["mash_ani"] = 1 + (1 / K) * np.log(2 * df.jaccard / (1 + df.jaccard))
df = df.sort_values("jaccard", ascending=False)
df.to_csv(os.path.join(ROOT, "data", "near_identical_pairs.csv"), index=False)
print("pairs with Jaccard >= 0.5:", len(df))
for thr in (0.99, 0.95, 0.9, 0.8, 0.7):
    print(f"  Jaccard >= {thr}: {(df.jaccard >= thr).sum()}")
print(df.head(25).to_string(index=False))
