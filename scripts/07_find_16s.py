"""Step 7: identify the species of each genome from its 16S rRNA genes, as a second check independent of the k-mer comparison.
16S-like regions are found by 15-mer matches to the 16S genes of the reference genomes, then aligned to them.
Reads:  genomes/, reference_genomes/
Writes: work/genome_16s.jsonl   (resumable)
Usage:  python scripts/07_find_16s.py [n_processes]"""
import os, sys, json, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclib as Q
from Bio import Align

GEN = Q.GENOMES; RES = Q.WORK
K16 = 15
REFS = None; REF15 = None; AL = None
def canon_kmers(seq: bytes, k):
    c = Q.LUT[np.frombuffer(seq, dtype=np.uint8)]; n = len(c)
    if n < k: return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int64)
    invalid = c == 255; cs = np.cumsum(invalid)
    win = cs[k - 1:] - np.concatenate(([0], cs[:n - k])); bad = win > 0
    cc = np.where(invalid, 0, c).astype(np.uint64); m = n - k + 1
    fw = np.zeros(m, dtype=np.uint64); rc = np.zeros(m, dtype=np.uint64)
    for j in range(k):
        x = cc[j:m + j]; fw = (fw << np.uint64(2)) | x
        rc = (rc >> np.uint64(2)) | ((np.uint64(3) - x) << np.uint64(2 * (k - 1)))
    can = np.minimum(fw, rc); pos = np.arange(m)
    return can[~bad], pos[~bad]
def init():
    global REFS, REF15, AL
    REFS = Q.load_16s_refs()
    parts = [canon_kmers(r["seq"].encode(), K16)[0] for r in REFS]
    REF15 = np.unique(np.concatenate(parts))
    AL = Align.PairwiseAligner(); AL.mode = "local"; AL.match_score, AL.mismatch_score, AL.open_gap_score, AL.extend_gap_score = 1, -2, -4, -2
def local_ident(q, r):
    aln = AL.align(q, r)[0]
    m = sum(sum(a == b for a, b in zip(q[s1:e1], r[s2:e2])) for (s1, e1), (s2, e2) in zip(*aln.aligned))
    L = sum(e1 - s1 for (s1, e1), _ in zip(*aln.aligned))
    return m / max(1, L), L
def analyse(gid):
    contigs = Q.contigs_from_text(Q.read_fasta_text(os.path.join(GEN, gid + ".fna")))
    regions = []
    for name, s in contigs:
        if len(s) < 800: continue
        km, pos = canon_kmers(s.encode(), K16)
        if len(km) == 0: continue
        idx = np.searchsorted(REF15, km); idx[idx >= len(REF15)] = len(REF15) - 1
        hp = pos[REF15[idx] == km]
        if len(hp) < 60: continue
        # split hit positions into clusters separated by gaps > 400 bp
        br = np.where(np.diff(hp) > 400)[0] + 1
        for cl in np.split(hp, br):
            if len(cl) < 60: continue
            a, b = max(0, int(cl.min()) - 40), min(len(s), int(cl.max()) + K16 + 40)
            if b - a < 700 or b - a > 2600: continue
            frag = s[a:b]
            for strand, f in (("+", frag), ("-", Q.rc_str(frag))):
                q = Q.kmer_set_str(f)
                top = sorted(((len(q & Q.kmer_set_str(r["seq"])), j) for j, r in enumerate(REFS)), reverse=True)[:5]
                best = (0.0, 0, None, None)
                for _, j in top:
                    ident, L = local_ident(f, REFS[j]["seq"])
                    if L >= 600 and ident > best[0]: best = (ident, L, REFS[j]["group"], REFS[j]["acc"])
                if strand == "+": plus = best
                else: minus = best
            best = plus if plus[0] >= minus[0] else minus
            if best[2] is not None:
                regions.append({"contig": name, "region_len": b - a, "identity": round(best[0], 4), "aln_len": int(best[1]), "best_group": best[2], "best_ref": best[3]})
    return {"id": gid, "n_regions": len(regions), "regions": regions}
def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    out = os.path.join(RES, "genome_16s.jsonl"); done = set()
    if os.path.exists(out):
        for l in open(out):
            try: done.add(json.loads(l)["id"])
            except Exception: pass
    ids = sorted(f[:-4] for f in os.listdir(GEN) if f.endswith(".fna")); todo = [i for i in ids if i not in done]
    print(f"16S pass: {len(ids)} genomes | done {len(done)} | todo {len(todo)}", flush=True); t0 = time.time()
    with Pool(n, initializer=init) as pool, open(out, "a") as fh:
        for k, rec in enumerate(pool.imap_unordered(analyse, todo), 1):
            fh.write(json.dumps(rec, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); fh.flush()
            if k % 20 == 0: print(f"{k}/{len(todo)} | {time.time()-t0:.0f}s", flush=True)
    print("16S DONE", flush=True)
if __name__ == "__main__": main()
