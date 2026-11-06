"""Step 3: measure the quality of every genome assembly.
For each genome: size, number of contigs, N50, GC content, coverage, how much of the sequence looks like E. coli
(by k-mer comparison with the reference genomes), duplicated content, completeness against E. coli K-12,
and the species of any 16S rRNA genes found. Also saves a small k-mer sketch of each genome for script 08.
Reads:  genomes/, reference_genomes/
Writes: work/genome_qc.jsonl, work/sketches/   (resumable)
Usage:  python scripts/03_measure_genomes.py [n_processes]"""
import os, sys, re, json, time
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclib as Q
from Bio import Align

GEN = Q.GENOMES; RES = Q.WORK; SKT = os.path.join(Q.WORK, "sketches")
os.makedirs(RES, exist_ok=True); os.makedirs(SKT, exist_ok=True)
SETS = None; REFS16 = None; ALIGNER = None
def init():
    global SETS, REFS16, ALIGNER
    SETS = Q.load_ref_sets(); REFS16 = Q.load_16s_refs()
    ALIGNER = Align.PairwiseAligner(); ALIGNER.mode = "global"
    ALIGNER.match_score, ALIGNER.mismatch_score, ALIGNER.open_gap_score, ALIGNER.extend_gap_score = 1, -1, -3, -1

def analyse(gid):
    txt = Q.read_fasta_text(os.path.join(GEN, gid + ".fna"))
    raw = [r for r in txt.split(">")[1:]]
    contigs, covs = [], []
    for rec in raw:
        head, _, body = rec.partition("\n")
        m = re.search(r"length_(\d+)_cov_([\d.]+)", head)
        contigs.append((head.split()[0], body.replace("\n", "").upper()))
        covs.append(float(m.group(2)) if m else np.nan)
    lens = np.array([len(s) for _, s in contigs]); order = np.argsort(-lens)
    tot = int(lens.sum()); cum = np.cumsum(lens[order]); n50 = int(lens[order][np.searchsorted(cum, tot / 2)]); l50 = int(np.searchsorted(cum, tot / 2) + 1)
    allseq = "".join(s for _, s in contigs); nN = allseq.count("N")
    gc = (allseq.count("G") + allseq.count("C")) / max(1, tot - nN)
    covs = np.array(covs)
    big = lens >= 10000
    med_cov = float(np.median(covs[big])) if big.any() and not np.isnan(covs[big]).all() else float("nan")
    rec = dict(id=gid, total_bp=tot, n_contigs=len(lens), n50=n50, l50=l50, largest=int(lens.max()), gc=round(gc, 4), n_bases=nN,
               contigs_ge1kb=int((lens >= 1000).sum()), bp_lt1kb=int(lens[lens < 1000].sum()), bp_lt5kb=int(lens[lens < 5000].sum()),
               median_cov_ge10kb=round(med_cov, 1))
    # ---- per-contig similarity (contigs >= 5 kb) and GC heterogeneity (>= 20 kb)
    cls = {}; nonec_big = 0; nonec_largest = 0; unl_bp = 0; gc_dev_bp = 0; bp20 = 0; cov_dev_bp = 0; bp10 = 0
    nonec_gcs = []
    for (name, s), L, cv in zip(contigs, lens, covs):
        if L >= 20000:
            g = (s.count("G") + s.count("C")) / max(1, L - s.count("N")); bp20 += L
            if abs(g - 0.508) > 0.05: gc_dev_bp += L
        if L >= 10000 and med_cov == med_cov and cv == cv:
            bp10 += L
            if cv < 0.3 * med_cov or cv > 3 * med_cov: cov_dev_bp += L
        if L < 5000: continue
        km = Q.kmers_sub(s.encode())
        se = Q.share(km, SETS["Ecoli"])
        if se >= 0.5: lab = "Ecoli_core"
        elif se >= 0.2: lab = "Ecoli_accessory"
        else:
            best, bs = None, 0.0
            for grp, ref in SETS.items():
                if grp in ("Ecoli", "K12"): continue
                sh = Q.share(km, ref)
                if sh > bs: best, bs = grp, sh
            if best is not None and bs >= 0.05 and bs > se:
                lab = "EschRelated" if best in Q.ESCH_CLADE else "other:" + best
            else: lab = "unassigned"
            if L >= 20000 and lab != "EschRelated":
                nonec_big += L; nonec_largest = max(nonec_largest, int(L))
                nonec_gcs.append(((s.count("G") + s.count("C")) / max(1, L - s.count("N")), L))
        cls[lab] = cls.get(lab, 0) + int(L)
    rec["class_bp"] = cls
    rec["nonecoli_bp_ge20kb"] = nonec_big; rec["nonecoli_largest_contig"] = nonec_largest
    rec["frac_bp_gc_dev_ge20kb"] = round(gc_dev_bp / max(1, bp20), 4)
    rec["frac_bp_cov_dev_ge10kb"] = round(cov_dev_bp / max(1, bp10), 4) if bp10 else None
    # ---- k-mer redundancy, K-12 core completeness, sketch (contigs >= 1 kb)
    kms = [Q.kmers_sub(s.encode()) for _, s in contigs if len(s) >= 1000]
    km_all = np.concatenate(kms) if kms else np.empty(0, dtype=np.uint64)
    uniq = np.unique(km_all)
    rec["kmer_redundancy"] = round(1 - len(uniq) / max(1, len(km_all)), 4)
    k12 = SETS["K12"]; idx = np.searchsorted(uniq, k12); idx[idx >= len(uniq)] = max(0, len(uniq) - 1)
    rec["k12_recovery"] = round(float((uniq[idx] == k12).mean()), 4) if len(uniq) else 0.0
    np.save(os.path.join(SKT, gid + ".npy"), uniq[((uniq * Q.C2) >> np.uint64(57)) == 0])
    # ---- 16S
    hits = Q.find_16s(contigs)
    rec["n_16s_regions"] = len(hits); s16 = []
    for cname, strand, frag in hits[:6]:
        ident, grp, acc = Q.best_16s(frag, REFS16, ALIGNER)
        s16.append({"contig": cname, "identity": round(ident, 4), "best_group": grp, "best_ref": acc})
    rec["s16"] = s16
    return rec

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    out = os.path.join(RES, "genome_qc.jsonl")
    done = set()
    if os.path.exists(out):
        for l in open(out):
            try: done.add(json.loads(l)["id"])
            except Exception: pass
    ids = sorted(f[:-4] for f in os.listdir(GEN) if f.endswith(".fna"))
    todo = [i for i in ids if i not in done]
    print(f"genomes: {len(ids)} | already done: {len(done)} | to do: {len(todo)}", flush=True)
    t0 = time.time()
    with Pool(n, initializer=init) as pool, open(out, "a") as fh:
        for k, rec in enumerate(pool.imap_unordered(analyse, todo), 1):
            fh.write(json.dumps(rec, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n"); fh.flush()
            if k % 10 == 0: print(f"{k}/{len(todo)} | {time.time()-t0:.0f}s", flush=True)
    print("QC DONE", flush=True)
if __name__ == "__main__":
    main()
