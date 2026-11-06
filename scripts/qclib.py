"""Shared functions for the genome quality checks (used by scripts 03, 06 and 07).
Species attribution compares 17-letter DNA words (k-mers) in each genome with 22 reference genomes
(downloaded by script 02). 16S rRNA genes are found and compared with the same references."""
import os, re, glob, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMES = os.path.join(ROOT, "genomes"); REFS = os.path.join(ROOT, "reference_genomes"); WORK = os.path.join(ROOT, "work")
os.makedirs(WORK, exist_ok=True)
K = 17
C1 = np.uint64(0x9E3779B97F4A7C15)       # hash multiplier for the 25% subsample
C2 = np.uint64(0xC2B2AE3D27D4EB4F)       # second multiplier for the sketch
LUT = np.full(256, 255, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"): LUT[_ch] = _i

GROUPS = {
    "Ecoli": ["NC_000913.3", "NC_004431.1", "NC_002695.2", "NC_011750.1", "NC_011751.1", "NC_011415.1", "NC_010468.1"],
    "Shigella": ["NC_004337.2", "NC_007384.1"],
    "E_fergusonii": ["NC_011740.1"],
    "Klebsiella": ["NC_016845.1", "NC_009648.1"],
    "Salmonella": ["NC_003197.2"],
    "Enterobacter": ["NC_014121.1"],
    "Citrobacter": ["NC_009792.1"],
    "Proteus": ["NC_010554.1"],
    "Pseudomonas": ["NC_002516.2", "NC_002947.4"],
    "Acinetobacter": ["NC_009085.1"],
    "Stenotrophomonas": ["NC_010943.1"],
    "Staphylococcus": ["NC_007795.1"],
    "Enterococcus": ["NC_004668.1"],
}
ESCH_CLADE = {"Ecoli", "Shigella", "E_fergusonii"}

def read_fasta_text(path):
    return open(path, encoding="ascii", errors="replace").read()

def contigs_from_text(txt):
    out = []
    for rec in txt.split(">")[1:]:
        head, _, body = rec.partition("\n")
        out.append((head.split()[0], body.replace("\n", "").upper()))
    return out

def kmers_sub(seq: bytes, sub=True):
    """Canonical 17-mers (uint64), optionally the 25% hash-subsample."""
    c = LUT[np.frombuffer(seq, dtype=np.uint8)]
    n = len(c)
    if n < K: return np.empty(0, dtype=np.uint64)
    invalid = (c == 255)
    cs = np.cumsum(invalid)
    win = cs[K - 1:] - np.concatenate(([0], cs[:n - K]))
    bad = win > 0
    cc = np.where(invalid, 0, c).astype(np.uint64)
    m = n - K + 1
    fw = np.zeros(m, dtype=np.uint64); rc = np.zeros(m, dtype=np.uint64)
    for j in range(K):
        x = cc[j:m + j]
        fw = (fw << np.uint64(2)) | x
        rc = (rc >> np.uint64(2)) | ((np.uint64(3) - x) << np.uint64(2 * (K - 1)))
    can = np.minimum(fw, rc)[~bad]
    if sub:
        can = can[((can * C1) >> np.uint64(62)) == 0]
    return can

def load_ref_sets():
    sets = {}
    cache = os.path.join(WORK, "refsets.npz")
    if os.path.exists(cache):
        z = np.load(cache); return {k: z[k] for k in z.files}
    for g, accs in GROUPS.items():
        parts = []
        for a in accs:
            p = os.path.join(REFS, a + ".fna")
            if not os.path.exists(p): continue
            txt = read_fasta_text(p)
            seq = "".join(l for l in txt.split("\n") if l and not l.startswith(">")).upper().encode()
            parts.append(np.unique(kmers_sub(seq)))
        if parts: sets[g] = np.unique(np.concatenate(parts))
    # also K-12 alone for the completeness measure
    p = os.path.join(REFS, "NC_000913.3.fna")
    seq = "".join(l for l in read_fasta_text(p).split("\n") if l and not l.startswith(">")).upper().encode()
    sets["K12"] = np.unique(kmers_sub(seq))
    np.savez(cache, **sets)
    return sets

def share(km, ref):
    if len(km) == 0 or len(ref) == 0: return float("nan")
    idx = np.searchsorted(ref, km); idx[idx >= len(ref)] = len(ref) - 1
    return float((ref[idx] == km).mean())

# ---------------- 16S ----------------
_comp = str.maketrans("ACGT", "TGCA")
rc_str = lambda s: s.translate(_comp)[::-1]
F27 = re.compile(r"AGAGTTTGATC[AC]TGGCTCAG")
def find_16s(contigs, minlen=1400):
    out, seen = [], set()
    for name, seq in contigs:
        if len(seq) < 1500: continue
        for strand, s in (("+", seq), ("-", rc_str(seq))):
            for m in F27.finditer(s):
                frag = s[m.start(): m.start() + 1500]
                if len(frag) >= minlen and "N" not in frag and frag[:300] not in seen:
                    seen.add(frag[:300]); out.append((name, strand, frag))
    return out

def kmer_set_str(s, k=12):
    return {s[i:i + k] for i in range(len(s) - k + 1)}

def load_16s_refs():
    p = os.path.join(WORK, "refs16s.json")
    if os.path.exists(p): return json.load(open(p))
    refs = []
    for g, accs in GROUPS.items():
        for a in accs:
            fp = os.path.join(REFS, a + ".fna")
            if not os.path.exists(fp): continue
            contigs = contigs_from_text(read_fasta_text(fp))
            for i, (_, strand, frag) in enumerate(find_16s(contigs)):
                refs.append({"group": g, "acc": a, "copy": i, "seq": frag})
    json.dump(refs, open(p, "w"))
    return refs

def best_16s(frag, refs, aligner):
    q = kmer_set_str(frag)
    scored = sorted(((len(q & kmer_set_str(r["seq"])), j) for j, r in enumerate(refs)), reverse=True)[:6]
    best = (0.0, None)
    for _, j in scored:
        aln = aligner.align(frag, refs[j]["seq"])[0]
        matches = sum(sum(a == b for a, b in zip(frag[s1:e1], refs[j]["seq"][s2:e2])) for (s1, e1), (s2, e2) in zip(*aln.aligned))
        ident = matches / max(len(frag), len(refs[j]["seq"]))
        if ident > best[0]: best = (ident, refs[j]["group"], refs[j]["acc"])
    return best
