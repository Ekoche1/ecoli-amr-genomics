"""Step 15: make two alternative versions of the gene calls, used to check the main result.
Some hits are named only to gene-family level (for example "blaCTX-M" with no allele number). The main rules count a
complete family-level hit as POSITIVE. The two alternatives count it as INDETERMINATE instead:
  S10a: family-level blaCTX-M hits
  S10b: family-level blaCTX-M and blaVIM hits
Everything else is unchanged.
Reads:  data/gene_calls_749.csv, data/gene_hits_749.csv
Writes: data/gene_calls_749_S10a.csv, data/gene_calls_749_S10b.csv
Usage:  python scripts/15_alternative_gene_calls.py"""
import os
import pandas as pd
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_rows", 100)
Hs = pd.read_csv(os.path.join(DATA, "gene_hits_749.csv"), dtype={"genome": str})
G = pd.read_csv(os.path.join(DATA, "gene_calls_749.csv"), dtype={"genome": str})
def s10(G, symbols):
    """A hit whose symbol is in `symbols` (family-level, no allele number) is INDETERMINATE for its component even if complete (BLASTX, >=90% coverage). Everything else unchanged."""
    H2 = Hs.copy(); m = H2.sym.isin(symbols) & (H2.state == "POSITIVE"); H2.loc[m, "state"] = "INDETERMINATE"
    comp = {}
    for c in ("G1", "G2", "G3", "G4", "G5"):
        sub = H2[H2.component == c]; pos = set(sub[sub.state == "POSITIVE"].genome); ind = set(sub[sub.state == "INDETERMINATE"].genome)
        comp[c] = {g: ("POSITIVE" if g in pos else ("INDETERMINATE" if g in ind else "NEGATIVE")) for g in G.genome}
    out = G[["genome"]].copy()
    for c in comp: out[f"C_{c}"] = out.genome.map(comp[c])
    for k in range(1, 6):
        cols = [f"C_G{i}" for i in range(1, k + 1)]
        out[f"R{k}"] = out[cols].apply(lambda r: "POSITIVE" if "POSITIVE" in r.values else ("INDETERMINATE" if "INDETERMINATE" in r.values else "NEGATIVE"), axis=1)
    return out
for label, name, syms in (("S10a: family-level blaCTX-M hits counted as indeterminate", "S10a", ["blaCTX-M"]),
                          ("S10b: family-level blaCTX-M and blaVIM hits counted as indeterminate", "S10b", ["blaCTX-M", "blaVIM"])):
    S = s10(G, syms); print("\n====", label)
    cols = [f"C_G{i}" for i in range(1, 6)] + [f"R{k}" for k in range(1, 6)]
    # check: with no symbols changed, the per-hit table gives back exactly the main gene calls
    P = s10(G, [])
    print("main gene calls re-derived from the per-hit table match data/gene_calls_749.csv in all 10 columns:", all((P[c] == G[c]).all() for c in cols))
    for c in cols:
        cnt = lambda d: (int((d[c] == "POSITIVE").sum()), int((d[c] == "INDETERMINATE").sum()), int((d[c] == "NEGATIVE").sum()))
        ch = G[c] != S[c]
        print(f"  {c:5s} main P/I/N {cnt(G)} -> {name} {cnt(S)} | genomes changed: {S[ch].genome.tolist()} ({G[ch][c].tolist()} -> {S[ch][c].tolist()})")
    S.to_csv(os.path.join(DATA, f"gene_calls_749_{name}.csv"), index=False)
    ev = lambda d, k: d[f"R{k}"].isin(["POSITIVE", "NEGATIVE"])
    print("  evaluable R1:", int(ev(G, 1).sum()), "->", int(ev(S, 1).sum()), "| R4:", int(ev(G, 4).sum()), "->", int(ev(S, 4).sum()), "| common R1&R4:", int((ev(G, 1) & ev(G, 4)).sum()), "->", int((ev(S, 1) & ev(S, 4)).sum()))
for g in ("E71", "H12", "I76"):
    print(g, G[G.genome == g].iloc[0].to_dict())
