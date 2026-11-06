"""The gene-calling rules (used by scripts 14 and 15).
A hit is POSITIVE when it is complete (AMRFinderPlus method EXACTX, ALLELEX or BLASTX, and >= 90% coverage) and its allele
belongs to a gene group; INDETERMINATE when it is partial or the allele's effect is unclear. Gene groups (components):
G1 CTX-M; G2 other ESBL; G3 acquired AmpC-family; G4 ampC promoter variants (C-42T, T-32A); G5 carbapenemase families.
Rules are cumulative: G1, G1+G2, G1-G3, G1-G4, and G1-G5 (supplementary)."""
import os, re
import pandas as pd
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
COMPLETE_METHODS = {"EXACTX", "ALLELEX", "BLASTX"}
COV_MIN = 90.0
FAM = re.compile(r"^bla(CTX-M|SHV|TEM|OXA|GES|PER|VEB|CMY|DHA|ACT|FOX|MOX|MIR|ACC|CFE|LAT|MOR|NDM|KPC|VIM|IMP|CARB|LAP)")
FAM_G1 = {"blaCTX-M"}
FAM_G3 = {"blaCMY", "blaDHA", "blaACT", "blaFOX", "blaMOX", "blaMIR", "blaACC", "blaCFE", "blaLAT", "blaMOR"}
FAM_G5 = {"blaNDM", "blaKPC", "blaVIM", "blaIMP"}
FAM_G2_AMBIG = {"blaSHV", "blaTEM", "blaOXA", "blaGES", "blaPER", "blaVEB"}       # families containing both ESBL and non-ESBL alleles
PROM_INCLUDED = {"ampC_C-42T", "ampC_T-32A"}
CLASS_TO_COMP = {"G1_CTXM_ESBL": "G1", "G2_other_ESBL": "G2", "G3_acquired_AmpC_family": "G3", "G5_carbapenemase_family": "G5"}
COLS = {  # (symbol, method, coverage, type) for the two file formats
    "new": ("Element symbol", "Method", "% Coverage of reference", "Type"),
    "supplied": ("Gene symbol", "Method", "% Coverage of reference sequence", "Element type"),
}
def load_table(release):
    t = pd.read_csv(os.path.join(DATA, f"allele_classes_{release}.csv"), dtype=str).fillna("")
    t = t[t.release == release].drop_duplicates("allele")
    return t.set_index("allele").plan_class.to_dict()
def family_of(sym):
    m = re.match(r"^(bla[A-Za-z]+(?:-M)?)", sym)
    return m.group(1) if m else sym
def hits_to_states(df, fmt, table):
    """Return (per-hit DataFrame, per-isolate component states) for one genome's annotation table."""
    sc, mc, cc, tc = COLS[fmt]
    rows = []
    for _, r in df.iterrows():
        sym, meth, typ = str(r[sc]), str(r[mc]), str(r[tc])
        cov = pd.to_numeric(r[cc], errors="coerce")
        if typ != "AMR": continue
        if sym.startswith("ampC_"):
            comp = "G4"; state = "POSITIVE" if sym in PROM_INCLUDED else "INDETERMINATE"
            rows.append((sym, meth, cov, comp, state)); continue
        if not FAM.match(sym): continue
        complete = (meth in COMPLETE_METHODS) and (cov >= COV_MIN)
        if sym in table:
            pc = table[sym]
            if pc in CLASS_TO_COMP:
                comp = CLASS_TO_COMP[pc]; state = "POSITIVE" if complete else "INDETERMINATE"
            elif pc == "INDETERMINATE":
                comp = "G1" if sym.startswith("blaCTX-M") else "G2"; state = "INDETERMINATE"
            else:  # documented non-determinant: only an incomplete hit in an ambiguous family stays uncertain
                if complete: continue
                fam = family_of(sym); comp = "G2" if fam in FAM_G2_AMBIG else None
                if comp is None: continue
                state = "INDETERMINATE"
        else:  # family-level symbol without an allele number, or an allele missing from the table
            fam = family_of(sym)
            if fam in FAM_G1: comp, state = "G1", ("POSITIVE" if complete else "INDETERMINATE")
            elif fam in FAM_G3: comp, state = "G3", ("POSITIVE" if complete else "INDETERMINATE")
            elif fam in FAM_G5: comp, state = "G5", ("POSITIVE" if complete else "INDETERMINATE")
            elif fam in FAM_G2_AMBIG: comp, state = "G2", "INDETERMINATE"
            else: comp, state = "G2", "INDETERMINATE"
        rows.append((sym, meth, cov, comp, state))
    hits = pd.DataFrame(rows, columns=["sym", "method", "coverage", "component", "state"])
    comp_state = {}
    for c in ("G1", "G2", "G3", "G4", "G5"):
        s = set(hits[hits.component == c].state)
        comp_state[c] = "POSITIVE" if "POSITIVE" in s else ("INDETERMINATE" if "INDETERMINATE" in s else "NEGATIVE")
    return hits, comp_state
def rule_states(comp_state):
    """Cumulative rule states: Gk = components G1..Gk (G5 supplementary). POSITIVE if any component positive; else INDETERMINATE if any indeterminate; else NEGATIVE."""
    out = {}
    for k in (1, 2, 3, 4, 5):
        comps = [comp_state[f"G{i}"] for i in range(1, k + 1)]
        out[f"R{k}"] = "POSITIVE" if "POSITIVE" in comps else ("INDETERMINATE" if "INDETERMINATE" in comps else "NEGATIVE")
    return out
