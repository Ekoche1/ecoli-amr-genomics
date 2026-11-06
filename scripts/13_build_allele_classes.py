"""Step 13: sort every beta-lactamase allele in the AMRFinderPlus gene catalogue into the gene groups used here.
Classes: G1_CTXM_ESBL, G2_other_ESBL, G3_acquired_AmpC_family, G5_carbapenemase_family, NONE_documented,
NONE_intrinsic_excluded, INDETERMINATE. Each class comes from the catalogue's own gene family and product name
(for example, whether the product name says "extended-spectrum"). No lab results are used.
Reads:  ReferenceGeneCatalog.txt and version.txt in the AMRFinderPlus database folder (database 2026-08-07.1)
Writes: data/allele_classes_2026-08-07.1.csv
Usage:  python scripts/13_build_allele_classes.py <AMRFinderPlus database folder>
        (the folder is shown by: amrfinder --database_version)"""
import os, re, sys
import pandas as pd
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
RELEASE = "2026-08-07.1"
DB = sys.argv[1]
version = open(os.path.join(DB, "version.txt")).read().strip()
if version != RELEASE: sys.exit(f"this database is version {version}; the analysis used {RELEASE}")
pd.set_option("display.width", 250)
ACQ_AMPC = {"blaCMY", "blaDHA", "blaACT", "blaFOX", "blaMOX", "blaMIR", "blaACC", "blaCFE", "blaLAT", "blaMOR"}
CARBA = {"blaNDM", "blaKPC", "blaVIM", "blaIMP"}
ES = re.compile(r"extended[- ]spectrum", re.I)
OXA163 = {"blaOXA-163"}          # documented exception: hydrolyses broad-spectrum cephalosporins (Poirel 2012)
def classify(r):
    fam, name, sub, allele = r.gene_family, (r.product_name or ""), (r.subclass or ""), r.allele
    if fam == "blaEC": return "NONE_intrinsic_excluded", "intrinsic chromosomal E. coli ampC (blaEC): never counted"
    if fam == "blaCTX-M":
        return ("G1_CTXM_ESBL", "catalogue product name states extended-spectrum") if ES.search(name) else ("INDETERMINATE", "CTX-M allele whose catalogue product name does not state extended-spectrum")
    if allele in OXA163: return "INDETERMINATE", "OXA-48-family exception documented to hydrolyse broad-spectrum cephalosporins (Poirel 2012)"
    if fam in ACQ_AMPC: return "G3_acquired_AmpC_family", "acquired AmpC-family gene by family; plasmid location NOT confirmed"
    if fam in CARBA: return "G5_carbapenemase_family", "carbapenemase family (supplementary rule G5)"
    if "promoter" in name.lower(): return "NONE_documented", "promoter-region entry, not a coding determinant"
    if fam == "blaOXA":
        if ES.search(name): return "G2_other_ESBL", "catalogue product name states extended-spectrum class D"
        if "carbapenem-hydrolyzing" in name.lower() or "carbapenem-hydrolysing" in name.lower(): return "NONE_documented", "catalogue: carbapenem-hydrolysing OXA (OXA-48 family etc.); not counted as cephalosporinase"
        if "CEPHALOSPORIN" in sub and not ES.search(name): return "NONE_documented", "descriptor-subclass conflict: not extended-spectrum by product name although subclass says CEPHALOSPORIN (S5c moves these)"
        return "NONE_documented", "OXA without extended-spectrum descriptor"
    if re.search(r"class A", name):
        if ES.search(name): return "G2_other_ESBL", "catalogue product name states extended-spectrum class A"
        if re.search(r"broad-spectrum", name, re.I): return "NONE_documented", "catalogue product name states broad-spectrum (non-ESBL)"
        if "carbapenem-hydrolyzing" in name.lower(): return "INDETERMINATE", "class A carbapenem-hydrolysing (e.g. GES types); ceftazidime effect not stated in descriptor"
        return "INDETERMINATE", "class A beta-lactamase with no functional descriptor in the catalogue"
    return "NONE_documented", "not a class A/C/D determinant covered by the rules"
df = pd.read_csv(os.path.join(DB, "ReferenceGeneCatalog.txt"), sep="\t", dtype=str).fillna("")
df = df[(df["class"].str.contains("BETA-LACTAM") | df.subclass.str.contains("CEPHALOSPORIN|CARBAPENEM|BETA-LACTAM")) & df.gene_family.str.startswith("bla")]
cl = df.apply(classify, axis=1, result_type="expand"); df = df.assign(plan_class=cl[0], rationale=cl[1], release=RELEASE)
tab = df[["release", "allele", "gene_family", "product_name", "subclass", "plan_class", "rationale"]]
print(f"=== database {RELEASE}: {len(df)} beta-lactam-family alleles ===")
print(df.plan_class.value_counts().to_string())
tab.to_csv(os.path.join(DATA, f"allele_classes_{RELEASE}.csv"), index=False)
