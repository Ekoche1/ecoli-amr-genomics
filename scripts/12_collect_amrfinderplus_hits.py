"""Step 12: check the AMRFinderPlus outputs and combine them into one table.
Checks that every retained genome has an output with the expected columns, and that every log shows the same
software and database version and finished without errors.
Reads:  amrfinder_output/, data/genome_qc_762.csv
Writes: data/amrfinderplus_hits_749.csv   (one row per AMRFinderPlus hit; a genome with no hits has no rows)
Usage:  python scripts/12_collect_amrfinderplus_hits.py"""
import os, re, collections
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); AF = os.path.join(ROOT, "amrfinder_output"); DATA = os.path.join(ROOT, "data")
EXPECTED_COLS = ['Protein id', 'Contig id', 'Start', 'Stop', 'Strand', 'Element symbol', 'Element name', 'Scope', 'Type', 'Subtype', 'Class', 'Subclass', 'Method', 'Target length', 'Reference sequence length', '% Coverage of reference', '% Identity to reference', 'Alignment length', 'Closest reference accession', 'Closest reference name', 'HMM accession', 'HMM description']
d = pd.read_csv(os.path.join(DATA, "genome_qc_762.csv"), dtype={"id": str}); retained = set(d[d.tier != "CONFIRMED"].id)
ids = sorted(f[:-4] for f in os.listdir(AF) if f.endswith(".tsv") and f != "timing.tsv")
print("=== 1. ID ACCOUNTING ===")
print(f"retained genomes expected: {len(retained)} | outputs found: {len(ids)}")
print("outputs for genomes that are not retained:", sorted(set(ids) - retained))
print("retained genomes without an output:", sorted(retained - set(ids)))
print("all annotated ids have a genome file:", all(os.path.exists(os.path.join(ROOT, "genomes", i + ".fna")) for i in ids))
print("=== 2. FILE STRUCTURE ===")
rows = []; bad_cols = []; header_only = []; dup_rows = 0
for i in ids:
    df = pd.read_csv(os.path.join(AF, i + ".tsv"), sep="\t", dtype=str)
    if list(df.columns) != EXPECTED_COLS: bad_cols.append(i)
    if len(df) == 0: header_only.append(i)
    if df.duplicated().any(): dup_rows += int(df.duplicated().sum())
    df.insert(0, "genome", i); rows.append(df)
allhits = pd.concat(rows, ignore_index=True)
allhits = allhits[allhits.genome.isin(retained)]
allhits.to_csv(os.path.join(DATA, "amrfinderplus_hits_749.csv"), index=False)
print("files with unexpected columns:", bad_cols, "| hit rows (retained genomes, all types):", len(allhits), "| duplicate rows:", dup_rows)
print("genomes with no hits at all:", len(header_only), header_only[:10])
print("Type values:", dict(collections.Counter(allhits.Type)), "| Scope values:", dict(collections.Counter(allhits.Scope)))
print("=== 3. SOFTWARE / DATABASE VERSION IN EVERY LOG ===")
ok = collections.Counter(); problems = []
for i in ids:
    t = open(os.path.join(AF, i + ".log"), encoding="utf-8", errors="replace").read()
    sw = re.search(r"Software version: (\S+)", t); db = re.search(r"Database version: (\S+)", t)
    ok[(sw.group(1) if sw else None, db.group(1) if db else None)] += 1
    if "took" not in t or re.search(r"\bERROR\b|\bWARNING\b|Error:", t): problems.append(i)
    if "--organism Escherichia" not in t and "-O Escherichia" not in t: problems.append(i + ":organism-option-missing-in-log")
print("(software, database) over", sum(ok.values()), "logs:", dict(ok))
print("logs with errors/warnings, no completion line, or missing organism option:", problems[:10], len(problems))
t = pd.read_csv(os.path.join(AF, "timing.tsv"), sep="\t", header=None, names=["id", "rc", "seconds"])
print(f"timing.tsv: rows {len(t)} | non-zero exit codes {int((t.rc != 0).sum())} | duplicate ids {int(t.id.duplicated().sum())}")
print("--- run manifest"); print(open(os.path.join(AF, "RUN_MANIFEST.txt")).read().strip()[:700])
