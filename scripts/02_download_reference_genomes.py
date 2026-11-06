"""Step 2: download the 22 reference genomes used by the quality checks, from NCBI RefSeq.

They cover E. coli and its close relatives (Shigella, E. fergusonii) plus the other genera most likely to
contaminate a clinical sample. The list of accessions is in qclib.py (GROUPS).

Writes: reference_genomes/<accession>.fna   (about 100 MB; not stored in the repository)
Usage:  python scripts/02_download_reference_genomes.py
"""
import os, sys, time, urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qclib as Q

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text"
os.makedirs(Q.REFS, exist_ok=True)
accessions = [a for accs in Q.GROUPS.values() for a in accs]
print("reference genomes:", len(accessions), flush=True)
for acc in accessions:
    dest = os.path.join(Q.REFS, acc + ".fna")
    if os.path.exists(dest):
        print("already present", acc); continue
    for attempt in range(10):
        try:
            data = urllib.request.urlopen(EFETCH.format(acc), timeout=300).read()
            if not data.startswith((">" + acc).encode()): raise IOError("unexpected response")
            open(dest + ".part", "wb").write(data); os.replace(dest + ".part", dest)
            print("downloaded", acc, round(len(data) / 1e6, 1), "MB", flush=True)
            break
        except Exception as ex:
            print("retry", acc, ex, flush=True); time.sleep(5 + 5 * attempt)
    else:
        print("FAILED", acc, flush=True)
    time.sleep(0.5)   # stay well inside NCBI's limit of 3 requests per second
print("finished:", sum(os.path.exists(os.path.join(Q.REFS, a + ".fna")) for a in accessions), "of", len(accessions))
