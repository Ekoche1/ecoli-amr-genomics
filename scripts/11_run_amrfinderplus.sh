#!/bin/bash
# Step 11: find resistance genes in the 749 retained genomes with AMRFinderPlus.
#
# Needs AMRFinderPlus 4.2.7 with database 2026-08-07.1, with `amrfinder` on the PATH, for example:
#   conda create -n amr -c conda-forge -c bioconda ncbi-amrfinderplus=4.2.7
#   conda activate amr
#   amrfinder_update --database_version 2026-08-07.1   (check with: amrfinder --database_version)
# Runs on Linux, macOS or Windows (WSL). The command used for every genome is:
#   amrfinder -n <genome.fna> --organism Escherichia --threads <T> -o <genome.tsv>
#
# Reads:  data/genome_qc_762.csv (to list the retained genomes), genomes/
# Writes: amrfinder_output/<id>.tsv and <id>.log, timing.tsv, RUN_MANIFEST.txt   (resumable; not stored in the repository)
# Usage (from the project folder): bash scripts/11_run_amrfinderplus.sh [parallel_jobs] [threads_per_job]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GEN="$ROOT/genomes"; OUT="$ROOT/amrfinder_output"; P="${1:-2}"; T="${2:-2}"
mkdir -p "$OUT"
IDS="$OUT/retained_ids.txt"
# retained genomes = every genome in the QC table whose tier is not CONFIRMED
python3 - "$ROOT/data/genome_qc_762.csv" > "$IDS" <<'EOF'
import csv, sys
for row in csv.DictReader(open(sys.argv[1], newline="")):
    if row["tier"] != "CONFIRMED":
        print(row["id"])
EOF
echo "retained genomes to annotate: $(wc -l < "$IDS")"
export GEN OUT T
if [ ! -s "$OUT/RUN_MANIFEST.txt" ]; then
  {
    echo "run_started: $(date -Is)"
    echo "command_template: amrfinder -n <genome.fna> --organism Escherichia --threads $T -o <out.tsv>"
    amrfinder --version 2>&1 | grep -v "screen size" | sed 's/^/amrfinder_version: /'
    amrfinder --database_version 2>&1 | grep -E "Software version|Database version" | sed 's/^/  /'
    echo "parallel_jobs: $P ; threads_per_job: $T"
  } > "$OUT/RUN_MANIFEST.txt"
fi
one() {
  id="$1"
  if [ -s "$OUT/$id.tsv" ] && grep -q "took" "$OUT/$id.log" 2>/dev/null; then echo "skip $id"; return; fi
  s=$(date +%s)
  amrfinder -n "$GEN/$id.fna" --organism Escherichia --threads "$T" -o "$OUT/$id.tmp" > "$OUT/$id.log" 2>&1
  rc=$?
  e=$(date +%s)
  if [ $rc -eq 0 ] && [ -s "$OUT/$id.tmp" ]; then mv "$OUT/$id.tmp" "$OUT/$id.tsv"; else echo "FAILED rc=$rc" >> "$OUT/$id.log"; fi
  printf "%s\t%s\t%s\n" "$id" "$rc" "$((e-s))" >> "$OUT/timing.tsv"
  echo "done $id rc=$rc $((e-s))s"
}
export -f one
tr -d '\r' < "$IDS" | grep -v '^$' | xargs -P "$P" -I{} bash -c 'one {}'
echo "RUN_FINISHED $(date -Is)" >> "$OUT/RUN_MANIFEST.txt"
echo ALL_DONE
