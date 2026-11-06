# Data

| File | What it is | Where it came from |
|---|---|---|
| `meta_data_tidy.txt` | Laboratory results for 736 isolates: disk diffusion results (S/I/R) and MICs for several antibiotics. The analysis uses only `mic_CAZ`, the ceftazidime MIC in mg/L. | Liverpool dataset, unchanged |
| `liverpool_dataset_README.md` | The dataset authors' own description of their data | Liverpool dataset, unchanged |
| `genome_qc_762.csv` | Quality-check measurements and the decision for all 762 genomes. `tier` is RETAIN, SUSPECTED (kept, but flagged) or CONFIRMED (excluded). `problem` and `evidence` explain each exclusion. | Scripts 03–10 |
| `near_identical_pairs.csv` | Pairs of genomes that share at least half of their sampled DNA k-mers (`jaccard` ≥ 0.5), with an estimate of how identical they are (`mash_ani`). Pairs at 0.99 or above are treated as near-identical. | Script 08 |
| `amrfinderplus_hits_749.csv` | Every AMRFinderPlus hit in the 749 genomes that passed the quality checks, one row per hit, with AMRFinderPlus's own columns. 35 genomes had no hits. | Scripts 11–12 |
| `allele_classes_2026-08-07.1.csv` | Every beta-lactamase variant in the AMRFinderPlus catalogue, with its gene group (`plan_class`) and the reason for it (`rationale`) | Script 13 |
| `gene_hits_749.csv` | The hits that count towards the gene groups: gene (`sym`), method, coverage, group (`component`) and call (`state`) | Script 14 |
| `gene_calls_749.csv` | **The main gene calls.** One row per genome: positive, negative or indeterminate for each gene group (`C_G1`–`C_G5`) and each cumulative step (`R1`–`R5`). | Script 14 |
| `gene_calls_749_S10a.csv`, `gene_calls_749_S10b.csv` | Two alternative versions of the gene calls, used as a check. Some hits are named only to gene family, for example "blaCTX-M" with no variant number. The main calls count a complete family-level hit as positive. S10a counts family-level CTX-M hits as indeterminate. S10b does the same for family-level CTX-M and VIM hits. Between them they change the calls for 3 genomes (E71, H12 and I76). | Script 15 |

**Not included:** the 762 genomes (about 4 GB). Script `01_download_genomes.py` downloads them from the Liverpool dataset.

**Licence:** the Liverpool files are CC BY 4.0 (Gerada, Harper and Hope, 2025). The tables I produced are also shared under CC BY 4.0.
