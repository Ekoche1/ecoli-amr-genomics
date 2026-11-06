# Scripts

The scripts are numbered in the order they run. Run them from the main project folder, for example: `python scripts/16_link_and_analyse.py`

**Quick rerun:** run scripts 16 to 19 only. They use the saved tables in `data/` and take about a minute.

**Full rerun:** run all of them, in order.

| Stage | Script | What it does | Saves to |
|---|---|---|---|
| **Get the data** | `01_download_genomes.py` | Downloads the 762 genomes from the Liverpool dataset and checks every file | `genomes/` |
| | `02_download_reference_genomes.py` | Downloads 22 reference genomes from NCBI, used for the species checks | `reference_genomes/` |
| **Quality checks** | `03_measure_genomes.py` | Measures each genome: size, fragmentation, GC content, species match and 16S genes | `work/` |
| | `04_build_qc_table.py` | Puts the measurements into one table | `work/` |
| | `05_low_coverage.py` | Measures low-coverage sequence, a possible sign of a second organism | `work/` |
| | `06_low_coverage_species.py` | Checks what that low-coverage sequence resembles | `work/` |
| | `07_find_16s.py` | Identifies the species from 16S rRNA genes, as a second check | `work/` |
| | `08_near_identical_pairs.py` | Finds near-identical genomes | `data/near_identical_pairs.csv` |
| | `09_assign_qc_tiers.py` | Decides which genomes to exclude, using genome measurements only | `work/` |
| | `10_final_qc_table.py` | Saves the final quality-check table | `data/genome_qc_762.csv` |
| **Gene search** | `11_run_amrfinderplus.sh` | Runs AMRFinderPlus on the 749 genomes that passed | `amrfinder_output/` |
| | `12_collect_amrfinderplus_hits.py` | Checks the AMRFinderPlus outputs and combines them into one table | `data/amrfinderplus_hits_749.csv` |
| **Gene classification** | `13_build_allele_classes.py` | Sorts every beta-lactamase variant in the AMRFinderPlus catalogue into the gene groups | `data/allele_classes_2026-08-07.1.csv` |
| | `14_call_genes.py` | Calls each genome positive, negative or indeterminate for G1–G5 | `data/gene_calls_749.csv`, `data/gene_hits_749.csv` |
| | `15_alternative_gene_calls.py` | Makes two alternative versions of the gene calls, used as a check | `data/gene_calls_749_S10a.csv`, `data/gene_calls_749_S10b.csv` |
| **Analysis** | `16_link_and_analyse.py` | Links the gene calls to the ceftazidime MICs and calculates sensitivity, specificity and confidence intervals | `results/` |
| | `17_genes_behind_each_step.py` | Shows which genes make the difference at each step | `results/step_additions_by_determinant.csv` |
| | `18_summary_tables.py` | Makes the headline table and the best and worst case for indeterminate calls | `results/results_summary.csv`, `results/indeterminate_bounds.csv` |
| **Figure** | `19_make_figure.py` | Draws the main figure: sensitivity and specificity at each step | `figures/sensitivity_specificity_by_step.png` |

Two helper files are used by the scripts above: `qclib.py` (shared quality-check functions) and `gene_rules.py` (the gene-calling rules).

## What you need

- **Python 3.12** with the packages in `requirements.txt`: `pip install -r requirements.txt`
- **For the full rerun only:** AMRFinderPlus 4.2.7 with database version 2026-08-07.1, on Linux, macOS or Windows with WSL:
  ```bash
  conda create -n amr -c conda-forge -c bioconda ncbi-amrfinderplus=4.2.7
  conda activate amr
  amrfinder_update --database_version 2026-08-07.1
  ```
  Script 13 needs the database folder, which `amrfinder --database_version` shows:
  `python scripts/13_build_allele_classes.py <database folder>`

## Good to know

- **The `genomes/`, `reference_genomes/`, `amrfinder_output/` and `work/` folders are created when you run the scripts.** They are not stored in this repository because they are large.
- **A full rerun rewrites some files in `data/`.** Their content will be the same as the saved versions, but the rows may be in a different order.
- **Script 11 (AMRFinderPlus):** the full 749-genome run was done once. After tidying, script 11 was tested on 3 genomes and matched the original run exactly.
- **Rough times on a standard laptop:** 2 hours to download, 2 hours for the quality checks, 6 hours for AMRFinderPlus, and a few minutes for everything else.
