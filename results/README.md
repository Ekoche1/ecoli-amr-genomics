# Results

All files here are written by scripts 16–18. **Start with `results_summary.csv`.**

**Terms used in the tables**

- **G1–G5** are the gene groups from the main README. Some tables call them **R1–R5** (the rule for G1 is R1, and so on).
- **TP** (true positive): resistant in the lab, flagged by the genes. **FN** (false negative): resistant, not flagged. **FP** (false positive): susceptible, flagged. **TN** (true negative): susceptible, not flagged.
- **Analysis sets:**
  - *own evaluable set*: the isolates where that gene group gave a definite answer.
  - *common set C(R1,R4)*: the 623 isolates where both G1 and G4 gave a definite answer. This is the fair G1 vs G4 comparison.
  - *common set, evaluable for all of R1–R4*: the 622 isolates where G1, G2, G3 and G4 all gave a definite answer.
  - *common set C(R4,R5)*: the 622 isolates where G4 and G5 both gave a definite answer.
- **Confidence intervals:** the main ones are *cluster bootstrap* intervals, which treat near-identical isolates as one group so they are not over-counted. *Wilson* intervals are a simpler version that ignores this, included for comparison.
- **PPV / NPV** (positive / negative predictive value) depend on how common resistance is in this collection. Because the collection is not a random sample, they should not be applied to a routine laboratory population.

**The tables**

| File | What it shows |
|---|---|
| `results_summary.csv` | **The headline table.** Sensitivity and specificity for each gene group, with 95% confidence intervals, in each analysis set. It also shows how many indeterminate calls were set aside. |
| `results_all_sets_long.csv` | The same results in full: every count and every confidence interval, including PPV and NPV. |
| `paired_contrasts.csv` | Each gene group compared directly with the next one on the same isolates. `d_sens_pp` and `d_spec_pp` are the changes in sensitivity and specificity, in percentage points, with confidence intervals. `resistant_newly_detected` and `susceptible_newly_flagged` count the isolates gained and the new false positives. `positives_lost` is 0 throughout, which means adding genes never removed a positive call. The exact McNemar p-values are unadjusted. |
| `step_additions_by_determinant.csv` | The isolates newly flagged at each step, by gene or mutation, with how many were resistant (R), intermediate (I), susceptible (S) or had no ceftazidime MIC. |
| `indeterminate_bounds.csv` | Sensitivity and specificity if all indeterminate calls had gone one way or the other: the best and worst case. |
| `reconciliation_state_by_caz_category.csv` | Every one of the 749 genomes counted, by gene call (positive, indeterminate, negative) and lab category (R, I, S, no MIC, no lab record), so that no isolate goes missing. |
| `linked_isolates_749.csv` | One row per genome: the raw table that all the results above are calculated from. Its columns are listed below. |
| `analysis_log.txt` | The printed record of the analysis run, including the 7 built-in checks (all pass). |
| `analysis_run_record.json` | The random seed, the check results and a digital fingerprint (SHA-256) of each input file, so a rerun can confirm it used the same data. |

**Columns in `linked_isolates_749.csv`**

| Column | Meaning |
|---|---|
| `genome_id` | The isolate's ID in the Liverpool dataset |
| `in_pheno_file` | Whether the isolate has a laboratory record |
| `caz_mic_raw`, `caz_mic_numeric` | Ceftazidime MIC (mg/L), as supplied and as a number |
| `caz_category` | R (above 4), S (1 or below), I (2–4), missing, or no_record |
| `primary_analysable` | True if the isolate is clearly R or S |
| `C_G1` … `C_G5` | The call for each gene group on its own |
| `R1` … `R5` | The call for each cumulative step (G1, G1+G2, … G1–G5) |
| `S10a_…`, `S10b_…` | The same calls under the two alternative versions (see `data/README.md`) |
| `qc_tier` | RETAIN, or SUSPECTED (kept, but flagged in the quality checks) |
| `clone_discordant_flag` | True for 5 isolates whose lab results disagree with near-identical genomes on 3 or more antibiotics. This is for information only: they were kept. |
| `unit_id`, `unit_size` | The group of near-identical genomes the isolate belongs to, used for the confidence intervals |
