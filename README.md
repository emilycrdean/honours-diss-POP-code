# Characterising rare variants in pelvic organ prolapse

Analysis code for the BSc (Hons) Genetics dissertation *Characterising rare variants in pelvic organ prolapse* (GENE490, University of Otago, 2026).

The project asked two questions:

1. **Are heritable connective tissue disorders associated with pelvic organ prolapse (POP)?** Ehlers–Danlos syndrome and Marfan syndrome (primary), with hypermobility syndrome and osteogenesis imperfecta (exploratory), in the All of Us Research Program. Marfan syndrome was also analysed in UK Biobank.
2. **Are rare pathogenic or likely pathogenic variants in connective tissue genes enriched in early-onset POP?** Women diagnosed with POP at age 18–35 without a connective tissue disorder diagnosis were compared with three comparison groups in All of Us whole-genome sequence data. An exploratory polygenic score analysis is also included.

**No data are included in this repository.** All analyses were run inside the All of Us Researcher Workbench or on an approved UK Biobank extract, and the scripts only work with authorised access.

## Repository structure

| Folder | Contents |
|---|---|
| `01_ukb_marfan/` | UK Biobank Marfan syndrome analysis (Aim 1) |
| `02_aou_ctd/` | All of Us connective tissue disorder analyses (Aim 1) |
| `03_aou_rare_variants/` | All of Us rare variant analysis (Aim 2) |
| `04_prs/` | Exploratory polygenic score analysis (Supplementary Methods 1) |


### 01_ukb_marfan

- `marfan_analysis.Rmd` is the complete analysis as an R Markdown document. It covers:
  - reading and repairing the extract (2,491 control lines with a misplaced field);
  - cleaning, then 1:4 propensity score matching on age and BMI within birth cohort;
  - conditional logistic regression and Fisher's exact test;
  - sensitivity analyses: adjustment for age at first birth, 1:10 matching, full-cohort adjusted logistic regression, and 500 random draws of controls within the caliper;
  - exploratory Kaplan–Meier and Cox models of age at POP, and a parity interaction.
- `ukb_extra_summaries.R` gives supporting summaries: cases vs matched controls vs the whole control pool, and POP by number of live births and birth cohort. Run it after the matching chunks of `marfan_analysis.Rmd`.

### 02_aou_ctd

These scripts are pasted as cells into Workbench R notebooks, one notebook per cohort, in the numbered order. They go below the cells the Dataset Builder generates.

Cohorts were built in the All of Us Cohort Builder using the concept sets listed in Supplementary Table S2 of the thesis. The Dataset Builder generates the code that loads each cohort, and that code is specific to each workspace, so it is not included

- `01_build_analysis_dataset.R`: one row per participant. This step:
  - classifies POP and menopause codes;
  - takes the median BMI and adds deprivation index;
  - saves each cohort's dataset and the controls' smoking status;
  - combines cases with the control pool and derives race/ethnicity.
- `02_run_all.R`: removes participants with hypermobility syndrome or osteogenesis imperfecta from the control pool, then:
  - runs the primary analysis (1:4 matching on age and BMI, caliper 0.2 SD, conditional logistic regression);
  - runs the main sensitivity analyses (age-only matching, matched and adjusted, crude Fisher's exact test).
- `03_additional_sensitivity.R`: the additional sensitivity analyses in Supplementary Table Y:
  - full-cohort logistic regression (adjusted, then adding the menopause code, then smoking);
  - Fisher's exact test in the matched sample;
  - exact matching on race/ethnicity;
  - additional matching on deprivation.
- `04_describe.R`: descriptive tables for matched cases and controls, plus the full-cohort covariate model.
- `05_exploratory_ctds.R`: the exploratory exposures (osteogenesis imperfecta, hypermobility syndrome, and cutis laxa, which was considered but not analysed).
- `06_marfan_resampling.R`: Marfan syndrome only. Refits the model on 500 random sets of controls drawn from within the caliper, and summarises the distribution of the odds ratio.

### 03_aou_rare_variants

This is a Python notebook pipeline. Each `# %%` block is one cell.

- `vpa_gene_list.py` builds `vpa_gene_list.csv`, the gene panels. The 167-gene revised panel is primary.
- `vpa_mechanisms.csv` assigns genes to the nine mechanism groups.
- `vpa_pipeline.py` contains the main cells:
  - CELL 0–1: locate workspace resources and set the configuration.
  - CELL 2: extract Variant Annotation Table rows for the panel genes.
  - CELL 3: define qualifying variants by ClinVar tier, with a per-gene audit of every filter step.
  - CELL 4–5: load the cohorts and identify carriers.
  - CELL 6: carrier comparisons with Fisher's exact test.
  - CELL 6b: mechanism-level results, Benjamini–Hochberg corrected.
  - CELL 7: sensitivity analyses.
  - CELL 8: ancestry-adjusted regression.
- `vpa_breakdown_cell.py` (CELL 9) gives the per-panel, per-mechanism and per-gene breakdown.
- `vpa_table_cell.py` (CELL 10) produces the thesis results tables.
- `vpa_followup_cells.py` contains three follow-up cells:
  - CELL 11: which variants account for most comparison-group carriers.
  - CELL 12: mechanism results against all four comparison groups, with and without variants that have benign submissions.
  - CELL 13: gene-level results.

### 04_prs

- `prs_pipeline.py`:
  - PRS-1 to PRS-3: the 29-SNP score built from the Pujol-Gualdo et al. (2022) lead variants, and its association with early-onset and later-onset POP.
  - PRS-G1 and PRS-G2: download PGS002288 and run a single-chromosome pilot.
- `prs_genomewide.py` (PRS-G3 and PRS-G4) builds the genome-wide PGS002288 score chromosome by chromosome in Hail, then runs the association models (age plus PC1–10).
- `prs_gw_diagnose.py` (PRS-G5) is a quality-control check. It compares the score's weights with the published odds ratios at the lead SNPs. The weights were oriented opposite to the published effects at every lead SNP, so the score sign was reversed before analysis.
- `pujol_gualdo_30_lead_variants.csv` holds the published lead variants (public summary statistics).

### 05_figures_tables

- `thesis_theme.R`: shared colours, theme and plotting helpers. Run it first.
- `aim1_figures.R`, `flow_diagrams.R`, `pva_results_figures.R` and `aou_or_bars.R`: thesis figures, drawn from aggregate results only.
- `make_*.js`: build the Word tables with the `docx` npm package (`node make_aim1_tables.js` and so on). They write to `thesis_figures/`.


Workspace-specific identifiers have been replaced with placeholders. Before running, replace these with your own workspace's values:

- `YOUR-WORKSPACE-ID` (the workspace Google project, used in bucket names);
- `YOUR-CDR-PROJECT`.

## Data and Statistics Dissemination Policy

In line with the All of Us Data and Statistics Dissemination Policy, the scripts replace any participant count from 1 to 20 with "<20" in output intended for export. They also withhold estimates from which such counts could be derived.

- Individual-level data and scores never leave the Workbench.
- Intermediate files are written to `/tmp` or to the workspace bucket, never to the notebook folder.

## Software

- **R (≥ 4.3):** tidyverse, dplyr, lubridate, MatchIt, survival, survminer, broom, tableone.
- **Python 3:** pandas, numpy, scipy, statsmodels, hail, google-cloud-bigquery.
- **Node.js:** docx (tables only).
