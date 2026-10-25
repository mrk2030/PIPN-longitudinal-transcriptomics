# PIPN longitudinal transcriptomics

Analysis code and selected reproducibility materials associated with the manuscript:

**Temporal Transcriptomic Dynamics Reveal Molecular Attenuation and Directional Pathway Reorganization in Paclitaxel-Induced Peripheral Neuropathy**

## Overview

This repository contains the R scripts used to reanalyze public dorsal root ganglion (DRG) transcriptomic datasets in paclitaxel-induced peripheral neuropathy (PIPN). The study examines treatment-associated transcriptional responses across behaviorally defined stages of neuropathy development and recovery, classifies longitudinal gene trajectories, evaluates Hallmark pathway dynamics, and assesses selective pathway-level reproducibility in an independent DRG dataset.

## Release scope and reproducibility status

The current public GitHub repository contains the initial R analysis workflow and selected manuscript-supporting results. Additional analyses supporting the revised manuscript, including the targeted day-23-versus-day-11 treatment interaction, DRG09 sample-exclusion sensitivity analysis, Hallmark implementation-sensitivity assessment, and expanded GSE113941 reproducibility verification, have been numerically verified in a separate pre-release environment.

The corresponding additional scripts, analytical outputs, and supporting documentation are being prepared for the updated Zenodo version 1.1.0 release. Until that release is finalized, the public GitHub repository should not be considered a complete reproduction package for every analysis reported in the revised manuscript.

## Public datasets

- **GSE185084**: primary longitudinal mouse DRG RNA-sequencing dataset.
- **GSE113941**: independent mouse DRG dataset used for pathway-level reproducibility analysis.
- **GSE249643**: candidate external dataset screened for quality control. The three deposited paclitaxel DRG CEL files were found to be exactly identical at the raw-intensity level and were therefore excluded from formal external evaluation.

An exploratory analysis of GSE230454/GSE230459 was performed during project development but did not contribute to the final manuscript and is intentionally not included in this release.

## Repository structure

```text
PIPN-longitudinal-transcriptomics/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── environment/
│   └── package_versions.txt
├── metadata/
│   └── GSE185084_DRG_sample_manifest_TEMPLATE.csv
├── results/
│   └── selected manuscript-supporting CSV outputs
├── scripts/
│   ├── 01_setup_project.R
│   ├── 02_download_GSE185084_metadata.R
│   ├── 03_download_GSE185084_counts.R
│   ├── 04_build_GSE185084_DRG_manifest.R
│   ├── 05_GSE185084_QC_and_DESeq2_setup.R
│   ├── 06_GSE185084_differential_expression.R
│   ├── 07_GSE185084_pathway_analysis.R
│   ├── 08_GSE185084_trajectory_enrichment.R
│   ├── 10_GSE249643_microarray_QC.R
│   └── 11_GSE113941_RNAseq_validation.R
└── supplementary/
    └── Supplementary_Tables_FINAL.xlsx
```

## Software environment

The analysis was performed in **R 4.6.1** using RStudio on Windows. Principal package versions are listed in `environment/package_versions.txt`.

Key packages include DESeq2, fgsea, msigdbr, GEOquery, limma, clusterProfiler, org.Mm.eg.db, oligo, ggplot2, and ggrepel.

## Analysis workflow

1. `01_setup_project.R` creates the expected directory structure.
2. `02_download_GSE185084_metadata.R` downloads and inspects GSE185084 metadata.
3. `03_download_GSE185084_counts.R` downloads GSE185084 supplementary files.
4. `04_build_GSE185084_DRG_manifest.R` extracts DRG metadata for manual curation of the final sample manifest.
5. `05_GSE185084_QC_and_DESeq2_setup.R` performs count filtering, DESeq2 setup, variance-stabilizing transformation, PCA, sample-distance analysis, and QC.
6. `06_GSE185084_differential_expression.R` performs time-specific differential expression, treatment-by-time testing, and longitudinal trajectory classification.
7. `07_GSE185084_pathway_analysis.R` performs Hallmark preranked gene-set enrichment analysis across days 4, 7, 11, and 23.
8. `08_GSE185084_trajectory_enrichment.R` performs Gene Ontology analysis of resolving and unresolved trajectory classes and summarizes unresolved driver genes.
9. `10_GSE249643_microarray_QC.R` checks raw CEL-file independence in GSE249643 and documents the duplicate-array exclusion.
10. `11_GSE113941_RNAseq_validation.R` performs the independent GSE113941 gene- and pathway-level reproducibility analysis.

## Important note on the GSE185084 sample manifest

The GEO metadata required manual curation to resolve sample labels and construct the final analysis manifest. The repository includes a manifest template but does not invent or reconstruct metadata that were not embedded directly in the archived scripts. The final manuscript-supporting outputs are included in `results/` for verification.

For exact reproduction from raw data, users should curate the DRG manifest from the GEO metadata using the same treatment and collection-day assignments reported in the manuscript before running Scripts 05 onward.

## Principal analysis definitions

- Low-count filter: at least 10 raw counts in at least 4 samples.
- Differential-expression model: `~ day * treatment`.
- Statistical significance: Benjamini-Hochberg adjusted P < 0.05 unless otherwise specified.
- Strongly resolving: day-11 significant and absolute day-23 effect <50% of the day-11 effect.
- Residual: day-11 significant, same direction at day 23, at least 50% of the day-11 effect retained, but day-23 FDR not significant.
- Persistent: day-11 significant, same direction at day 23, at least 50% retained, and day-23 FDR <0.05.
- Hallmark GSEA: signed differential-expression statistic, gene-set size 10-500, random seed 123.

## Expected headline results

The scripts and included result tables support the following manuscript-level findings:

- 0, 4, 62, and 16 FDR-significant genes at days 4, 7, 11, and 23, respectively.
- Among 62 day-11-responsive genes: 50 strongly resolving, 8 residual, and 4 persistent.
- 19, 13, 24, and 6 significant Hallmark pathways at days 4, 7, 11, and 23.
- Three day-11 Hallmark programs remained significant at day 23 with reversed enrichment direction: Myogenesis, Epithelial-Mesenchymal Transition, and Protein Secretion.
- GSE113941 showed selective pathway-level reproducibility rather than global replication.

## Data availability

Primary transcriptomic data are publicly available through NCBI GEO under accession numbers GSE185084, GSE113941, and GSE249643.

Selected processed outputs that support the manuscript are included in this repository. The supplementary tables supplied with the manuscript are included under `supplementary/`.

## Reproducibility notes

A random seed of 123 was used where stochastic procedures were involved. The scripts use relative project paths rather than user-specific absolute paths.

Because the work was developed interactively during exploratory analysis, this public release removes malformed backups and exploratory analyses that did not contribute to the final manuscript. Analytical definitions and statistical procedures supporting the reported results are unchanged.

## Citation

The previously archived analysis code is referenced as:

Khan, M. R., & Alhosaini, K. A. (2026). PIPN longitudinal transcriptomics analysis code (Version 1.0.1) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23390467

An updated reproducibility package (version 1.1.0) is being finalized. Its release information and citation will be added after the public Zenodo record is verified.

## Contact

Mohammad Rashid Khan  
Department of Pharmacology and Toxicology  
College of Pharmacy, King Saud University  
Riyadh, Saudi Arabia
