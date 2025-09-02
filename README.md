# Data Availability – Microbiome Meta-analysis

This folder contains all files associated with the microbiome meta-analysis study.  
It includes ASV and taxonomy tables, metadata, alpha/beta diversity outputs, PICRUSt2 predictions, and all analysis scripts.

---

## 1. ASV and Taxonomy Tables
- `seqtab_feature_table.tsv`, `seqtab.nochim_feature_table.tsv`, `seqtab.nochim.matched_feature_table.tsv`  
  → Feature tables (ASV abundances across samples).  
- `taxa.tsv`, `taxa_new.tsv`  
  → Taxonomic assignments for each ASV.

---

## 2. Metadata
- `16S_metadata.txt`, `16S_metadata.xlsx`  
- `ITS_metadata.txt`, `ITS_metadata.xlsx`  
- `all_projects_16S_sample-metadata.tsv`  
- `matched_metadata.txt`, `sample-metadata.tsv`  
  → Curated metadata describing each sample (study ID, treatment, location, etc.).

---

## 3. Alpha Diversity
- `all_alpha_diversity.txt`  
  → Results of meta-analysis for alpha diversity across studies.

---

## 4. PICRUSt2 Functional Predictions
Folder: `picrust2_output/`  
- `KO_predicted.tsv.gz` → Predicted KEGG ortholog abundances  
- `EC_predicted.tsv.gz` → Predicted Enzyme Commission abundances  
- `pathways_out/path_abun_unstrat.tsv.gz` → Predicted MetaCyc pathways  
- `marker_predicted_and_nsti.tsv.gz` → Marker genes and NSTI scores  
- `out.tre` → Reference phylogenetic tree  

---

## 5. Reproducible Scripts and Workflow
Folder: `code/`  

### R scripts  
- `alpha_diversity.R` → Calculation of Shannon, Simpson, richness metrics  
- `beta_diversity.R` → CLR transformation, Aitchison distances, dbRDA, ordinations  
- `differential_abundance.R` → DESeq2 analysis for differentially abundant taxa  
- `PICRUSt2_analysis.R` → R workflow for processing PICRUSt2 KO/EC/pathway tables  

### Bash scripts  
- `16SPICRUSt2.sh` → Pipeline for running PICRUSt2 on pooled 16S data  
- `download_all.sh`, `gsa_download.sh` → Scripts for retrieving raw data from repositories  

---

## 6. Supplementary Files
- `Boolean_search_strings.tsv` → Boolean queries used for literature search  
- `Supplementary_Data_S1.xlsx` → Characteristics of included studies (n = 41)  

---

## 7. External Raw Data
- The 16S rRNA gene amplicon data were retrieved from the NCBI SRA and GSA repositories, and can be accessed using the accession IDs included in the download scripts (download_all.sh, gsa_download.sh) provided in the repository.
The Fungal ITS amplicon data were also obtained from the NCBI SRA and GSA repositories, with accession IDs available in the same download scripts.

---

### Notes
- All `.R` scripts are reproducible with R ≥ 4.1 and packages: `phyloseq`, `vegan`, `metafor`, `DESeq2`, `ggplot2`.  
- PICRUSt2 was run with version v2.6.2
- Please refer to the manuscript for methodological details.

---

