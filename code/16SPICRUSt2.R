#!/bin/bash
# ============================================================
# 16S Functional Potential Prediction with PICRUSt2
# # Purpose: Pool all 16S projects, run PICRUSt2, 
#          and generate KO and MetaCyc pathway predictions
# ============================================================

# ---- Step 1: Define input pooled files (generated in R/phyloseq export) ----
FEATURE_TABLE="/home/rstudio/project/all_projects_16S_feature-table.tsv"
SEQS="/home/rstudio/project/all_projects_16S_ASVs.fasta"
OUTDIR="/home/rstudio/project/picrust2_output"

# ---- Step 2: Run PICRUSt2 pipeline ----
# -i = ASV representative sequences
# -t = feature/OTU abundance table
# -o = output directory
# -p = number of threads (adjust to your server)
picrust2_pipeline.py \
-s $SEQS \
-i $FEATURE_TABLE \
-o $OUTDIR \
-p 4

# ---- Step 3: Main outputs ----
# $OUTDIR/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz   → KO abundances
# $OUTDIR/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz   → Enzyme class abundances
# $OUTDIR/pathways_out/path_abun_unstrat.tsv.gz              → MetaCyc pathway abundances
# $OUTDIR/marker_predicted_and_nsti.tsv.gz                   → NSTI values (prediction accuracy)

# ---- Step 4: Notes ----
# • Input pooled files were exported from pooled phyloseq objects in R.
# • KO and MetaCyc outputs were used for downstream statistical analyses in R.
# • NSTI values should be checked to interpret prediction confidence.
