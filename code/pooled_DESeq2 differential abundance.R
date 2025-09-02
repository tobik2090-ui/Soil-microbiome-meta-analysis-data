##############################################################
# Meta-analysis of soil microbiome responses to fungicide
#  - 16S (bacteria/archaea) and ITS (fungi)
#  - Pooled across 5 projects
#  - DESeq2 differential abundance (barplots)
#  - PCoA (Bray–Curtis) beta diversity analysis
##############################################################

# ---- Load packages ----
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(plyr)
library(vegan)

# ---- Helper function to merge projects ----
merge_amplicon <- function(amplicon = "16S") {
  ps_list <- list()
  for (i in 1:5) {
    ps <- readRDS(paste0("project", i, "/", amplicon, "/phyloseq.rds"))
    taxa <- readRDS(paste0("project", i, "/", amplicon, "/taxa_new.rds"))
    tax_table(ps) <- tax_table(as.matrix(taxa))
    ps_list[[i]] <- ps
  }
  do.call(merge_phyloseq, ps_list)
}

# ---- 1. Merge phyloseq objects ----
ps_16S_all <- merge_amplicon("16S")
ps_ITS_all <- merge_amplicon("ITS")

# ---- 2. PCoA plots (Bray–Curtis) ----
make_pcoa <- function(ps, title, file) {
  ord <- ordinate(ps, method = "PCoA", distance = "bray")
  cb_palette <- c("Control" = "#E69F00", "Treated" = "#56B4E9")
  p <- plot_ordination(ps, ord, color = "Treatment") +
    geom_point(size = 4, shape = 21, stroke = 1, color = "black", aes(fill = Treatment)) +
    scale_fill_manual(values = cb_palette) +
    scale_color_manual(values = cb_palette) +
    labs(
      title = title,
      fill = "Treatment",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    )
  ggsave(file, p, dpi = 900, width = 8, height = 6)
}

make_pcoa(ps_16S_all, "Pooled PCoA of 16S community structure (Bray–Curtis)", "PCoA_16S_Treatment.png")
make_pcoa(ps_ITS_all, "Pooled PCoA of ITS community structure (Bray–Curtis)", "PCoA_ITS_Treatment.png")

# ---- 3. PERMANOVA ----
adonis_test <- function(ps) {
  bc <- phyloseq::distance(ps, method = "bray")
  md <- data.frame(sample_data(ps))
  adonis2(bc ~ Treatment + Project, data = md, permutations = 999)
}
adonis_16S <- adonis_test(ps_16S_all)
adonis_ITS <- adonis_test(ps_ITS_all)
print(adonis_16S)
print(adonis_ITS)

# ---- 4. DESeq2 differential abundance ----
make_deseq <- function(ps, amplicon, p_filter = 0.05, lfc_cutoff = NULL, fdr = FALSE) {
  ps_genus <- tax_glom(ps, "Genus")
  keep <- colSums(otu_table(ps_genus) > 0) >= 5
  ps_genus <- prune_taxa(keep, ps_genus)
  ps_genus <- prune_taxa(colSums(otu_table(ps_genus)) > 0, ps_genus)
  ps_genus <- prune_samples(rowSums(otu_table(ps_genus)) > 0, ps_genus)
  
  # DESeq2
  dds <- phyloseq_to_deseq2(ps_genus, ~ Treatment)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, fitType = "parametric")
  res <- results(dds, contrast = c("Treatment", "Treated", "Control"))
  
  # Annotate taxonomy
  taxa_tab <- as.data.frame(tax_table(ps_genus))
  res_df <- as.data.frame(res)
  res_df$ASV <- rownames(res_df)
  res_df <- cbind(res_df, taxa_tab[res_df$ASV, c("Genus", "Phylum")])
  res_df$Phylum <- gsub("^p__", "", res_df$Phylum)
  res_df$Genus  <- gsub("^g__", "", res_df$Genus)
  res_df$Genus  <- gsub("^g_", "", res_df$Genus)
  
  # Apply filters
  if (fdr) {
    sig <- res_df[!is.na(res_df$padj) & res_df$padj < p_filter & !is.na(res_df$Genus) & res_df$Genus != "", ]
  } else {
    sig <- res_df[!is.na(res_df$pvalue) & res_df$pvalue < p_filter & !is.na(res_df$Genus) & res_df$Genus != "", ]
  }
  if (!is.null(lfc_cutoff)) {
    sig <- sig[abs(sig$log2FoldChange) > lfc_cutoff, ]
  }
  sig$Genus <- make.unique(sig$Genus)
  
  # Plot
  p <- ggplot(sig, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = Phylum)) +
    geom_bar(stat = "identity", color = "black") +
    coord_flip() +
    labs(x = "Genus", y = "Log2 Fold Change") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(face = "bold", size = 17),
      legend.text  = element_text(face = "bold", size = 15)
    )
  ggsave(paste0(amplicon, "_DESeq2.png"), p, dpi = 900, width = 8, height = 6)
  
  return(sig)
}

# 16S: loose filter, unadjusted p < 0.05, |log2FC| > 1.5
sig16 <- make_deseq(ps_16S_all, "16S", p_filter = 0.05, lfc_cutoff = 1.5, fdr = FALSE)

# ITS: FDR strict (padj < 0.05)
sigITS <- make_deseq(ps_ITS_all, "ITS", p_filter = 0.05, lfc_cutoff = NULL, fdr = TRUE)

# ITS: loose (unadjusted p < 0.05)
sigITS_loose <- make_deseq(ps_ITS_all, "ITS_loose", p_filter = 0.05, lfc_cutoff = NULL, fdr = FALSE)

##############################################################
# End of Script
##############################################################
