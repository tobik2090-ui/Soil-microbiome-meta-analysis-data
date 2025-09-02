#############################################
# Beta Diversity Meta-analysis (16S & ITS)
# dbRDA with Aitchison distance (CLR)
# Output: dbRDA plots (PNG, 600 dpi)
#############################################

# ---- 1. Load packages ----
library(phyloseq)
library(vegan)
library(compositions)   # for CLR transform
library(ggplot2)
library(ggforce)        # for half-violin plots
library(dplyr)

# ---- 2. Helper function: CLR transform phyloseq ----
clr_transform <- function(phy) {
  mat <- as.data.frame(as(otu_table(phy), "matrix"))
  if(taxa_are_rows(phy)) mat <- t(mat)
  mat <- clr(mat + 1)  # pseudocount
  return(mat)
}

# ---- 3. dbRDA analysis function ----
run_dbrda <- function(phy, dataset_label) {
  # CLR + Aitchison distance
  mat <- clr_transform(phy)
  dist_mat <- dist(mat, method = "euclidean")
  
  # Metadata
  meta <- data.frame(sample_data(phy))
  
  # Ensure factors
  meta$Treatment <- factor(meta$Treatment)
  meta$Study <- factor(meta$Project)  # assumes Project = Study
  
  # dbRDA model
  mod <- capscale(dist_mat ~ Treatment + Condition(Study), data = meta)
  
  # Permutation test (stratified by Study)
  anova_res <- anova(mod, permutations = how(nperm = 9999, strata = meta$Study))
  
  # Partial R²
  r2 <- RsquareAdj(mod)$adj.r.squared
  
  # Extract CAP1 scores
  site_scores <- scores(mod, display = "sites", choices = 1)$sites
  df_plot <- cbind(meta, CAP1 = site_scores)
  
  list(model = mod, anova = anova_res, r2 = r2, plotdata = df_plot, dataset = dataset_label)
}

# ---- 4. Load phyloseq objects ----
# Adjust paths
s16 <- readRDS("data/Project_metaanlysis/combined/16S/phyloseq.rds")
sITS <- readRDS("data/Project_metaanlysis/combined/ITS/phyloseq.rds")

# ---- 5. Run dbRDA ----
res16 <- run_dbrda(s16, "16S")
resITS <- run_dbrda(sITS, "ITS")

# ---- 6. Plotting function ----
plot_cap1 <- function(res) {
  ggplot(res$plotdata, aes(x = Treatment, y = CAP1, fill = Treatment)) +
    geom_half_violin(side = "l", alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    theme_minimal(base_size = 12) +
    labs(title = paste0(res$dataset, ": dbRDA CAP1 Scores"),
         subtitle = paste0("Treatment effect: p = ", signif(res$anova$`Pr(>F)`[1], 3),
                           ", adj R² = ", round(res$r2, 3)),
         x = "", y = "CAP1") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
}

# ---- 7. Export plots ----
p16 <- plot_cap1(res16)
pITS <- plot_cap1(resITS)

ggsave("16S_dbRDA_CAP1.png", p16, width = 6, height = 4, dpi = 600)
ggsave("ITS_dbRDA_CAP1.png", pITS, width = 6, height = 4, dpi = 600)
