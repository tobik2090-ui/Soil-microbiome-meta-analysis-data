#############################################
# Alpha Diversity Meta-analysis (16S & ITS)
# Fungicide treatment vs control
# Output: Forest plots (PNG, 600 dpi)
#############################################

# ---- 1. Load packages ----
library(phyloseq)
library(microbiome)
library(dplyr)
library(metafor)

# ---- 2. Load phyloseq objects ----
# Adjust file paths if needed
s16_1  <- readRDS("data/Project_metaanlysis/project1/16S/phyloseq.rds")
s16_2  <- readRDS("data/Project_metaanlysis/project2/16S/phyloseq.rds")
s16_3  <- readRDS("data/Project_metaanlysis/project3/16S/phyloseq.rds")
s16_4  <- readRDS("data/Project_metaanlysis/project4/16S/phyloseq.rds")
s16_5  <- readRDS("data/Project_metaanlysis/project5/16S/phyloseq.rds")

sITS_1 <- readRDS("data/Project_metaanlysis/project1/ITS/phyloseq.rds")
sITS_2 <- readRDS("data/Project_metaanlysis/project2/ITS/phyloseq.rds")
sITS_3 <- readRDS("data/Project_metaanlysis/project3/ITS/phyloseq.rds")
sITS_4 <- readRDS("data/Project_metaanlysis/project4/ITS/phyloseq.rds")
sITS_5 <- readRDS("data/Project_metaanlysis/project5/ITS/phyloseq.rds")

# ---- 3. Extract Shannon diversity ----
get_shannon <- function(phy, project_id) {
  alpha_vals <- microbiome::alpha(phy)
  shannon_vals <- alpha_vals$diversity_shannon
  meta <- data.frame(sample_data(phy))
  df <- cbind(meta, Shannon = shannon_vals)
  df$Project <- project_id
  return(df)
}

df16_1  <- get_shannon(s16_1, "Project1_16S")
df16_2  <- get_shannon(s16_2, "Project2_16S")
df16_3  <- get_shannon(s16_3, "Project3_16S")
df16_4  <- get_shannon(s16_4, "Project4_16S")
df16_5  <- get_shannon(s16_5, "Project5_16S")

dfITS_1 <- get_shannon(sITS_1, "Project1_ITS")
dfITS_2 <- get_shannon(sITS_2, "Project2_ITS")
dfITS_3 <- get_shannon(sITS_3, "Project3_ITS")
dfITS_4 <- get_shannon(sITS_4, "Project4_ITS")
dfITS_5 <- get_shannon(sITS_5, "Project5_ITS")

# ---- 4. Compute effect sizes (Hedges’ g) ----
compute_effect <- function(df) {
  ctrl <- df %>% filter(Treatment == "Control") %>% pull(Shannon)
  trt  <- df %>% filter(Treatment == "Treated") %>% pull(Shannon)
  
  n1 <- length(ctrl); n2 <- length(trt)
  m1 <- mean(ctrl, na.rm = TRUE); m2 <- mean(trt, na.rm = TRUE)
  sd1 <- sd(ctrl, na.rm = TRUE);  sd2 <- sd(trt, na.rm = TRUE)
  
  esc <- escalc(measure = "SMD",
                m1i = m1, sd1i = sd1, n1i = n1,
                m2i = m2, sd2i = sd2, n2i = n2)
  
  data.frame(Project = unique(df$Project), yi = esc$yi, vi = esc$vi)
}

# Author lookup (shared for ITS and 16S)
author_lookup <- data.frame(
  Project = c("Project1_16S", "Project2_16S", "Project3_16S", "Project4_16S", "Project5_16S",
              "Project1_ITS", "Project2_ITS", "Project3_ITS", "Project4_ITS", "Project5_ITS"),
  Author  = c("Chuang et al., 2021", "Ezazi et al., 2021", "Sim et al., 2023",
              "Streletskii et al., 2022", "Sliti et al., 2024",
              "Chuang et al., 2021", "Ezazi et al., 2021", "Sim et al., 2023",
              "Streletskii et al., 2022", "Sliti et al., 2024")
)

# Combine all projects
all_dfs <- list(df16_1, df16_2, df16_3, df16_4, df16_5,
                dfITS_1, dfITS_2, dfITS_3, dfITS_4, dfITS_5)

effect_list <- lapply(all_dfs, compute_effect)
effect_sizes <- bind_rows(effect_list)
effect_sizes <- merge(effect_sizes, author_lookup, by = "Project")

# Split 16S and ITS
effect_16S <- subset(effect_sizes, grepl("16S", Project))
effect_ITS <- subset(effect_sizes, grepl("ITS", Project))

# ---- 5. Run meta-analysis ----
res_16S <- rma(yi, vi, data = effect_16S, method = "REML")
res_ITS <- rma(yi, vi, data = effect_ITS, method = "REML")

# ---- 6. Export forest plots (PNG, 600 dpi) ----

# ITS
png("ITS_ForestPlot_AlphaDiversity.png",
    width = 7, height = 5, units = "in", res = 600)

forest(res_ITS,
       slab = effect_ITS$Author,
       xlab = "Effect size (Hedges' g)",
       main = "ITS: Forest Plot of Fungicide Effect on Alpha Diversity",
       cex = 0.9,
       fonts = "sans",
       mlab = "Random effects model",
       col = "darkgreen")

pval_fmt <- formatC(res_ITS$pval, format = "f", digits = 4)
mtext(paste0("p = ", pval_fmt,
             ", I² = ", round(res_ITS$I2, 1), "%"),
      side = 1, line = 5, cex = 0.9)

dev.off()

# 16S
png("16S_ForestPlot_AlphaDiversity.png",
    width = 7, height = 5, units = "in", res = 600)

forest(res_16S,
       slab = effect_16S$Author,
       xlab = "Effect size (Hedges' g)",
       main = "16S: Forest Plot of Fungicide Effect on Alpha Diversity",
       cex = 0.9,
       fonts = "sans",
       mlab = "Random effects model",
       col = "darkgreen")

pval_fmt <- formatC(res_16S$pval, format = "f", digits = 4)
mtext(paste0("p = ", pval_fmt,
             ", I² = ", round(res_16S$I2, 1), "%"),
      side = 1, line = 5, cex = 0.9)

dev.off()
