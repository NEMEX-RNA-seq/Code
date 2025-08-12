# --- Loading in Required Packages ----------------------
library(tidyverse)
library(tidybulk)
library(readxl)
library(SummarizedExperiment)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(boot)

# --- Save Summarized Experiment ------------------
se <- readRDS(file = "./results/se.rds")

# --- Differential Expression ---------------
# Define a function to take a summarized experiment object and performs differential expression using DESeq2
compare_to <- function(
    se,
    formula) {
  se |> 
    test_differential_abundance(formula,               #testing differential abundance with deseq2 in tidybulk package
                                method = "deseq2",
                                fitType = "local",     #fitType = local is method to estimate dispersions
                                action = "get")  |>  #action = get means we will return the DE results
    as.data.frame() |>                                 #converts results into a data frame that we can use later
    inner_join(rowData(se) |>                         # matches DE results with gene annotations in Summarized Experiment
                 as.data.frame(), 
               by = c("transcript"="gene_id")) |> 
    mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down"))  #adding direction
}

## --- Bootstrap - Differential Expression ----------------------
bootstrap_de <- function(
    se, 
    formula, nboot) {
  total_res <- data.frame()
  for(i in 1:nboot) {
    newsamples <- sample(colnames(se), replace = TRUE)
    se_boot <- se[,newsamples]
    de_res <- compare_to(se_boot, formula) %>%
      mutate(bootstrap_number = i)
    total_res <- rbind(total_res, de_res)
  } 
  total_res <- total_res %>%
    group_by(transcript) %>%
    mutate(presence = mean(pvalue <0.05 & abs(log2FoldChange)>log2(1.5))) %>% #gives you how often the gene was significant (DEG)
    ungroup()
}

## --- Bootstrap comparisons --------------------
###using bootstrapping to test how often a gene is significant while sampling different samples with replacement (given by the presence column)
###going to run and save all of these so that we can bring them back for comparisons later
M_As_M_Un <- bootstrap_de(
  se=se[,se$sex=="Male" & 
          grepl("As|Un", se$metal)] %>%
    (\(se){
      se$metal = factor(se$metal, levels=c("Un", "As"));
      se})(), 
  formula = ~metal, 
  nboot = 50
)

hist(M_As_M_Un$presence)
saveRDS(M_As_M_Un, file = "./results/M_As_M_Un.rds")

# Male Cd vs Un
M_Cd_M_Un <- bootstrap_de(
  se=se[,se$sex=="Male" & 
          grepl("Cd|Un", se$metal)] %>%
    (\(se){
      se$metal = factor(se$metal, levels=c("Un", "Cd"));
      se})(), 
  formula = ~metal, 
  nboot = 50
)
hist(M_Cd_M_Un$presence)
saveRDS(M_Cd_M_Un, file = "./results/M_Cd_M_Un.rds")

# Male Pb vs Un
M_Pb_M_Un <- bootstrap_de(
  se = se[, se$sex == "Male" & grepl("Pb|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(M_Pb_M_Un$presence)
saveRDS(M_Pb_M_Un, file = "./results/M_Pb_M_Un.rds")

# Male Cr vs Un
M_Cr_M_Un <- bootstrap_de(
  se = se[, se$sex == "Male" & grepl("Cr|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(M_Cr_M_Un$presence)
saveRDS(M_Cr_M_Un, file = "./results/M_Cr_M_Un.rds")

# Male Mix vs Un
M_Mix_M_Un <- bootstrap_de(
  se = se[, se$sex == "Male" & grepl("Mx|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(M_Mix_M_Un$presence)
saveRDS(M_Mix_M_Un, file = "./results/M_Mix_M_Un.rds")

# Female As vs Un
F_As_F_Un <- bootstrap_de(
  se = se[, se$sex == "Female" & grepl("As|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "As")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(F_As_F_Un$presence)
saveRDS(F_As_F_Un, file = "./results/F_As_F_Un.rds")

# Female Cd vs Un
F_Cd_F_Un <- bootstrap_de(
  se = se[, se$sex == "Female" & grepl("Cd|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cd")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(F_Cd_F_Un$presence)
saveRDS(F_Cd_F_Un, file = "./results/F_Cd_F_Un.rds")

# Female Pb vs Un
F_Pb_F_Un <- bootstrap_de(
  se = se[, se$sex == "Female" & grepl("Pb|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(F_Pb_F_Un$presence)
saveRDS(F_Pb_F_Un, file = "./results/F_Pb_F_Un.rds")

# Female Cr vs Un
F_Cr_F_Un <- bootstrap_de(
  se = se[, se$sex == "Female" & grepl("Cr|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })(),
  formula = ~metal,
  nboot = 50
)
hist(F_Cr_F_Un$presence)
saveRDS(F_Cr_F_Un, file = "./results/F_Cr_F_Un.rds")

# Female Mix vs Un
F_Mix_F_Un <- bootstrap_de(
  se = se[, se$sex == "Female" & grepl("Mx|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })(),
  formula = ~metal, 
  nboot = 50
)
hist(F_Mix_F_Un$presence)
saveRDS(F_Mix_F_Un, file = "./results/F_Mix_F_Un.rds")

# --- Load in bootstrapping results here! -----------
M_As_M_Un <- readRDS(file = "./results/M_As_M_Un.rds")
M_Cd_M_Un <- readRDS(file = "./results/M_Cd_M_Un.rds")
M_Pb_M_Un <- readRDS(file = "./results/M_Pb_M_Un.rds")
M_Cr_M_Un <- readRDS(file = "./results/M_Cr_M_Un.rds")
M_Mix_M_Un <- readRDS(file = "./results/M_Mix_M_Un.rds")
F_As_F_Un <- readRDS(file = "./results/F_As_F_Un.rds")
F_Cd_F_Un <- readRDS(file = "./results/F_Cd_F_Un.rds")
F_Pb_F_Un <- readRDS(file = "./results/F_Pb_F_Un.rds")
F_Cr_F_Un <- readRDS(file = "./results/F_Cr_F_Un.rds")
F_Mix_F_Un <- readRDS(file = "./results/F_Mix_F_Un.rds")

# --- Bootstrapping Significant DEGs --------------
boot_df <- list(
  "Male As vs Un" = M_As_M_Un,
  "Male Cd vs Un" = M_Cd_M_Un,
  "Male Pb vs Un" = M_Pb_M_Un,
  "Male Cr vs Un" = M_Cr_M_Un,
  "Male Mix vs Un" = M_Mix_M_Un,
  "Female As vs Un" = F_As_F_Un,
  "Female Cd vs Un" = F_Cd_F_Un,
  "Female Pb vs Un" = F_Pb_F_Un,
  "Female Cr vs Un" = F_Cr_F_Un,
  "Female Mix vs Un" = F_Mix_F_Un
) %>%
  bind_rows(.id = "condition") %>%
  filter(presence > 0.3) %>% #filtering out genes that were significant more than 30% of the time
  dplyr::select(condition, transcript) %>%
  distinct()

# --- CompareCluster Results After Bootstrapping -------------
metalcluster_enrres_boot <- compareCluster(
  transcript ~ condition,
  data = boot_df,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

##arsenic comparisons
dotplot(metalcluster_enrres_boot %>%
          filter(grepl("As", Cluster)), 
        showCategory = 10) + ggtitle("Bootstrapped GO Enrichment - Female vs Male As")
###no significant male degs

#cadmium comparisons
dotplot(metalcluster_enrres_boot %>%
          filter(grepl("Cd", Cluster)), 
        showCategory = 10) + ggtitle("Bootstrapped GO Enrichment - Female vs Male Cd")
###no significant male degs

#lead comparisons
dotplot(metalcluster_enrres_boot %>%
          filter(grepl("Pb", Cluster)), 
        showCategory = 10) + ggtitle("Bootstrapped GO Enrichment - Female vs Male Pb")

#chromium comparisons
dotplot(metalcluster_enrres_boot %>%
          filter(grepl("Cr", Cluster)), 
        showCategory = 10) + ggtitle("Bootstrapped GO Enrichment - Female vs Male Cr")

#mix comparisons
dotplot(metalcluster_enrres_boot %>%
          filter(grepl("Mix", Cluster)), 
        showCategory = 10) + ggtitle("Bootstrapped GO Enrichment - Female vs Male Mix")

###splitting male and female sig degs
sig_boot_df_female <- boot_df %>% filter(str_detect(condition, "^Female"))
sig_boot_df_male <- boot_df %>% filter(str_detect(condition, "^Male"))

# Female enrichment
female_boot_cluster <- compareCluster(
  transcript ~ condition,
  data = sig_boot_df_female,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# Male enrichment
male_boot_cluster <- compareCluster(
  transcript ~ condition,
  data = sig_boot_df_male,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType= "ENSEMBL", 
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)
###new male and female dotplots
dotplot(female_boot_cluster, showCategory = 10) + ggtitle("GO Enrichment - Female Comparisons")
dotplot(male_boot_cluster, showCategory = 10) + ggtitle("GO Enrichment - Male Comparisons")

# --- Volcano plots ----------------
#change to median log2foldchange and medianpvalue for all
##female comparisons first
EnhancedVolcano::EnhancedVolcano(
  toptable = (F_As_F_Un %>%                               # name of results df
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue))),             
  x = "median_logFC",                                      # name of log2FC column
  y="median_pvalue",                                       # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=((F_As_F_Un %>%                               # name of results df
          group_by(gene_name) %>%
          reframe(median_logFC = median(log2FoldChange), 
                  median_pvalue = median(pvalue))))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F As vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Cd_F_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(F_Cd_F_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Cd vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Cr_F_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_Cr_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Cr vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Pb_F_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(F_Pb_F_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Pb vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Mix_F_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(F_Mix_F_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Mix vs Un RNA-Seq Data") # subtitle

##male comparisons
EnhancedVolcano::EnhancedVolcano(
  toptable = M_As_M_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(M_As_M_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M As vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Cd_M_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(M_Cd_M_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Cd vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Cr_M_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(M_Cr_M_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Cr vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Pb_M_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(M_Pb_M_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Pb vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Mix_M_Un %>%
    group_by(gene_name) %>%
    reframe(median_logFC = median(log2FoldChange), 
            median_pvalue = median(pvalue)),                                   # name of results df
  x = "median_logFC",                                    # name of log2FC column
  y="median_pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=(M_Mix_M_Un %>%
         group_by(gene_name) %>%
         reframe(median_logFC = median(log2FoldChange), 
                 median_pvalue = median(pvalue)))$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Mix vs Un RNA-Seq Data") # subtitle

# --- Deconvolution Test ----------------
deconv_test <- se %>%
  deconvolve_cellularity(
    sample = sample, 
    transcript = gene_id, 
    method = "cibersort", 
    action="get")










