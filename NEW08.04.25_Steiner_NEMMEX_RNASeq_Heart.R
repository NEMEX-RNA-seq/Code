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

# --- Pre-Process ------------------
counts <- readxl::read_excel("./gene_count.xls.xlsx")
exp <- counts %>%
  dplyr::select(gene_id, starts_with("HT"), starts_with("L"), starts_with("P")) %>%
  tibble::column_to_rownames("gene_id")

fdata <- counts %>% #gives you information about your genes
  dplyr::select(-c(
    dplyr::starts_with("HT"),
    dplyr::starts_with("PF"),
    dplyr::starts_with("L")
    ))
rownames(fdata) <- fdata$gene_id #making row names the gene IDs to match exp data

meta <- data.frame(
  sample = colnames(exp)
) %>%
  mutate(tissue=case_when(
    grepl("^L", sample) ~ "Lung", 
    grepl("^P", sample) ~ "Pre-Frontal Cortex", 
    grepl("^HT", sample) ~ "Heart"
  )) %>%
  separate(sample, into = c("tissue_short", "metal", "id"), sep = "_", remove = FALSE)%>%
  mutate(sex = case_when(
    grepl("^M", id) ~ "Male", 
    grepl("^F", id) ~ "Female"
  ))
rownames(meta) <- meta$sample

###checking to make sure rownames are the same
# ensure sample order matches and stop if it does not
stopifnot(all(rownames(meta) == colnames(exp)))

# ensure that feature order matches and stop if it does not
stopifnot(all(rownames(fdata) == rownames(exp)))

# Create SummarizedExperiment Object
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(exp)), 
  colData = meta, 
  rowData = fdata)
# Filter out to keep heart data
se <- se[, grepl("HT", se$tissue_short)]

# --- PCA --------------
se <- se %>%
  identify_abundant(factor_of_interest = metal) %>% #looking at what are the abundant genes within each metal
  scale_abundance() %>%
  identify_abundant(minimum_counts = 5, 
                    minimum_proportion = 0.5) %>%
  keep_abundant() #only keeping most abundant genes, filtering out lowly expressed genes (Novogene did not do this), eliminates some noise

##make PCA object
se <- se %>%
  reduce_dimensions(method = "PCA", 
                    .dims = 20, 
                    top = 20000)
## PCA Plot
se %>%
  pivot_sample() %>%
  ggplot(aes(PC1, PC2, color = metal, label=se$sample)) +
  geom_point() +
  geom_text_repel()+
  theme_pubr(legend = "right") +
  scale_color_npg() + 
  labs(color = "metal")

##removing outliers
se <- se[, !grepl("HT_Pb_M6", se$sample)]
se <- se[, !grepl("HT_Pb_F5", se$sample)]
se <- se[, !grepl("HT_Un_F4", se$sample)]

##checking to make sure outliers are removed
##running PCA again
se <- se %>%
  identify_abundant(factor_of_interest = metal) %>%
  scale_abundance() %>%
  identify_abundant(minimum_counts = 5, 
                    minimum_proportion = 0.5) %>%
  keep_abundant()

se <- se %>%
  reduce_dimensions(method = "PCA", 
                    .dims = 20, 
                    top = 20000)
## PCA Plot without outliers
se %>%
  pivot_sample() %>%
  ggplot(aes(PC1, PC2, color = metal)) +
  geom_point() +
  theme_pubr(legend = "right") +
  scale_color_npg() + 
  ggtitle("PCA Plot of Heart RNASeq Samples without Outliers")

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

## --- comparisons between metals/mix and untreated, need to do by each individual sex/metal --------------------
se[,se$sex=="Male" &                #making a new SummarizedExperiment with only Male samples and As/Un for metals
     grepl("As|Un", se$metal)] %>%
  (\(se){                           #making a new function (\ indicates function), where we are reordering the levels to make Un the reference level
    se$metal = factor(se$metal, levels=c("Un", "As"));
    se})()
                            
M_As_M_Un <- compare_to(
  se=se[,se$sex=="Male" & 
          grepl("As|Un", se$metal)] %>%
    (\(se){
    se$metal = factor(se$metal, levels=c("Un", "As"));
    se})(), 
  formula = ~metal
)

se[,se$sex=="Male" & 
     grepl("Cd|Un", se$metal)] %>%
  (\(se){
    se$metal = factor(se$metal, levels=c("Un", "Cd"));
    se})()

M_Cd_M_Un <- compare_to(
  se=se[,se$sex=="Male" & 
          grepl("Cd|Un", se$metal)] %>%
    (\(se){
      se$metal = factor(se$metal, levels=c("Un", "Cd"));
      se})(), 
  formula = ~metal
)

# Male Pb vs Un
se[, se$sex == "Male" & grepl("Pb|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })()

M_Pb_M_Un <- compare_to(
  se = se[, se$sex == "Male" & grepl("Pb|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })(),
  formula = ~metal
)

# Male Cr vs Un
se[, se$sex == "Male" & grepl("Cr|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })()

M_Cr_M_Un <- compare_to(
  se = se[, se$sex == "Male" & grepl("Cr|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })(),
  formula = ~metal
)

# Male Mix vs Un
se[, se$sex == "Male" & grepl("Mx|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })()

M_Mix_M_Un <- compare_to(
  se = se[, se$sex == "Male" & grepl("Mx|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })(),
  formula = ~metal
)

# Female As vs Un
se_F_As_Un <- se[, se$sex == "Female" & grepl("As|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "As")); se })()

F_As_F_Un <- compare_to(
  se = se_F_As_Un,
  formula = ~metal
)

# Female Cd vs Un
se_F_Cd_Un <- se[, se$sex == "Female" & grepl("Cd|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cd")); se })()

F_Cd_F_Un <- compare_to(
  se = se_F_Cd_Un,
  formula = ~metal
)

# Female Pb vs Un
se[, se$sex == "Female" & grepl("Pb|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })()

F_Pb_F_Un <- compare_to(
  se = se[, se$sex == "Female" & grepl("Pb|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Pb")); se })(),
  formula = ~metal
)

# Female Cr vs Un
se[, se$sex == "Female" & grepl("Cr|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })()

F_Cr_F_Un <- compare_to(
  se = se[, se$sex == "Female" & grepl("Cr|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Cr")); se })(),
  formula = ~metal
)

# Female Mix vs Un
se[, se$sex == "Female" & grepl("Mx|Un", se$metal)] %>%
  (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })()

F_Mix_F_Un <- compare_to(
  se = se[, se$sex == "Female" & grepl("Mx|Un", se$metal)] %>%
    (\(se) { se$metal <- factor(se$metal, levels = c("Un", "Mx")); se })(),
  formula = ~metal
)

deg_df <- list(
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
  bind_rows(.id = "condition")

# Optional filtering
sig_df <- deg_df %>%
  filter(pvalue < 0.05, abs(log2FoldChange) > log2(1.5)) #alternate between pvalue and padj

# --- compareCluster and dotplots --------------
#use compareCluster to make all of the comparisons between metals for males and females (instead of enrichGO) by condition
##cannot use enrichGO because we are making so many comparisons
metalcluster <- compareCluster(
  transcript ~ condition,
  data = sig_df,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,  
  keyType = "ENSEMBL",   
  ont = "BP",           
  pAdjustMethod = "BH", #controlling FDR
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
##visualizing results
dotplot(metalcluster, showCategory = 5) + ggtitle("GO Enrichment by Condition")
###too long of a dotplot and hard to read, going to split by female and male 

###splitting male and female sig degs
sig_df_female <- sig_df %>% filter(str_detect(condition, "^Female"))
sig_df_male <- sig_df %>% filter(str_detect(condition, "^Male"))

# Female enrichment
female_cluster <- compareCluster(
  transcript ~ condition,
  data = sig_df_female,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# Male enrichment
male_cluster <- compareCluster(
  transcript ~ condition,
  data = sig_df_male,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType= "ENSEMBL", 
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

###new male and female dotplots
dotplot(female_cluster, showCategory = 10) + ggtitle("GO Enrichment - Female Comparisons")
dotplot(male_cluster, showCategory = 10) + ggtitle("GO Enrichment - Male Comparisons")

# --- sex differences comparecluster ----------------------------
#making the dotplots for each metal with male and female comparisons together
metalcluster_enrres <- compareCluster(
  transcript ~ condition,
  data = sig_df,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

##dotplot for arsenic
dotplot(metalcluster_enrres %>%
          filter(grepl("As", Cluster)), 
        showCategory = 10) + ggtitle("GO Enrichment - Female vs Male As")
###no significant male degs

##dotplot for cadmium - N/A

##dotplot for chromium
dotplot(metalcluster_enrres %>%
          filter(grepl("Cr", Cluster)), 
        showCategory = 10) + ggtitle("GO Enrichment - Female vs Male Cr")

##dotplot for lead
dotplot(metalcluster_enrres %>%
          filter(grepl("Pb", Cluster)), 
        showCategory = 10) + ggtitle("GO Enrichment - Female vs Male Pb")

##dotplot for mix
dotplot(metalcluster_enrres %>%
          filter(grepl("Mix", Cluster)), 
        showCategory = 10) + ggtitle("GO Enrichment - Female vs Male Mix")
###no significant female degs

##dotplot for untreated
dotplot(metalcluster_enrres %>%
          filter(grepl("Un", Cluster)), 
        showCategory = 5) + ggtitle("GO Enrichment - Female vs Male Un")

##see what terms overlap for all comparisons
commonterms <- metalcluster_enrres@compareClusterResult %>%
  group_by(Description) %>%
  reframe(n=length(unique(Cluster))) %>%
  filter(n>1) %>%
  pull(Description)

#only showing top 10 overlapping terms for common biology for all comparisons
dotplot(metalcluster_enrres %>%
          filter(Description %in% commonterms), 
        showCategory = 10) + ggtitle("GO Enrichment - Common Biology")

# --- Volcano plots ----------------
##female comparisons first
EnhancedVolcano::EnhancedVolcano(
  toptable = F_As_F_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_As_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F As vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Cd_F_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_Cd_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Cd vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Cr_F_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_Cr_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Cr vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Pb_F_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_Pb_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Pb vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = F_Mix_F_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=F_Mix_F_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of F Mix vs Un RNA-Seq Data") # subtitle

##male comparisons
EnhancedVolcano::EnhancedVolcano(
  toptable = M_As_M_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=M_As_M_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M As vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Cd_M_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=M_Cd_M_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Cd vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Cr_M_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=M_Cr_M_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Cr vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Pb_M_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=M_Pb_M_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Pb vs Un RNA-Seq Data") # subtitle

EnhancedVolcano::EnhancedVolcano(
  toptable = M_Mix_M_Un,                                   # name of results df
  x = "log2FoldChange",                                    # name of log2FC column
  y="pvalue",                                                # name of p-value column
  pCutoff = 0.05,                                          # p-value cutoff
  FCcutoff = log2(2),                                      # fold change cutoff
  lab=M_Mix_M_Un$gene_name,                                # gene names as label                    
  title = "Volcano Plot",                                  # title
  col = c('black', 'pink', 'purple', 'red3'),              # pallete colors
  drawConnectors = TRUE,                                   # box the labels
  subtitle = "Differential Expression of M Mix vs Un RNA-Seq Data") # subtitle











