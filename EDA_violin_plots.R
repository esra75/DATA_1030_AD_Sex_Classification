library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gt)
library(dplyr)
library(tidyverse)
library(edgeR)
library(uwot) 
library(ggplot2)
library(ggfortify)
library(knitr)
library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1. Load and annotate metadata
setwd("/oscar/data/larschan/etaner1/DATA1030_ML/tsv")
df <- read.delim("experiment_report_2025_120.tsv", header=TRUE, sep="\t", quote="", check.names=FALSE, skip=1)
df$Sex <- ifelse(grepl("female", df$`Biosample summary`, ignore.case = TRUE), "female", "male")
df$diagnosis <- ifelse(grepl("Alzheimer", df$`Biosample summary`, ignore.case = TRUE), "Alzheimer",
                       ifelse(grepl("mild cognitive impairment", df$`Biosample summary`, ignore.case = TRUE), "MCI",
                              ifelse(grepl("Cognitive impairment", df$`Biosample summary`, ignore.case = TRUE), "Cognitive impairment", "Control")))
metadata <- df %>%
  select(Accession, Sex, diagnosis) %>%
  arrange(Accession)

# 2. Read and merge raw count data (120 columns)
setwd("/oscar/data/larschan/etaner1/DATA1030_ML/count_matrix_120")
files <- list.files(pattern="\\.tsv$")
read_tsv_counts <- function(f) {
  dat <- read.delim(f)
  dat <- dat[, c("gene_id", "expected_count")]
  colnames(dat)[2] <- gsub(".tsv", "", f)
  return(dat)
}
merged <- reduce(lapply(files, read_tsv_counts), full_join, by="gene_id")
rownames(merged) <- merged$gene_id
merged <- merged[, -1] # drop gene_id column
colnames(merged) <- metadata$Accession

# 3. Filter, normalize, and calculate logCPM (NO collapse step)
dge <- DGEList(counts = merged)
keep <- filterByExpr(dge, group = metadata$Sex)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE)
gene_ids_noversion <- sub("\\..*", "", rownames(logCPM))
rownames(logCPM) <- gene_ids_noversion

# 4. Map Ensembl IDs, etc.
symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids_noversion,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)
symbols <- symbols[!is.na(symbols$SYMBOL), ]
sex_genes <- c("XIST", "UTY", "DDX3Y")
ids_sex_genes <- symbols$ENSEMBL[symbols$SYMBOL %in% sex_genes]

# 5. Extract expression for these genes and make plot dataframe
expr_interest <- as.data.frame(t(logCPM[rownames(logCPM) %in% ids_sex_genes, ]))
colnames(expr_interest) <- symbols$SYMBOL[match(colnames(expr_interest), symbols$ENSEMBL)]
expr_interest$Sex <- metadata$Sex
expr_long <- expr_interest %>%
  pivot_longer(cols = sex_genes, names_to = "Gene", values_to = "Expression")

# 6. Plot
ggplot(expr_long, aes(x = Sex, y = Expression, fill = Sex)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = c("female" = "#E94B8A", "male" = "#4A90E2")) +
  labs(
    title = "Expression of Canonical Sex-Linked Genes by Sex",
    y = "logCPM (Normalized Expression)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 12)
  )


# sex gene summary stats
sex_gene_summary <- expr_long %>%
  group_by(Gene, Sex) %>%
  summarise(
    n = n(),
    mean = mean(Expression, na.rm = TRUE),
    median = median(Expression, na.rm = TRUE),
    sd = sd(Expression, na.rm = TRUE),
    min = min(Expression, na.rm = TRUE),
    q25 = quantile(Expression, 0.25, na.rm = TRUE),
    q75 = quantile(Expression, 0.75, na.rm = TRUE),
    max = max(Expression, na.rm = TRUE)
  ) %>%
  ungroup()

print(sex_gene_summary)


# write out 
kable(
  sex_gene_summary,
  digits = 2,
  caption = "Summary of Sex-Linked Gene Expression"
)

write.csv(
  sex_gene_summary,
  "sex_gene_summary.csv",
  row.names = FALSE
)



# PCA & UMAP for extra validation 
# make sure logCPM has genes as rows and samples as columns:
cat(dim(logCPM), "\n")


# 1. pca using ALL expressed genes 
pca_all <- prcomp(t(logCPM), scale. = TRUE)

autoplot(
  pca_all,
  data = metadata,
  colour = "Sex",
  shape = "diagnosis",
  size = 3
) +
  ggtitle("PCA: All Genes Colored by Sex") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# 2. PCA USING ONLY SEX-LINKED GENES
# get expression for the sex-linked genes
sex_expr <- t(logCPM[rownames(logCPM) %in% symbols$ENSEMBL[symbols$SYMBOL %in% sex_genes], ])

pca_sex <- prcomp(sex_expr, scale. = TRUE)
autoplot(
  pca_sex,
  data = metadata,
  colour = "Sex",
  size = 3
) +
  ggtitle("PCA: Canonical Sex-Linked Genes (XIST, UTY, DDX3Y)") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# 3. UMAP using ALL genes
set.seed(123)
umap_all <- umap(t(logCPM), n_neighbors = 15, min_dist = 0.2, metric = "cosine")

umap_all_df <- as.data.frame(umap_all)
colnames(umap_all_df) <- c("UMAP1", "UMAP2")
umap_all_df$Sex <- metadata$Sex
umap_all_df$Diagnosis <- metadata$diagnosis

ggplot(umap_all_df, aes(UMAP1, UMAP2, color = Sex, shape = Diagnosis)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  ggtitle("UMAP: All Genes Colored by Sex") +
  theme(text = element_text(size = 14), plot.title = element_text(face = "bold", hjust = 0.5))


# 4. UMAP USING ONLY SEX-LINKED GENES
set.seed(123)
umap_sex <- umap(sex_expr, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

umap_sex_df <- as.data.frame(umap_sex)
colnames(umap_sex_df) <- c("UMAP1", "UMAP2")
umap_sex_df$Sex <- metadata$Sex

ggplot(umap_sex_df, aes(UMAP1, UMAP2, color = Sex)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  ggtitle("UMAP: Sex-Linked Genes Only (XIST, UTY, DDX3Y)") +
  theme(text = element_text(size = 14), plot.title = element_text(face = "bold", hjust = 0.5))


# oof these looks not good :( 
# modified Drosophila code --> plot only X chr exp ? 
# ugh the ENSEMBL database won't link 
