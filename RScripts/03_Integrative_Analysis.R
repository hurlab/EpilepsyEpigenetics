################################################################################
# Integrative Analysis: RNA-Seq and RRBS Data Integration
################################################################################
# This script integrates differential gene expression (DEG) and differential
# methylation (DMG) data to identify genes with coordinated changes in
# methylation and expression across two epileptogenesis models.
#
# Analysis workflow:
#   1. Load DEG and DMG results from previous analyses
#   2. Within-model DEG-DMG overlaps (KI and KA separately)
#   3. Cross-model concordant genes (shared DEGs and DMGs)
#   4. Scatter plots of methylation vs expression changes
#   5. Directional categorization (Hypo-Down, Hypo-Up, Hyper-Down, Hyper-Up)
#   6. CpG-site level visualization for key genes
#   7. Model-specific correlation analysis
#   8. Permutation testing for overlap significance
#
# Author: [Your Name]
# Date: December 2025
################################################################################

################################################################################
# 1. Initialize Workspace
################################################################################

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (nchar(script_path) > 0) {
    setwd(dirname(script_path))
  }
}

cat("Working directory:", getwd(), "\n")

################################################################################
# 2. Load Required Packages
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  # Data manipulation
  dplyr, tidyr, data.table, stringr,

  # Visualization
  ggplot2, ggpubr, pheatmap, RColorBrewer,
  gridExtra, grid, ComplexHeatmap,

  # Genomic data
  GenomicRanges, genomation,

  # Statistics
  VennDetail, UpSetR,

  # File I/O
  openxlsx, readxl
)

################################################################################
# 3. Define Global Variables and Directories
################################################################################

analysis_date <- '250502'
padj_threshold <- 0.05
qvalue_threshold <- 0.01
meth_diff_threshold <- 25

# Directories
output_dir <- 'output_Integrative/'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input directories from previous analyses
rnaseq_dir <- 'output_RNASeq/'
rrbs_dir <- 'output_RRBS/'
dmc_dir <- 'output_250805/'

# Data files
sample_file <- '../GEOData/samples_nooutlier.txt'
fpkm_file <- '../GEOData/fpkmDatawithMean_220205.txt'
beta_file <- '../GEOData/beta_values.txt.gz'

################################################################################
# 4. Load Sample Information
################################################################################

cat("\n=== Loading sample metadata ===\n")

samples <- read.delim(sample_file, sep = "\t", header = TRUE)
samples <- samples[!grepl("5ADs|1stlimbic", samples$group), ]

# Create sample mapping for models
samples$model <- ifelse(grepl("3ClassV1|Sham$", samples$group), "KI", "KA")
samples$condition <- ifelse(grepl("ClassV", samples$group), "Seizure", "Control")

cat("Sample summary:\n")
print(table(samples$model, samples$condition))

################################################################################
# 5. Load DEG Results
################################################################################

cat("\n=== Loading DEG results ===\n")

# Load DEG tables
ki_deg_file <- file.path(rnaseq_dir, "KI_3ClassV_vs_KI_Sham_DEtable_padj0.05_250502.txt")
ka_deg_file <- file.path(rnaseq_dir, "KA_3ClassV_vs_KA_ShamSaline_DEtable_padj0.05_250502.txt")

# Check alternative locations if not found in output_RNASeq
if (!file.exists(ki_deg_file)) {
  ki_deg_file <- "../220205_Madhu_RNASeq_ReAnalysis_KI-KA-DBonly/output_250818/KI_3ClassV_vs_KI_Sham_DEtable_padj0.05_250502.txt"
}
if (!file.exists(ka_deg_file)) {
  ka_deg_file <- "../220205_Madhu_RNASeq_ReAnalysis_KI-KA-DBonly/output_250818/KA_3ClassV_vs_KA_ShamSaline_DEtable_padj0.05_250502.txt"
}

# Load DEG data
ki_deg <- read.delim(ki_deg_file, sep = "\t", header = TRUE)
ka_deg <- read.delim(ka_deg_file, sep = "\t", header = TRUE)

cat("KI DEGs:", nrow(ki_deg), "\n")
cat("KA DEGs:", nrow(ka_deg), "\n")

# Extract gene lists
ki_deg_genes <- ki_deg$SYMBOL
ka_deg_genes <- ka_deg$SYMBOL

################################################################################
# 6. Load DMG Results
################################################################################

cat("\n=== Loading DMG results ===\n")

# Load annotated DMC data with gene information
ki_dmc_file <- file.path(dmc_dir, "KI_3ClassV_vs_KI_Sham_DECpG-diff25q001.txt")
ka_dmc_file <- file.path(dmc_dir, "KA_3ClassV_vs_KA_ShamSaline_DECpG-diff25q001.txt")

# Check if files exist
if (!file.exists(ki_dmc_file)) {
  stop("KI DMC file not found: ", ki_dmc_file)
}
if (!file.exists(ka_dmc_file)) {
  stop("KA DMC file not found: ", ka_dmc_file)
}

ki_dmc <- read.delim(ki_dmc_file, sep = "\t", header = TRUE)
ka_dmc <- read.delim(ka_dmc_file, sep = "\t", header = TRUE)

cat("KI DMCs:", nrow(ki_dmc), "\n")
cat("KA DMCs:", nrow(ka_dmc), "\n")

# Load gene annotation objects
if (file.exists("gene.obj.rdata")) {
  load("gene.obj.rdata")
} else {
  stop("gene.obj.rdata not found. Please run 02_RRBS_Analysis.R first.")
}

if (file.exists("ensembl.gene.rdata")) {
  load("ensembl.gene.rdata")
} else {
  stop("ensembl.gene.rdata not found. Please run 02_RRBS_Analysis.R first.")
}

################################################################################
# 7. Annotate DMCs with Genes
################################################################################

cat("\n=== Annotating DMCs with genes ===\n")

# Function to annotate CpGs (same as in 02_RRBS_Analysis.R)
annotate_cpgs_with_genes <- function(cpg_data, gene.obj, ensembl.gene) {
  cat("  Annotating", nrow(cpg_data), "DMCpGs with genes...\n")

  cpg_gr <- GRanges(
    seqnames = cpg_data$chr,
    ranges = IRanges(start = cpg_data$start, end = cpg_data$end),
    strand = "*"
  )

  suppressWarnings({
    gene_annot <- genomation::annotateWithGeneParts(cpg_gr, gene.obj)
    tss_annot <- genomation::getAssociationWithTSS(gene_annot)
  })

  annot_df <- as.data.frame(gene_annot@members)
  cpg_data$prom <- as.integer(annot_df$prom)
  cpg_data$exon <- as.integer(annot_df$exon)
  cpg_data$intron <- as.integer(annot_df$intron)

  tss_df <- as.data.frame(tss_annot)
  cpg_data$feature.name <- tss_df$feature.name
  cpg_data$dist.to.feature <- tss_df$dist.to.feature

  gene_id_to_symbol <- setNames(
    ensembl.gene$external_gene_name,
    ensembl.gene$ensembl_gene_id
  )

  cpg_data$external_gene_name <- gene_id_to_symbol[cpg_data$feature.name]
  cpg_data$external_gene_name <- ifelse(
    is.na(cpg_data$external_gene_name),
    cpg_data$feature.name,
    cpg_data$external_gene_name
  )

  cat("  Annotated", sum(!is.na(cpg_data$external_gene_name)), "CpG-gene associations\n")

  return(cpg_data)
}

# Annotate DMCs
ki_dmc_annotated <- annotate_cpgs_with_genes(ki_dmc, gene.obj, ensembl.gene)
ka_dmc_annotated <- annotate_cpgs_with_genes(ka_dmc, gene.obj, ensembl.gene)

# Create DMG lists (unique genes from annotated CpGs)
ki_dmg_genes <- unique(ki_dmc_annotated$external_gene_name[
  !is.na(ki_dmc_annotated$external_gene_name) &
  ki_dmc_annotated$external_gene_name != ""
])

ka_dmg_genes <- unique(ka_dmc_annotated$external_gene_name[
  !is.na(ka_dmc_annotated$external_gene_name) &
  ka_dmc_annotated$external_gene_name != ""
])

cat("KI DMGs:", length(ki_dmg_genes), "\n")
cat("KA DMGs:", length(ka_dmg_genes), "\n")

################################################################################
# 8. Within-Model DEG-DMG Overlap Analysis
################################################################################

cat("\n=== Analyzing within-model DEG-DMG overlaps ===\n")

# KI model overlap
ki_shared_genes <- intersect(ki_deg_genes, ki_dmg_genes)
cat("KI: Genes with both DEG and DMG:", length(ki_shared_genes), "\n")

# KA model overlap
ka_shared_genes <- intersect(ka_deg_genes, ka_dmg_genes)
cat("KA: Genes with both DEG and DMG:", length(ka_shared_genes), "\n")

# Create Venn diagrams
pdf(file.path(output_dir, "01_Venn_DEG_DMG_KI.pdf"), width = 8, height = 8)
venn_ki <- venndetail(list(DEG = ki_deg_genes, DMG = ki_dmg_genes))
plot(venn_ki, type = "venn")
dev.off()

pdf(file.path(output_dir, "01_Venn_DEG_DMG_KA.pdf"), width = 8, height = 8)
venn_ka <- venndetail(list(DEG = ka_deg_genes, DMG = ka_dmg_genes))
plot(venn_ka, type = "venn")
dev.off()

################################################################################
# 9. Extract Methylation and Expression Data for Shared Genes
################################################################################

cat("\n=== Extracting data for shared genes ===\n")

# Function to get gene-level methylation summary
get_gene_meth_summary <- function(cpg_annotated, genes) {
  cpg_subset <- cpg_annotated[cpg_annotated$external_gene_name %in% genes, ]

  gene_summary <- cpg_subset %>%
    group_by(external_gene_name) %>%
    summarise(
      mean_meth_diff = mean(meth.diff, na.rm = TRUE),
      n_cpgs = n(),
      min_qvalue = min(qvalue, na.rm = TRUE),
      .groups = "drop"
    )

  return(gene_summary)
}

# KI model
ki_shared_data <- data.frame(SYMBOL = ki_shared_genes)

# Add expression data
ki_deg_subset <- ki_deg[ki_deg$SYMBOL %in% ki_shared_genes, ]
ki_shared_data <- ki_shared_data %>%
  left_join(
    ki_deg_subset %>% dplyr::select(SYMBOL, log2FoldChange, padj),
    by = "SYMBOL"
  )

# Add methylation data
ki_meth_summary <- get_gene_meth_summary(ki_dmc_annotated, ki_shared_genes)
ki_shared_data <- ki_shared_data %>%
  left_join(
    ki_meth_summary,
    by = c("SYMBOL" = "external_gene_name")
  )

# KA model
ka_shared_data <- data.frame(SYMBOL = ka_shared_genes)

ka_deg_subset <- ka_deg[ka_deg$SYMBOL %in% ka_shared_genes, ]
ka_shared_data <- ka_shared_data %>%
  left_join(
    ka_deg_subset %>% dplyr::select(SYMBOL, log2FoldChange, padj),
    by = "SYMBOL"
  )

ka_meth_summary <- get_gene_meth_summary(ka_dmc_annotated, ka_shared_genes)
ka_shared_data <- ka_shared_data %>%
  left_join(
    ka_meth_summary,
    by = c("SYMBOL" = "external_gene_name")
  )

################################################################################
# 10. Directional Categorization
################################################################################

cat("\n=== Categorizing genes by direction of change ===\n")

# Function to categorize genes
categorize_direction <- function(data) {
  data$category <- ifelse(
    data$mean_meth_diff > 0 & data$log2FoldChange > 0,
    "Hyper-Up",
    ifelse(
      data$mean_meth_diff > 0 & data$log2FoldChange < 0,
      "Hyper-Down",
      ifelse(
        data$mean_meth_diff < 0 & data$log2FoldChange > 0,
        "Hypo-Up",
        "Hypo-Down"
      )
    )
  )
  return(data)
}

ki_shared_data <- categorize_direction(ki_shared_data)
ka_shared_data <- categorize_direction(ka_shared_data)

cat("KI model categories:\n")
print(table(ki_shared_data$category))

cat("\nKA model categories:\n")
print(table(ka_shared_data$category))

# Save shared gene data
write.table(
  ki_shared_data,
  file = file.path(output_dir, "KI_Shared_DEG_DMG.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  ka_shared_data,
  file = file.path(output_dir, "KA_Shared_DEG_DMG.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 11. Scatter Plots of Methylation vs Expression
################################################################################

cat("\n=== Creating scatter plots ===\n")

# KI model scatter plot
p_ki_scatter <- ggplot(ki_shared_data,
                       aes(x = mean_meth_diff, y = log2FoldChange, color = category)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Hyper-Down" = "blue",
    "Hypo-Up" = "blue",
    "Hyper-Up" = "orange",
    "Hypo-Down" = "orange"
  )) +
  labs(
    x = "Mean Methylation Difference (%)",
    y = "Log2 Fold Change (Expression)",
    title = "KI Model: Genes with Both DMG and DEG",
    color = "Category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "02_Scatter_KI_Meth_vs_Expr.pdf"), p_ki_scatter,
       width = 10, height = 8)

# KA model scatter plot
p_ka_scatter <- ggplot(ka_shared_data,
                       aes(x = mean_meth_diff, y = log2FoldChange, color = category)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c(
    "Hyper-Down" = "blue",
    "Hypo-Up" = "blue",
    "Hyper-Up" = "orange",
    "Hypo-Down" = "orange"
  )) +
  labs(
    x = "Mean Methylation Difference (%)",
    y = "Log2 Fold Change (Expression)",
    title = "KA Model: Genes with Both DMG and DEG",
    color = "Category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "02_Scatter_KA_Meth_vs_Expr.pdf"), p_ka_scatter,
       width = 10, height = 8)

################################################################################
# 12. Cross-Model Concordant Genes
################################################################################

cat("\n=== Identifying model-independent genes ===\n")

# Shared DEGs between models
shared_deg <- intersect(ki_deg_genes, ka_deg_genes)
cat("Shared DEGs between KI and KA:", length(shared_deg), "\n")

# Shared DMGs between models
shared_dmg <- intersect(ki_dmg_genes, ka_dmg_genes)
cat("Shared DMGs between KI and KA:", length(shared_dmg), "\n")

# Genes that are both DEG and DMG in BOTH models
model_independent_genes <- intersect(shared_deg, shared_dmg)
cat("Model-independent genes (DEG+DMG in both):", length(model_independent_genes), "\n")

if (length(model_independent_genes) > 0) {
  # Extract data for these genes
  model_indep_data <- data.frame(SYMBOL = model_independent_genes)

  # Add KI data
  model_indep_data <- model_indep_data %>%
    left_join(
      ki_deg %>% dplyr::select(SYMBOL, log2FoldChange, padj),
      by = "SYMBOL",
      suffix = c("", "_KI_expr")
    ) %>%
    dplyr::rename(KI_log2FC = log2FoldChange, KI_expr_padj = padj) %>%
    left_join(
      ki_meth_summary,
      by = c("SYMBOL" = "external_gene_name")
    ) %>%
    dplyr::rename(KI_meth_diff = mean_meth_diff, KI_n_cpgs = n_cpgs, KI_meth_qval = min_qvalue)

  # Add KA data
  model_indep_data <- model_indep_data %>%
    left_join(
      ka_deg %>% dplyr::select(SYMBOL, log2FoldChange, padj),
      by = "SYMBOL"
    ) %>%
    dplyr::rename(KA_log2FC = log2FoldChange, KA_expr_padj = padj) %>%
    left_join(
      ka_meth_summary,
      by = c("SYMBOL" = "external_gene_name")
    ) %>%
    dplyr::rename(KA_meth_diff = mean_meth_diff, KA_n_cpgs = n_cpgs, KA_meth_qval = min_qvalue)

  # Check concordance
  model_indep_data$expr_concordant <- sign(model_indep_data$KI_log2FC) ==
                                       sign(model_indep_data$KA_log2FC)
  model_indep_data$meth_concordant <- sign(model_indep_data$KI_meth_diff) ==
                                       sign(model_indep_data$KA_meth_diff)

  cat("\nExpression concordance:\n")
  print(table(model_indep_data$expr_concordant))

  cat("\nMethylation concordance:\n")
  print(table(model_indep_data$meth_concordant))

  # Save
  write.table(
    model_indep_data,
    file = file.path(output_dir, "Model_Independent_Genes.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  # Create visualization for top genes
  if (nrow(model_indep_data) >= 5) {
    # Sort by average absolute log2FC
    model_indep_data$avg_abs_log2FC <- (abs(model_indep_data$KI_log2FC) +
                                         abs(model_indep_data$KA_log2FC)) / 2
    top_genes <- model_indep_data[order(-model_indep_data$avg_abs_log2FC), ][1:min(10, nrow(model_indep_data)), ]

    cat("\nTop", nrow(top_genes), "model-independent genes:\n")
    print(top_genes[, c("SYMBOL", "KI_log2FC", "KA_log2FC", "KI_meth_diff", "KA_meth_diff")])
  }
}

# Create Venn diagram showing overlap structure
pdf(file.path(output_dir, "03_Venn_ModelIndependent.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# DEGs
library(VennDetail)
venn_deg_both <- venndetail(list(
  KI_DEG = ki_deg_genes,
  KA_DEG = ka_deg_genes
))
plot(venn_deg_both, type = "venn", main = "DEGs")

# DMGs
venn_dmg_both <- venndetail(list(
  KI_DMG = ki_dmg_genes,
  KA_DMG = ka_dmg_genes
))
plot(venn_dmg_both, type = "venn", main = "DMGs")

dev.off()

################################################################################
# 13. CpG-Site Level Visualization for Key Genes
################################################################################

cat("\n=== Creating CpG-site level visualizations ===\n")

# Load beta values for visualization
beta_data <- read.delim(beta_file, sep = "\t", header = TRUE)

# Function to visualize CpG sites for a gene
visualize_gene_cpgs <- function(gene_name, cpg_annotated_ki, cpg_annotated_ka,
                                 beta_data, samples, output_file) {

  # Get CpG sites for this gene
  gene_cpgs_ki <- cpg_annotated_ki[cpg_annotated_ki$external_gene_name == gene_name, ]
  gene_cpgs_ka <- cpg_annotated_ka[cpg_annotated_ka$external_gene_name == gene_name, ]

  if (nrow(gene_cpgs_ki) == 0 && nrow(gene_cpgs_ka) == 0) {
    cat("No CpG sites found for gene:", gene_name, "\n")
    return(NULL)
  }

  # Combine unique CpG sites
  all_cpg_sites <- unique(c(
    paste0(gene_cpgs_ki$chr, ":", gene_cpgs_ki$start),
    paste0(gene_cpgs_ka$chr, ":", gene_cpgs_ka$start)
  ))

  # Extract beta values for these sites
  beta_data$cpg_id <- paste0(beta_data$chr, ":", beta_data$start)
  gene_beta <- beta_data[beta_data$cpg_id %in% all_cpg_sites, ]

  if (nrow(gene_beta) == 0) {
    cat("No beta values found for gene:", gene_name, "\n")
    return(NULL)
  }

  # Prepare data for plotting
  sample_cols <- intersect(samples$sample, colnames(gene_beta))
  plot_data <- gene_beta[, c("chr", "start", "cpg_id", sample_cols)]

  # Reshape to long format
  plot_data_long <- plot_data %>%
    tidyr::pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "beta"
    )

  # Add group information
  plot_data_long <- plot_data_long %>%
    left_join(samples[, c("sample", "group", "model")], by = "sample")

  # Mark significant CpGs
  plot_data_long$significant_ki <- plot_data_long$cpg_id %in%
    paste0(gene_cpgs_ki$chr, ":", gene_cpgs_ki$start)
  plot_data_long$significant_ka <- plot_data_long$cpg_id %in%
    paste0(gene_cpgs_ka$chr, ":", gene_cpgs_ka$start)

  # Create plot
  p <- ggplot(plot_data_long, aes(x = as.factor(start), y = beta, color = group)) +
    geom_boxplot() +
    geom_rect(
      data = plot_data_long[plot_data_long$significant_ki, ],
      aes(xmin = as.numeric(as.factor(start)) - 0.5,
          xmax = as.numeric(as.factor(start)) + 0.5,
          ymin = -Inf, ymax = Inf),
      fill = "pink", alpha = 0.2, color = NA, inherit.aes = FALSE
    ) +
    labs(
      title = paste0("CpG Methylation for Gene: ", gene_name),
      x = "CpG Position",
      y = "Beta Value (Methylation)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  ggsave(output_file, p, width = 12, height = 6)

  return(plot_data_long)
}

# Visualize top model-independent genes (if any)
if (exists("model_indep_data") && nrow(model_indep_data) > 0) {
  top_genes_to_plot <- head(model_indep_data$SYMBOL, 5)

  for (gene in top_genes_to_plot) {
    output_file <- file.path(output_dir, paste0("04_CpG_Sites_", gene, ".pdf"))
    tryCatch({
      visualize_gene_cpgs(
        gene,
        ki_dmc_annotated,
        ka_dmc_annotated,
        beta_data,
        samples,
        output_file
      )
      cat("Created visualization for", gene, "\n")
    }, error = function(e) {
      cat("Error creating plot for", gene, ":", conditionMessage(e), "\n")
    })
  }
}

################################################################################
# 14. Model-Specific Correlation Analysis
################################################################################

cat("\n=== Performing model-specific correlation analysis ===\n")

# Load FPKM data for correlation analysis
# (Using FPKM rather than counts for correlations as they're already normalized)
if (file.exists(fpkm_file)) {
  cat("Loading FPKM data from:", fpkm_file, "\n")
  fpkm_data <- read.delim(fpkm_file, sep = "\t", header = TRUE, row.names = 1)
  fpkm_data <- fpkm_data[, !grepl("Mean", colnames(fpkm_data))]
  cat("FPKM data dimensions:", nrow(fpkm_data), "genes ×", ncol(fpkm_data), "samples\n")
} else {
  stop("FPKM file not found: ", fpkm_file, "\n",
       "FPKM values are needed for correlation analysis with methylation beta values.")
}

# Function to calculate correlations for a model
calculate_correlations <- function(shared_genes, cpg_annotated, fpkm_data,
                                   beta_data, samples_subset) {

  # Prepare gene-CpG pairs
  gene_cpg_pairs <- cpg_annotated %>%
    filter(external_gene_name %in% shared_genes) %>%
    dplyr::select(external_gene_name, chr, start, meth.diff, qvalue) %>%
    distinct()

  if (nrow(gene_cpg_pairs) == 0) {
    return(NULL)
  }

  # Match samples
  fpkm_samples <- intersect(samples_subset$sample, colnames(fpkm_data))
  beta_samples <- intersect(samples_subset$sample, colnames(beta_data))
  common_samples <- intersect(fpkm_samples, beta_samples)

  if (length(common_samples) < 3) {
    cat("Not enough samples for correlation\n")
    return(NULL)
  }

  # Calculate correlations
  correlations <- list()

  for (i in 1:nrow(gene_cpg_pairs)) {
    gene <- gene_cpg_pairs$external_gene_name[i]
    cpg_chr <- gene_cpg_pairs$chr[i]
    cpg_start <- gene_cpg_pairs$start[i]

    # Get expression values
    if (gene %in% rownames(fpkm_data)) {
      expr_vals <- as.numeric(fpkm_data[gene, common_samples])

      # Get methylation values
      beta_cpg <- beta_data[beta_data$chr == cpg_chr & beta_data$start == cpg_start, ]

      if (nrow(beta_cpg) > 0) {
        meth_vals <- as.numeric(beta_cpg[1, common_samples])

        # Remove missing values
        valid_idx <- !is.na(expr_vals) & !is.na(meth_vals)

        if (sum(valid_idx) >= 3) {
          # Calculate Spearman correlation
          cor_test <- cor.test(
            meth_vals[valid_idx],
            expr_vals[valid_idx],
            method = "spearman"
          )

          correlations[[length(correlations) + 1]] <- data.frame(
            gene = gene,
            cpg_chr = cpg_chr,
            cpg_start = cpg_start,
            correlation = cor_test$estimate,
            pvalue = cor_test$p.value,
            n_samples = sum(valid_idx)
          )
        }
      }
    }
  }

  if (length(correlations) > 0) {
    cor_df <- do.call(rbind, correlations)
    cor_df$padj <- p.adjust(cor_df$pvalue, method = "BH")
    return(cor_df)
  } else {
    return(NULL)
  }
}

# Prepare beta data for correlation
beta_data_cor <- beta_data
if (!"cpg.site" %in% colnames(beta_data_cor)) {
  beta_data_cor$cpg.site <- paste0(beta_data_cor$chr, "_", beta_data_cor$start)
}

# KI model correlations
ki_samples <- samples[samples$model == "KI", ]
ki_correlations <- calculate_correlations(
  ki_shared_genes,
  ki_dmc_annotated,
  fpkm_data,
  beta_data_cor,
  ki_samples
)

if (!is.null(ki_correlations)) {
  cat("KI model: Calculated", nrow(ki_correlations), "correlations\n")
  cat("  Significant (padj < 0.05):", sum(ki_correlations$padj < 0.05, na.rm = TRUE), "\n")

  write.table(
    ki_correlations,
    file = file.path(output_dir, "KI_CpG_Gene_Correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# KA model correlations
ka_samples <- samples[samples$model == "KA", ]
ka_correlations <- calculate_correlations(
  ka_shared_genes,
  ka_dmc_annotated,
  fpkm_data,
  beta_data_cor,
  ka_samples
)

if (!is.null(ka_correlations)) {
  cat("KA model: Calculated", nrow(ka_correlations), "correlations\n")
  cat("  Significant (padj < 0.05):", sum(ka_correlations$padj < 0.05, na.rm = TRUE), "\n")

  write.table(
    ka_correlations,
    file = file.path(output_dir, "KA_CpG_Gene_Correlations.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# Create correlation distribution plots
if (!is.null(ki_correlations) && !is.null(ka_correlations)) {
  combined_cor <- rbind(
    data.frame(model = "KI", correlation = ki_correlations$correlation),
    data.frame(model = "KA", correlation = ka_correlations$correlation)
  )

  p_cor_dist <- ggplot(combined_cor, aes(x = correlation, fill = model)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    labs(
      title = "Distribution of CpG-Gene Correlations by Model",
      x = "Spearman Correlation",
      y = "Count"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  ggsave(file.path(output_dir, "05_Correlation_Distribution.pdf"), p_cor_dist,
         width = 10, height = 6)
}

################################################################################
# 15. Permutation Testing for DEG-DMG Overlaps
################################################################################

cat("\n=== Performing permutation testing ===\n")

permutation_test_overlap <- function(set1, set2, universe, n_perm = 10000) {
  observed_overlap <- length(intersect(set1, set2))
  size1 <- length(set1)
  size2 <- length(set2)

  null_overlaps <- replicate(n_perm, {
    random_set1 <- sample(universe, size1)
    random_set2 <- sample(universe, size2)
    length(intersect(random_set1, random_set2))
  })

  mean_null <- mean(null_overlaps)
  sd_null <- sd(null_overlaps)
  z_score <- (observed_overlap - mean_null) / sd_null
  p_value <- sum(null_overlaps >= observed_overlap) / n_perm

  return(list(
    observed = observed_overlap,
    expected = mean_null,
    sd = sd_null,
    z_score = z_score,
    p_value = p_value
  ))
}

# Define universe (all genes with both expression and methylation data)
universe_expr <- rownames(fpkm_data)
universe_meth <- unique(unlist(c(
  ki_dmc_annotated$external_gene_name[!is.na(ki_dmc_annotated$external_gene_name)],
  ka_dmc_annotated$external_gene_name[!is.na(ka_dmc_annotated$external_gene_name)]
)))
universe <- intersect(universe_expr, universe_meth)

cat("Universe size (genes with both expr and meth data):", length(universe), "\n")

# Test KI DEG-DMG overlap
perm_ki <- permutation_test_overlap(ki_deg_genes, ki_dmg_genes, universe, n_perm = 10000)

cat("\nKI DEG-DMG Overlap Permutation Test:\n")
cat("  Observed:", perm_ki$observed, "\n")
cat("  Expected:", round(perm_ki$expected, 2), "\n")
cat("  Z-score:", round(perm_ki$z_score, 2), "\n")
cat("  P-value:", perm_ki$p_value, "\n")

# Test KA DEG-DMG overlap
perm_ka <- permutation_test_overlap(ka_deg_genes, ka_dmg_genes, universe, n_perm = 10000)

cat("\nKA DEG-DMG Overlap Permutation Test:\n")
cat("  Observed:", perm_ka$observed, "\n")
cat("  Expected:", round(perm_ka$expected, 2), "\n")
cat("  Z-score:", round(perm_ka$z_score, 2), "\n")
cat("  P-value:", perm_ka$p_value, "\n")

# Save results
perm_results <- data.frame(
  Model = c("KI", "KA"),
  Observed = c(perm_ki$observed, perm_ka$observed),
  Expected = c(perm_ki$expected, perm_ka$expected),
  SD = c(perm_ki$sd, perm_ka$sd),
  Z_score = c(perm_ki$z_score, perm_ka$z_score),
  P_value = c(perm_ki$p_value, perm_ka$p_value)
)

write.table(
  perm_results,
  file = file.path(output_dir, "Permutation_Test_DEG_DMG_Overlap.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# 16. Save Workspace
################################################################################

cat("\n=== Saving workspace ===\n")

save.image(file = file.path(output_dir, "Integrative_Analysis_Workspace.RData"))

cat("\n=== Integrative analysis complete ===\n")
cat("Output directory:", output_dir, "\n")

################################################################################
# End of Script
################################################################################
